#!/usr/bin/env Rscript
# Usage
# Rscript /mnt/18T/chibao/gliomas/code/04.upstream/4.2_scATAC/01_chromatin_assay.r --manifest "/mnt/18T/chibao/gliomas/data/upstream/atlas_atac/manifests/atac_manifest.tsv" --outdir "/mnt/18T/chibao/gliomas/data/upstream/atlas_atac" --overwrite FALSE

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(fs)
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

# ---------------------------
# CLI
# ---------------------------
opt_list <- list(
  make_option("--manifest", type="character", help="Path to Phase-0 manifest TSV."),
  make_option("--outdir",   type="character",
              default="/mnt/18T/chibao/gliomas/data/upstream/atlas_atac",
              help="Base dir containing manifests/, logs/, reports/, rds/. [default %default]"),
  make_option("--overwrite", type="logical", default=FALSE,
              help="Overwrite existing per-sample RDS if present. [default %default]"),
  make_option("--tile-size", type="integer", default=5000,
              help="Tile size (bp) when counting fragments-only samples. [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$manifest) || !file.exists(opt$manifest)) {
  stop("Manifest is required and must exist. Use --manifest path/to/atac_manifest.tsv")
}

# ---------------------------
# Paths & minimal outputs
# ---------------------------
base_dir <- opt$outdir
dir_create(file.path(base_dir, "rds", "raw"))
dir_create(file.path(base_dir, "logs"))
dir_create(file.path(base_dir, "reports"))

log_main    <- file.path(base_dir, "logs", "main.log")
log_error   <- file.path(base_dir, "logs", "error.log")
summary_tsv <- file.path(base_dir, "reports", "atac_import_qc_summary.tsv")

log_msg <- function(..., .error=FALSE){
  msg <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", paste(..., collapse=" "))
  cat(msg, "\n", file = if (.error) log_error else log_main, append = TRUE)
}

append_summary <- function(df_row){
  create_header <- !file.exists(summary_tsv)
  suppressWarnings(
    write_tsv(as.data.frame(df_row), summary_tsv, append = TRUE, col_names = create_header)
  )
}

`%||%` <- function(a,b) if (is.null(a) || is.na(a)) b else a

# ---------------------------
# Helpers
# ---------------------------

# Read optional singlecell.csv.gz and standardize a 'barcode' column
read_singlecell <- function(p){
  if (is.na(p) || !file.exists(p)) return(NULL)
  sc <- tryCatch({
    readr::read_csv(p, show_col_types = FALSE)
  }, error = function(e) NULL)
  if (is.null(sc)) return(NULL)

  nml <- tolower(names(sc))
  if (!"barcode" %in% nml) {
    bc_idx <- which(nml %in% c("barcode","barcodes","cell","cell_id","cellid"))
    if (length(bc_idx)) {
      names(sc)[bc_idx[1]] <- "barcode"
    } else {
      # default to first column
      names(sc)[1] <- "barcode"
    }
  } else {
    names(sc)[nml == "barcode"] <- "barcode"
  }
  sc
}

# Robustly read 10x h5 and return a dgCMatrix of ATAC peaks
read_10x_h5_peaks <- function(h5_path){
  mm <- tryCatch(Signac::Read10X_h5(h5_path), error = function(e) e)
  if (inherits(mm, "error")) {
    mm <- tryCatch(Seurat::Read10X_h5(h5_path), error = function(e) e)
  }
  if (inherits(mm, "error")) stop("Read10X_h5 failed: ", conditionMessage(mm))

  # If a list, pick the ATAC/Peaks matrix (or largest nnz as fallback)
  if (is.list(mm)) {
    nms <- names(mm) %||% rep("", length(mm))
    cand <- which(grepl("peaks|atac", nms, ignore.case = TRUE))
    if (length(cand) == 0L) {
      nnz <- vapply(mm, function(x) if (inherits(x, "dgCMatrix")) length(x@x) else 0L, 0L)
      if (all(nnz == 0L)) {
        stop("Read10X_h5 returned a list but no dgCMatrix-like element found.")
      }
      cand <- which.max(nnz)
    } else {
      cand <- cand[1]
    }
    mm <- mm[[cand]]
  }

  if (!inherits(mm, "dgCMatrix")) {
    stop("Expected dgCMatrix; got: ", paste(class(mm), collapse="/"))
  }
  if (is.null(rownames(mm))) {
    stop("Peaks matrix has no rownames; cannot build GRanges.")
  }

  # Normalize rownames if needed: "chr1_100_200" -> "chr1:100-200"
  if (all(!grepl(":", rownames(mm))) && all(grepl("^[^_]+_[0-9]+_[0-9]+$", rownames(mm)))) {
    rn <- gsub("^([^_]+)_([0-9]+)_([0-9]+)$", "\\1:\\2-\\3", rownames(mm))
    rownames(mm) <- rn
  }

  # Final sanity: require "chr:start-end"
  if (!all(grepl("^[^:]+:[0-9]+-[0-9]+$", rownames(mm)))) {
    stop("Peak rownames are not in 'chr:start-end' format; first few: ",
         paste(utils::head(rownames(mm), 3), collapse=", "))
  }

  mm
}

# Read triplet (mtx+peaks+barcodes) robustly
read_triplet_counts <- function(mtx, peaks, barcodes){
  mm  <- Matrix::readMM(mtx)
  bc  <- readr::read_tsv(barcodes, col_names = FALSE, show_col_types = FALSE)[[1]]
  pk  <- readr::read_tsv(peaks,    col_names = FALSE, show_col_types = FALSE)

  if (ncol(pk) >= 3) {
    feat <- paste0(pk[[1]], ":", pk[[2]], "-", pk[[3]])
    gr   <- GRanges(seqnames = pk[[1]], ranges = IRanges(start = pk[[2]] + 1, end = pk[[3]]))
  } else {
    feat <- pk[[1]]
    parts <- do.call(rbind, strsplit(feat, "[:-]"))
    gr <- GRanges(parts[,1], IRanges(as.integer(parts[,2]), as.integer(parts[,3])))
  }

  mm <- as(mm, "dgCMatrix")
  if (nrow(mm) != length(feat)) stop("Rows in mtx and peaks mismatch.")
  if (ncol(mm) != length(bc))   stop("Cols in mtx and barcodes mismatch.")

  rownames(mm) <- feat
  colnames(mm) <- bc
  list(counts = mm, gr = gr)
}

# Try to infer seqinfo (used for tiles/fragments-only)
seqinfo_for_genome <- function(genome, fragments_path = NULL){
  try_ensdb <- function(){
    if (genome == "hg38" && requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
      edb <- getNamespace("EnsDb.Hsapiens.v86")$EnsDb.Hsapiens.v86
      si <- seqinfo(GenomicFeatures::genes(edb))
      return(keepStandardChromosomes(si, pruning.mode = "coarse"))
    }
    if (genome == "hg19" && requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
      si <- seqinfo(GenomicFeatures::genes(txdb))
      return(keepStandardChromosomes(si, pruning.mode = "coarse"))
    }
    NULL
  }
  si <- try_ensdb()
  if (!is.null(si)) return(si)

  if (!is.null(fragments_path) && file.exists(fragments_path)) {
    log_msg("Seqinfo fallback via fragments for:", fragments_path)
    con <- gzfile(fragments_path, open="rt")
    on.exit(close(con), add = TRUE)
    nmax <- 5e5
    chr <- c(); end <- c(); i <- 0L
    repeat {
      ln <- readLines(con, n = 50000L)
      if (!length(ln)) break
      i <- i + length(ln)
      sp <- strsplit(ln, "\t", fixed = TRUE)
      chr <- c(chr, vapply(sp, `[[`, "", 1L))
      end <- c(end, as.integer(vapply(sp, `[[`, "", 3L)))
      if (i >= nmax) break
    }
    df <- tibble(chr = chr, end = end) %>%
      group_by(chr) %>% summarise(max_end = max(end, na.rm = TRUE), .groups="drop")
    df <- df %>% filter(chr %in% paste0("chr", c(1:22,"X","Y")))
    gr <- GRanges(df$chr, IRanges(1, df$max_end))
    return(seqinfo(gr))
  }
  stop("Unable to determine seqinfo for genome=", genome, ". Install EnsDb/TxDb or provide fragments.")
}

# Build counts from fragments into bins (placeholder for QC; real recount happens later)
build_counts_from_fragments <- function(fragments, genome, tile_size = 5000){
  if (!file.exists(fragments)) stop("fragments not found: ", fragments)
  si <- seqinfo_for_genome(genome, fragments)
  si_std <- keepStandardChromosomes(si, pruning.mode = "coarse")
  tiles <- unlist(tileGenome(seqlengths = seqlengths(si_std), tilewidth = tile_size, cut.last.tile.in.chrom = TRUE))
  tiles <- trim(tiles)
  mat <- FeatureMatrix(fragment.path = fragments, features = tiles, cells = NULL)
  rownames(mat) <- paste0(seqnames(tiles), ":", start(tiles), "-", end(tiles))
  list(counts = mat, gr = tiles)
}

# Attach optional singlecell metadata by barcode
attach_singlecell_meta <- function(obj, singlecell_path){
  sc <- read_singlecell(singlecell_path)
  if (is.null(sc) || !"barcode" %in% names(sc)) return(obj)
  common <- intersect(Cells(obj), sc$barcode)
  if (!length(common)) return(obj)
  sc_sub <- sc[match(common, sc$barcode), , drop = FALSE]
  meta <- sc_sub
  rownames(meta) <- meta$barcode
  meta$barcode <- NULL
  AddMetaData(obj, metadata = meta)
}

# RDS path
rds_path_for <- function(base_dir, project_id, sample_id){
  dir_create(file.path(base_dir, "rds", "raw", project_id))
  file.path(base_dir, "rds", "raw", project_id, paste0(sample_id, ".rds"))
}

# ---------------------------
# Load manifest
# ---------------------------
man <- suppressWarnings(readr::read_tsv(opt$manifest, show_col_types = FALSE)) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), NA, .)))

log_msg("Phase 1 | Starting import for", nrow(man), "samples from manifest:", opt$manifest)

# ---------------------------
# Per-sample loop
# ---------------------------
for (i in seq_len(nrow(man))) {
  row <- man[i, ]

  project_id <- row$project_id
  sample_id  <- row$sample_id
  type       <- row$type
  genome     <- row$genome
  h5         <- row$h5
  mtx        <- row$mtx
  peaks      <- row$peaks
  barcodes   <- row$barcodes
  fragments  <- row$fragments
  features   <- row$features
  singlecell <- row$singlecell

  out_rds <- rds_path_for(base_dir, project_id, sample_id)

  if (file.exists(out_rds) && !isTRUE(opt$overwrite)) {
    log_msg("Phase 1 | SKIP existing:", sample_id, "->", out_rds)
    next
  }

  log_msg("Phase 1 | Importing:", sample_id, "| type:", type, "| genome:", genome)

  # --- Build ChromatinAssay by type
  obj <- NULL
  err <- NULL

  if (type == "10x_h5") { 
    # Khởi tạo mat và gr là NULL trước khi thử
    mat <- NULL
    gr <- NULL
    
    ok <- try({
      mat <- read_10x_h5_peaks(h5)
      
      # Chỉ tiếp tục nếu mat được đọc thành công và có dữ liệu
      if(is.null(mat) || ncol(mat) == 0) {
        stop("Read10X_h5: Matrix is NULL or empty after reading.")
      }
      
      gr  <- StringToGRanges(rownames(mat), sep = c(":", "-"))
      
    }, silent = TRUE)

    # Nếu có lỗi (hoặc mat là NULL), gán err
    if (inherits(ok, "try-error")) {
      err <- conditionMessage(attr(ok, "condition"))
    } else if (is.null(mat) || ncol(mat) == 0) {
       err <- "Matrix is NULL or empty after reading 10x h5."
    }

    # Nếu không có lỗi, tiến hành tạo Assay & Object trong khối try mới
    if (is.null(err)) {
        ok_obj <- try({
            assay <- CreateChromatinAssay(
              counts = mat,
              ranges = gr,
              genome = genome,
              sep    = c(":","-")
            )
            # Dòng lỗi xảy ra trước đó là ở đây:
            obj <- CreateSeuratObject(assays = list(ATAC = assay), project = sample_id)
            obj <- attach_singlecell_meta(obj, singlecell)

            log_msg("Phase 1 | 10x_h5 parsed:", sample_id,
                    "rows=", nrow(mat), "cols=", ncol(mat))
        }, silent = TRUE)

        # Nếu có lỗi khi tạo object/assay, gán lỗi vào err
        if (inherits(ok_obj, "try-error")) {
            err <- conditionMessage(attr(ok_obj, "condition"))
        }
    }
  } else if (type == "mtx_triplet") {
    tp <- NULL # Khởi tạo an toàn
    ok <- try({
      tp <- read_triplet_counts(mtx, peaks, barcodes)
      if(is.null(tp) || ncol(tp$counts) == 0) {
        stop("read_triplet_counts failed or returned empty matrix.")
      }
    }, silent = TRUE)

    if (inherits(ok, "try-error")) {
      err <- conditionMessage(attr(ok, "condition"))
    } else if (is.null(tp) || ncol(tp$counts) == 0) {
       err <- "Matrix is NULL or empty after reading mtx triplet."
    }

    if (is.null(err)) {
        ok_obj <- try({
            assay <- CreateChromatinAssay(
                counts = tp$counts,
                ranges = tp$gr,
                genome = genome,
                sep    = c(":","-")
            )
            obj <- CreateSeuratObject(assays = list(ATAC = assay), project = sample_id)
            obj <- attach_singlecell_meta(obj, singlecell)
            log_msg("Phase 1 | mtx_triplet parsed:", sample_id,
                    "rows=", nrow(tp$counts), "cols=", ncol(tp$counts))
        }, silent = TRUE)

        if (inherits(ok_obj, "try-error")) {
            err <- conditionMessage(attr(ok_obj, "condition"))
        }
    }
  } else if (type == "fragments") {
    tc <- NULL # Khởi tạo an toàn
    ok <- try({
      tc <- build_counts_from_fragments(fragments, genome, tile_size = opt$`tile-size`)
      if(is.null(tc) || ncol(tc$counts) == 0) {
        stop("build_counts_from_fragments failed or returned empty matrix.")
      }
    }, silent = TRUE)

    if (inherits(ok, "try-error")) {
      err <- conditionMessage(attr(ok, "condition"))
    } else if (is.null(tc) || ncol(tc$counts) == 0) {
       err <- "Matrix is NULL or empty after tiling fragments."
    }
    
    if (is.null(err)) {
        ok_obj <- try({
            assay <- CreateChromatinAssay(
                counts    = tc$counts,
                ranges    = tc$gr,
                genome    = genome,
                sep       = c(":","-"),
                fragments = fragments
            )
            obj <- CreateSeuratObject(assays = list(ATAC = assay), project = sample_id)
            obj <- attach_singlecell_meta(obj, singlecell)
            log_msg("Phase 1 | fragments tiled:", sample_id,
                    "rows=", nrow(tc$counts), "cols=", ncol(tc$counts),
                    "tile_size=", opt$`tile-size`)
        }, silent = TRUE)
        
        if (inherits(ok_obj, "try-error")) {
            err <- conditionMessage(attr(ok_obj, "condition"))
        }
    }
  } else {
    err <- paste0("Unknown type: ", type)
  }

  if (!is.null(err)) {
    log_msg("Phase 1 | ERROR:", sample_id, "->", err, .error = TRUE)
    next
  }

  # Enforce UCSC style & standard chromosomes when possible
  try({
    seqlevelsStyle(obj[["ATAC"]]@ranges) <- "UCSC"
    obj[["ATAC"]]@ranges <- keepStandardChromosomes(obj[["ATAC"]]@ranges, pruning.mode = "coarse")
  }, silent = TRUE)

  # Basic import stats & save
  n_cells <- ncol(obj)
  n_feat  <- nrow(obj[["ATAC"]])

  saveRDS(obj, out_rds)
  log_msg("Phase 1 | SAVED:", sample_id, "cells=", n_cells, "peaks/bins=", n_feat, "->", out_rds)

  append_summary(tibble::tibble(
    phase = "phase1_import",
    project_id = project_id,
    sample_id  = sample_id,
    type       = type,
    genome     = genome,
    path_h5    = h5    %||% NA_character_,
    path_mtx   = mtx   %||% NA_character_,
    path_peaks = peaks %||% NA_character_,
    path_barc  = barcodes %||% NA_character_,
    path_frags = fragments %||% NA_character_,
    path_singlecell = singlecell %||% NA_character_,
    n_cells    = n_cells,
    n_features = n_feat,
    note       = ifelse(type=="fragments",
                        paste0("tile_size=", opt$`tile-size`),
                        NA_character_)
  ))
}

log_msg("Phase 1 | Completed import for", nrow(man), "samples.")

