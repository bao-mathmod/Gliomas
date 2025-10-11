#!/usr/bin/env Rscript

#Usage: 
# Rscript /mnt/18T/chibao/gliomas/code/04.upstream/4.2_scATAC/samplesheet.r \
#   --roots "/mnt/18T/chibao/gliomas/data/output_cell/snATAC,/mnt/18T/chibao/gliomas/data/output_cell/ATAC_multiome" \
#   --outdir "/mnt/18T/chibao/gliomas/atlas_atac" \
#   --genome-map "PRJNA578617=hg19"

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(fs)
})

# ---------------------------
# CLI
# ---------------------------
opt_list <- list(
  make_option("--roots", type="character", default=NULL,
              help="Comma-separated list of roots to scan (recursively)."),
  make_option("--outdir", type="character", default="/mnt/18T/chibao/gliomas/atlas_atac",
              help="Output directory for manifests/ and logs/. [default %default]"),
  make_option("--genome-map", type="character", default="PRJNA578617=hg19",
              help="Comma-separated project=genome overrides (e.g., 'PRJNA578617=hg19,PRJNAxxxxxx=hg19'). Others default to hg38.")
)

opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$roots) || opt$roots == "") {
  opt$roots <- paste(
    "/mnt/18T/chibao/gliomas/data/output_cell/snATAC",
    "/mnt/18T/chibao/gliomas/data/output_cell/ATAC_multiome",
    sep = ","
  )
}

roots <- strsplit(opt$roots, "\\s*,\\s*")[[1]] %>% unique()
outdir <- opt$outdir
genome_map_in <- opt$`genome-map`

# ---------------------------
# Setup outputs (minimal files)
# ---------------------------
dir_create(file.path(outdir, "manifests"))
dir_create(file.path(outdir, "logs"))

log_main  <- file.path(outdir, "logs", "main.log")
log_error <- file.path(outdir, "logs", "error.log")
manifest_path <- file.path(outdir, "manifests", "atac_manifest.tsv")

log_msg <- function(..., .error=FALSE){
  msg <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", paste(..., collapse=" "))
  cat(msg, "\n", file = if (.error) log_error else log_main, append = TRUE)
}

# ---------------------------
# Parse genome overrides
# ---------------------------
genome_overrides <- list()
if (!is.null(genome_map_in) && nchar(genome_map_in) > 0) {
  kv <- strsplit(genome_map_in, "\\s*,\\s*")[[1]]
  for (p in kv) {
    if (grepl("=", p)) {
      k <- sub("=.*$", "", p)
      v <- sub("^.*?=", "", p)
      genome_overrides[[k]] <- v
    }
  }
}

# Vectorized inference
infer_genome <- function(project_ids){
  vapply(project_ids, function(pid){
    if (!is.null(genome_overrides[[pid]])) genome_overrides[[pid]] else "hg38"
  }, FUN.VALUE = character(1))
}

# ---------------------------
# Helpers
# ---------------------------
guess_project <- function(path){
  parts <- strsplit(path, "/")[[1]]
  prj <- parts[grepl("^PRJNA\\d+$", parts)]
  if (length(prj) > 0) return(prj[[length(prj)]])
  return(basename(dirname(path)))
}

# Normalize to a per-sample "root" (filename stem without format-specific suffix)
sample_root_from_fname <- function(fname){
  root <- fname
  # 10x h5
  root <- gsub("_filtered_peak_bc_matrix\\.h5$", "", root)
  root <- gsub("\\.raw_peak_bc_matrix\\.h5$", "", root)
  # 10x singlecell
  root <- gsub("\\.singlecell\\.csv\\.gz$", "", root)
  # triplets
  root <- gsub("_matrix\\.mtx\\.gz$", "", root)
  root <- gsub("_peaks\\.bed\\.gz$", "", root)
  root <- gsub("_barcodes\\.tsv\\.gz$", "", root)
  root <- gsub("\\.features\\.tsv\\.gz$", "", root)
  # fragments
  root <- gsub("_atac_fragments\\.tsv\\.gz$", "", root)
  # generic strip common double extensions
  root <- gsub("\\.tsv\\.gz$", "", root)
  root <- gsub("\\.gz$", "", root)
  return(root)
}

# ---------------------------
# Scan one root for files
# ---------------------------
scan_root <- function(root){
  if (!dir_exists(root)) {
    log_msg("Root does not exist, skipping:", root, .error = TRUE)
    return(tibble())
  }
  all <- dir_ls(root, recurse = TRUE, type = "file", fail = FALSE)
  if (!length(all)) return(tibble())

  fi  <- file_info(all)

  tibble(path = all,
         size = fi$size,
         mtime = as.double(fi$modification_time)) %>%
    mutate(
      fname = path_file(path),
      dir   = path_dir(path),
      project_id = map_chr(path, guess_project),
      sample_root = sample_root_from_fname(fname),
      is_h5_filtered = grepl("filtered_peak_bc_matrix\\.h5$", fname),
      is_h5_raw      = grepl("\\.raw_peak_bc_matrix\\.h5$", fname),
      is_mtx         = grepl("matrix\\.mtx\\.gz$", fname),
      is_peaks       = grepl("peaks\\.bed\\.gz$", fname),
      is_barcodes    = grepl("barcodes\\.tsv\\.gz$", fname),
      is_features    = grepl("features\\.tsv\\.gz$", fname),
      is_fragments   = grepl("atac_fragments\\.tsv\\.gz$", fname),
      is_tbi         = grepl("\\.tbi\\.gz$", fname),
      is_singlecell  = grepl("singlecell\\.csv\\.gz$", fname)
    )
}

log_msg("Phase 0 | Scanning roots:", paste(roots, collapse=", "))
df <- map_dfr(roots, scan_root)

if (nrow(df) == 0) {
  log_msg("No files found under roots. Nothing to do.", .error = TRUE)
  quit(status = 1)
}

# ---------------------------
# Assemble candidates by sample_root (avoid dir-based joins)
# ---------------------------

# 10x_h5
h5_tbl <- df %>%
  filter(is_h5_filtered | is_h5_raw) %>%
  transmute(
    project_id, dir, sample_root,
    h5 = path,
    type = "10x_h5"
  ) %>%
  distinct()

# Singlecell: choose ONE per (project_id, sample_root) — largest file
singlecell_tbl <- df %>%
  filter(is_singlecell) %>%
  group_by(project_id, sample_root) %>%
  slice_max(order_by = size, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(project_id, sample_root, singlecell = path)

h5_tbl <- h5_tbl %>%
  left_join(singlecell_tbl, by = c("project_id","sample_root"))

# Triplets: require (mtx, peaks, barcodes) with SAME (project_id, sample_root)
mtx_tbl <- df %>%
  filter(is_mtx | is_peaks | is_barcodes | is_features | is_singlecell) %>%
  group_by(project_id, sample_root) %>%
  summarise(
    mtx       = path[is_mtx][1] %||% NA_character_,
    peaks     = path[is_peaks][1] %||% NA_character_,
    barcodes  = path[is_barcodes][1] %||% NA_character_,
    features  = path[is_features][1] %||% NA_character_,
    singlecell= path[is_singlecell][1] %||% NA_character_,
    .groups = "drop"
  ) %>%
  filter(!is.na(mtx) & !is.na(peaks) & !is.na(barcodes)) %>%
  mutate(type = "mtx_triplet")

# Fragments-only
frag_tbl <- df %>%
  filter(is_fragments) %>%
  group_by(project_id, sample_root) %>%
  summarise(
    fragments = path[1],
    has_index = any(is_tbi),
    .groups = "drop"
  ) %>%
  mutate(type = "fragments")

# ---------------------------
# Build normalized rows; use project_id+sample_root for sample_id
# Priority: 10x_h5 > mtx_triplet > fragments
# ---------------------------
norm_h5 <- h5_tbl %>%
  transmute(
    project_id,
    sample_root,
    sample_id = paste0(project_id, "_", sample_root),
    type,
    genome = infer_genome(project_id),
    h5 = h5,
    mtx = NA_character_,
    peaks = NA_character_,
    barcodes = NA_character_,
    fragments = NA_character_,
    features = NA_character_,
    singlecell = coalesce(singlecell, NA_character_),
    notes = NA_character_
  )

norm_mtx <- mtx_tbl %>%
  transmute(
    project_id,
    sample_root,
    sample_id = paste0(project_id, "_", sample_root),
    type,
    genome = infer_genome(project_id),
    h5 = NA_character_,
    mtx = mtx,
    peaks = peaks,
    barcodes = barcodes,
    fragments = NA_character_,
    features = coalesce(features, NA_character_),
    singlecell = coalesce(singlecell, NA_character_),
    notes = NA_character_
  )

norm_frag <- frag_tbl %>%
  transmute(
    project_id,
    sample_root,
    sample_id = paste0(project_id, "_", sample_root),
    type,
    genome = infer_genome(project_id),
    h5 = NA_character_,
    mtx = NA_character_,
    peaks = NA_character_,
    barcodes = NA_character_,
    fragments = fragments,
    features = NA_character_,
    singlecell = NA_character_,
    notes = ifelse(has_index, "fragments_index=present", "fragments_index=absent")
  )

combined <- bind_rows(norm_h5, norm_mtx, norm_frag) %>%
  arrange(project_id,
          factor(type, levels=c("10x_h5","mtx_triplet","fragments")),
          sample_id)

manifest <- combined %>%
  group_by(sample_id) %>%
  slice(1L) %>%                     # enforce priority order
  ungroup() %>%
  select(-sample_root) %>%           # no longer needed
  arrange(project_id, sample_id)

# ---------------------------
# Validation
# ---------------------------
validate_row <- function(row){
  errs <- c()
  if (row$type == "10x_h5") {
    if (is.na(row$h5) || !file_exists(row$h5)) errs <- c(errs, "missing_h5")
  } else if (row$type == "mtx_triplet") {
    if (is.na(row$mtx) || !file_exists(row$mtx)) errs <- c(errs, "missing_mtx")
    if (is.na(row$peaks) || !file_exists(row$peaks)) errs <- c(errs, "missing_peaks")
    if (is.na(row$barcodes) || !file_exists(row$barcodes)) errs <- c(errs, "missing_barcodes")
  } else if (row$type == "fragments") {
    if (is.na(row$fragments) || !file_exists(row$fragments)) errs <- c(errs, "missing_fragments")
  }
  errs
}

problems <- list()
if (nrow(manifest) > 0){
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i,]
    errs <- validate_row(row)
    if (length(errs)) {
      problems[[row$sample_id]] <- errs
      log_msg("VALIDATION", row$sample_id, "→", paste(errs, collapse=","), .error = TRUE)
    }
  }
}

if (length(problems)){
  log_msg("Validation found", length(problems), "sample(s) with issues (see error.log). Proceeding with manifest.")
}

# ---------------------------
# Write manifest
# ---------------------------
write_tsv(manifest, manifest_path)

# ---------------------------
# Summaries
# ---------------------------
counts <- manifest %>% count(type, name = "n")
total  <- nrow(manifest)

log_msg("Phase 0 | Detected samples total:", total)
if (nrow(counts)) {
  for (j in seq_len(nrow(counts))) {
    log_msg("Phase 0 |", counts$type[j], "=", counts$n[j])
  }
}
gcounts <- manifest %>% count(genome, name = "n")
for (j in seq_len(nrow(gcounts))) {
  log_msg("Phase 0 | genome", gcounts$genome[j], "=", gcounts$n[j])
}
log_msg("Phase 0 | Manifest written to:", manifest_path)

