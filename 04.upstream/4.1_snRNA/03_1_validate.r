#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(dplyr); library(ggplot2); library(stringr)
})

# Paths (adjust if needed)
before_dir <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds"        # original RDS
after_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds_hgnc_2"   # remapped RDS
audit_dir  <- file.path(dirname(after_dir), "hgnc_audit"); dir.create(audit_dir, showWarnings=FALSE, recursive=TRUE)
audit_tsv  <- file.path(audit_dir, "hgnc_audit_summary.tsv")

# markers to sanity-check (edit/extend as you like)
marker_symbols <- c(
  # mitochondria / ribo
  "MT-ND1","MT-CO1","RPLP0","RPS18",
  # immune / endothelium
  "CD3E","MS4A1","PECAM1","AIF1",
  # fibro/tumor-ish
  "COL1A1","FBLN1","MKI67"
)

# helper to fetch counts matrix from a Seurat object (RNA assay)
get_counts <- function(obj){
  # Hàm này chỉ dùng cú pháp slot cũ, tương thích với Seurat v4
  DefaultAssay(obj) <- "RNA"
  # Tự động lấy slot="counts". Seurat cũ không dùng layer, nên không cần LayerNames
  return(LayerData(obj, assay="RNA", layer="counts"))
}

# initialize audit TSV
if (!file.exists(audit_tsv)) {
  write.table(data.frame(
    sample_id=character(),
    cells_before=integer(), cells_after=integer(),
    genes_before=integer(), genes_after=integer(),
    mapped_n=integer(), unmapped_n=integer(),
    collapsed_symbol_n=integer(), max_collapse_size=integer(),
    counts_r2=double(), counts_max_abs_diff=double(),
    genes_r2=double(), genes_max_abs_diff=double(),
    missing_markers=character(),
    stringsAsFactors=FALSE
  ), file=audit_tsv, sep="\t", quote=FALSE, row.names=FALSE)
}

# find matching pairs by filename
after_files <- list.files(after_dir, pattern="\\.rds$", full.names=TRUE)
for (af in after_files) {
  bn <- basename(af)
  bf <- file.path(before_dir, bn)
  if (!file.exists(bf)) {
    message("Skip (no BEFORE match): ", bn)
    next
  }

  # load
  obj_before <- readRDS(bf)
  obj_after  <- readRDS(af)

  # sample id
  sid <- if (!is.null(obj_after$sample_id)) unique(as.character(obj_after$sample_id))[1] else tools::file_path_sans_ext(bn)

  # counts matrices
  M_before <- get_counts(obj_before)
  M_after  <- get_counts(obj_after)

  # align cells (some pipelines rename with -sample suffix; handle that)
  cells_b <- colnames(M_before); cells_a <- colnames(M_after)
  common_cells <- intersect(cells_b, cells_a)
  if (length(common_cells) == 0) {
    # try stripping sample suffix like "-<sid>"
    strip_suf <- function(x) sub("-.*$", "", x)
    cb2 <- strip_suf(cells_b); ca2 <- strip_suf(cells_a)
    common_cells <- intersect(cb2, ca2)
    if (length(common_cells) > 0) {
      # map back to original names
      idx_b <- match(common_cells, cb2); idx_a <- match(common_cells, ca2)
      M_before <- M_before[, idx_b, drop=FALSE]
      colnames(M_before) <- common_cells
      M_after  <- M_after[, idx_a, drop=FALSE]
      colnames(M_after) <- common_cells
    } else {
      message("No overlapping cells to compare for: ", bn)
      next
    }
  } else {
    M_before <- M_before[, common_cells, drop=FALSE]
    M_after  <- M_after[,  common_cells, drop=FALSE]
  }

  # summary numbers
  cells_before <- ncol(get_counts(obj_before))
  cells_after  <- ncol(get_counts(obj_after))
  genes_before <- nrow(get_counts(obj_before))
  genes_after  <- nrow(get_counts(obj_after))

  # mapping stats — infer from rownames (Ensembl → symbol)
  rn_b <- rownames(get_counts(obj_before))
  rn_a <- rownames(get_counts(obj_after))

  looks_ens <- grepl("^ENSG\\d+", rn_b)
  unmapped  <- sum(looks_ens & !(toupper(sub("\\.\\d+$","", rn_b)) %in% toupper(rn_a)))
  mapped    <- sum(looks_ens) - unmapped

  # collapse stats — how many symbols appear to be results of multiple Ensembl merged?
  # crude proxy: symbols that appeared >1 times in BEFORE when converted to symbol space
  # (we don’t recompute here; we can approximate by counting duplicated stripped Ensembl in BEFORE)
  # better: collect from your mapping log, but we estimate:
  ens_core <- toupper(sub("\\.\\d+$","", rn_b[looks_ens]))
  # symbols per core-id are unknown here; instead look at AFTER for potential "popular" symbols
  # estimate collapse size via duplicated symbols in AFTER is already resolved (unique); skip exact sizes.
  # But we can compute an upper bound by checking if genes_after << genes_before.
  collapsed_symbol_n <- max(0, genes_before - genes_after)

  # per-cell totals and detected genes
  csum_b <- Matrix::colSums(M_before)
  csum_a <- Matrix::colSums(M_after)
  gdet_b <- Matrix::colSums(M_before > 0)
  gdet_a <- Matrix::colSums(M_after  > 0)

  # correlations + diffs
  r2 <- function(x,y){
    if (length(x) < 2) return(NA_real_)
    stats::cor(x,y, method="pearson")^2
  }
  counts_r2 <- r2(csum_b, csum_a)
  counts_max_abs_diff <- max(abs(csum_b - csum_a))

  genes_r2 <- r2(gdet_b, gdet_a)
  genes_max_abs_diff <- max(abs(gdet_b - gdet_a))

  # marker sanity
  miss <- marker_symbols[!marker_symbols %in% rn_a]
  missing_markers <- if (length(miss)==0) "" else paste(miss, collapse=";")

  # write audit line
  line <- data.frame(
    sample_id = sid,
    cells_before = cells_before, cells_after = cells_after,
    genes_before = genes_before, genes_after = genes_after,
    mapped_n = mapped, unmapped_n = unmapped,
    collapsed_symbol_n = collapsed_symbol_n, max_collapse_size = NA_integer_,
    counts_r2 = counts_r2, counts_max_abs_diff = counts_max_abs_diff,
    genes_r2 = genes_r2, genes_max_abs_diff = genes_max_abs_diff,
    missing_markers = missing_markers,
    stringsAsFactors = FALSE
  )
  write.table(line, file=audit_tsv, sep="\t", quote=FALSE, row.names=FALSE, append=TRUE, col.names=!file.exists(audit_tsv))

  # quick plots (optional)
  p1 <- ggplot(data.frame(b=csum_b, a=csum_a), aes(b,a)) + geom_point(alpha=.3, size=.6) +
        ggtitle(paste0(sid, " counts/cell (before vs after); R2=", sprintf("%.4f", counts_r2))) +
        xlab("before total counts") + ylab("after total counts")
  ggsave(file.path(audit_dir, paste0(sid, "_counts_corr.png")), p1, width=5.5, height=4.5, dpi=150)

  p2 <- ggplot(data.frame(b=gdet_b, a=gdet_a), aes(b,a)) + geom_point(alpha=.3, size=.6) +
        ggtitle(paste0(sid, " genes detected/cell (before vs after); R2=", sprintf("%.4f", genes_r2))) +
        xlab("before genes>0") + ylab("after genes>0")
  ggsave(file.path(audit_dir, paste0(sid, "_genes_corr.png")), p2, width=5.5, height=4.5, dpi=150)

  # log a few examples of unmapped
  if (unmapped > 0) {
    ex <- head(rn_b[looks_ens & !(toupper(sub("\\.\\d+$","", rn_b)) %in% toupper(rn_a))], 20)
    writeLines(paste0("Unmapped examples (", sid, "): ", paste(ex, collapse=", ")),
               con = file.path(audit_dir, paste0(sid, "_unmapped_examples.txt")))
  }

  cat("Audited:", sid, "\n")
}
cat("Audit written to:", audit_tsv, "\n")


#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(dplyr); library(ggplot2); library(stringr)
})

# --- CONFIGURATION ---
before_dir <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds"        # original RDS
after_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds_hgnc_2"   # remapped RDS
audit_dir  <- file.path(dirname(after_dir), "hgnc_audit"); dir.create(audit_dir, showWarnings=FALSE, recursive=TRUE)
audit_tsv  <- file.path(audit_dir, "hgnc_audit_summary.tsv")

marker_symbols <- c(
  # mitochondria / ribo
  "MT-ND1","MT-CO1","RPLP0","RPS18",
  # immune / endothelium
  "CD3E","MS4A1","PECAM1","AIF1",
  # fibro/tumor-ish
  "COL1A1","FBLN1","MKI67"
)

# --- HELPER FUNCTIONS ---

# Lấy ma trận 'counts' từ RNA Assay
get_counts <- function(obj){
  # Hàm này chỉ dùng cú pháp slot cũ, tương thích với Seurat v4
  DefaultAssay(obj) <- "RNA"
  # Tự động lấy slot="counts". Seurat cũ không dùng layer, nên không cần LayerNames
  return(LayerData(obj, assay="RNA", layer="counts"))
}

# Hàm tính R-squared (Correlation squared)
r2 <- function(x,y){
  if (length(x) < 2) return(NA_real_)
  stats::cor(x,y, method="pearson")^2
}

# --- SECTION 1: DỮ LIỆU & CĂN CHỈNH MA TRẬN (ALIGNMENT) ---

align_and_prepare_matrices <- function(bf, af) {
  # 1. Tải dữ liệu
  obj_before <- readRDS(bf)
  obj_after  <- readRDS(af)
  
  M_before <- get_counts(obj_before)
  M_after  <- get_counts(obj_after)

  # 2. Đồng nhất Hạt nhân (Align Cells)
  cells_b <- colnames(M_before); cells_a <- colnames(M_after)
  common_cells <- intersect(cells_b, cells_a)
  
  if (length(common_cells) == 0) {
    # Thử loại bỏ hậu tố mẫu nếu không có tên trùng khớp trực tiếp
    strip_suf <- function(x) sub("-.*$", "", x)
    cb2 <- strip_suf(cells_b); ca2 <- strip_suf(cells_a)
    common_cells <- intersect(cb2, ca2)
    
    if (length(common_cells) > 0) {
      idx_b <- match(common_cells, cb2); idx_a <- match(common_cells, ca2)
      M_before <- M_before[, idx_b, drop=FALSE]; colnames(M_before) <- common_cells
      M_after  <- M_after[, idx_a, drop=FALSE]; colnames(M_after) <- common_cells
    } else {
      # Không có tế bào trùng khớp
      return(NULL)
    }
  } else {
    M_before <- M_before[, common_cells, drop=FALSE]
    M_after  <- M_after[,  common_cells, drop=FALSE]
  }

  # 3. Chuẩn bị metadata và ma trận đã căn chỉnh để chuyển sang hàm tính toán
  sid <- if (!is.null(obj_after$sample_id)) unique(as.character(obj_after$sample_id))[1] else tools::file_path_sans_ext(basename(af))
  
  return(list(
    sid = sid,
    M_before = M_before,
    M_after = M_after,
    cells_before_total = ncol(get_counts(obj_before)),
    genes_before_total = nrow(get_counts(obj_before)),
    genes_after_total = nrow(get_counts(obj_after))
  ))
}

# --- SECTION 2: TÍNH TOÁN VÀ GHI AUDIT (CALCULATE & WRITE) ---

calculate_and_write_audit <- function(prep_list, audit_tsv, audit_dir, marker_symbols) {
  if (is.null(prep_list)) return(NULL)

  # Trích xuất dữ liệu
  sid <- prep_list$sid
  M_before <- prep_list$M_before
  M_after <- prep_list$M_after
  cells_before <- prep_list$cells_before_total
  genes_before <- prep_list$genes_before_total
  genes_after <- prep_list$genes_after_total
  cells_after <- ncol(M_after) 

  # --- Thống kê Mapping Gen ---
  rn_b <- rownames(M_before); rn_a <- rownames(M_after)

  looks_ens <- grepl("^ENSG\\d+", rn_b)
  # Unmapped: Ensembl ID trong 'before' mà không tìm thấy symbol trong 'after'
  unmapped  <- sum(looks_ens & !(toupper(sub("\\.\\d+$","", rn_b)) %in% toupper(rn_a)))
  mapped    <- sum(looks_ens) - unmapped
  collapsed_symbol_n <- max(0, genes_before - genes_after) # ước tính

  # --- Tính toán Tương quan (Correlation) ---
  csum_b <- Matrix::colSums(M_before); csum_a <- Matrix::colSums(M_after)
  gdet_b <- Matrix::colSums(M_before > 0); gdet_a <- Matrix::colSums(M_after  > 0)
  
  counts_r2 <- r2(csum_b, csum_a)
  counts_max_abs_diff <- max(abs(csum_b - csum_a))
  genes_r2 <- r2(gdet_b, gdet_a)
  genes_max_abs_diff <- max(abs(gdet_b - gdet_a))

  # --- Kiểm tra Marker Genes ---
  miss <- marker_symbols[!marker_symbols %in% rn_a]
  missing_markers <- if (length(miss)==0) "" else paste(miss, collapse=";")

  # --- Ghi kết quả Audit (audit line) ---
  line <- data.frame(
    sample_id = sid,
    cells_before = cells_before, cells_after = cells_after,
    genes_before = genes_before, genes_after = genes_after,
    mapped_n = mapped, unmapped_n = unmapped,
    collapsed_symbol_n = collapsed_symbol_n, max_collapse_size = NA_integer_,
    counts_r2 = counts_r2, counts_max_abs_diff = counts_max_abs_diff,
    genes_r2 = genes_r2, genes_max_abs_diff = genes_max_abs_diff,
    missing_markers = missing_markers,
    stringsAsFactors = FALSE
  )
  # Ghi vào TSV (append=TRUE)
  write.table(line, file=audit_tsv, sep="\t", quote=FALSE, row.names=FALSE, append=TRUE, col.names=!file.exists(audit_tsv))

  # --- Vẽ biểu đồ và Ghi log ---
  # Plots
  p1 <- ggplot(data.frame(b=csum_b, a=csum_a), aes(b,a)) + geom_point(alpha=.3, size=.6) +
        ggtitle(paste0(sid, " counts/cell (before vs after); R2=", sprintf("%.4f", counts_r2))) +
        xlab("before total counts") + ylab("after total counts")
  ggsave(file.path(audit_dir, paste0(sid, "_counts_corr.png")), p1, width=5.5, height=4.5, dpi=150)

  p2 <- ggplot(data.frame(b=gdet_b, a=gdet_a), aes(b,a)) + geom_point(alpha=.3, size=.6) +
        ggtitle(paste0(sid, " genes detected/cell (before vs after); R2=", sprintf("%.4f", genes_r2))) +
        xlab("before genes>0") + ylab("after genes>0")
  ggsave(file.path(audit_dir, paste0(sid, "_genes_corr.png")), p2, width=5.5, height=4.5, dpi=150)

  # Log unmapped
  if (unmapped > 0) {
    ex <- head(rn_b[looks_ens & !(toupper(sub("\\.\\d+$","", rn_b)) %in% toupper(rn_a))], 20)
    writeLines(paste0("Unmapped examples (", sid, "): ", paste(ex, collapse=", ")),
               con = file.path(audit_dir, paste0(sid, "_unmapped_examples.txt")))
  }

  cat("Audited:", sid, "\n")
}

# --- MAIN EXECUTION ---

# 1. Khởi tạo header cho TSV nếu chưa tồn tại
if (!file.exists(audit_tsv)) {
  write.table(data.frame(
    sample_id=character(), cells_before=integer(), cells_after=integer(),
    genes_before=integer(), genes_after=integer(), mapped_n=integer(),
    unmapped_n=integer(), collapsed_symbol_n=integer(), max_collapse_size=integer(),
    counts_r2=double(), counts_max_abs_diff=double(), genes_r2=double(),
    genes_max_abs_diff=double(), missing_markers=character(),
    stringsAsFactors=FALSE
  ), file=audit_tsv, sep="\t", quote=FALSE, row.names=FALSE)
}

# 2. Vòng lặp chính điều phối
after_files <- list.files(after_dir, pattern="\\.rds$", full.names=TRUE)
for (af in after_files) {
  bn <- basename(af)
  bf <- file.path(before_dir, bn)
  
  if (!file.exists(bf)) {
    message("Skip (no BEFORE match): ", bn)
    next
  }
  
  # A. Chuẩn bị và căn chỉnh
  prep_data <- align_and_prepare_matrices(bf, af)
  
  if (is.null(prep_data)) {
    message("Skipping sample due to no common cells: ", bn)
    next
  }
  
  # B. Tính toán và ghi audit
  calculate_and_write_audit(prep_data, audit_tsv, audit_dir, marker_symbols)
}

cat("Audit written to:", audit_tsv, "\n")