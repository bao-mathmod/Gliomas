#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(FNN)       # kNN for metrics
  library(igraph)    # connected components
  library(mclust)    # adjustedRandIndex (no dist())
  library(uwot)      # Required for uwot/UMAP operations
})

## ===================== USER SETTINGS =====================
# NOTE: Replace 'scvi_obj' with the name of your actual Seurat object
# that contains the reduction specified by 'reduc'.
obj           <- scvi_obj
reduc         <- "integrated.scvi"      # name of integrated reduction (e.g., 'pca', 'harmony', 'integrated.scvi')
subset_n      <- 40000                  # 20k–50k is fine for tuning, balances speed/accuracy

ncol(Embeddings(obj, reduc))
# If you prefer to explicitly set dims, set dims_use_list <- list(1:30, 1:50, 1:70, 1:100)
dims_use_list <- NULL                     # auto-generate from available PCs if NULL
k_grid        <- c(50, 100, 150, 200)   # UMAP n.neighbors (more robust at atlas scale)
min_dist_grid <- c(0.15, 0.30, 0.50)    # UMAP min.dist
metric_grid   <- c("euclidean", "cosine") # try both, pick best
n_epochs      <- 500                    # more epochs helps large graphs/better convergence
seeds         <- c(1234, 4321)          # for stability vs seed
k_score       <- 15                     # k used in scoring metrics (e.g., 15-30)
n_eval_pairs  <- 5000                   # pairs for Jaccard metrics (subset for speed)
annoy_trees   <- 75                     # Annoy index trees (higher = better quality, slower index build)
annoy_search_k<- 6000                   # Annoy search_k (higher = better recall/accuracy, slower query)
force_threads <- NULL                    # e.g., 20; set NULL to let uwot autodetect
centers_guess <- NULL                    # kmeans centers; NULL = auto heuristic
set.seed(1234)
## =========================================================

## ---------- Sanity checks & auto dims ----------
if (!exists(deparse(substitute(obj)))) {
  stop("Object 'obj' (currently set to scvi_obj) is not defined. Please load your Seurat object.")
}

if (is.null(reduc) || !(reduc %in% Reductions(obj))) {
  stop(sprintf("Reduction '%s' not found in object. Available: %s",
              reduc, paste0(names(obj@reductions), collapse = ", ")))
}
Xr <- Embeddings(obj, reduc)
n_pcs_available <- ncol(Xr)

# Build dims_use_list automatically if not provided
if (is.null(dims_use_list)) {
  candidates <- c(30, 50, 70, 100)
  candidates <- candidates[candidates <= n_pcs_available]
  if (length(candidates) == 0) candidates <- min(30, n_pcs_available)
  
  # Ensure candidates are unique and sorted
  candidates <- unique(sort(candidates))

  dims_use_list <- lapply(candidates, function(p) 1:p)
  message("Auto dims_use_list: ", paste(sprintf("1:%d", sapply(dims_use_list, length)), collapse=", "))
}

## ---------- Helper: stratified sample over batches if present ----------
stratified_sample <- function(obj, n_total) {
  metas <- obj@meta.data
  by <- NULL
  if ("orig.ident" %in% colnames(metas)) by <- "orig.ident"
  else if ("batch" %in% colnames(metas)) by <- "batch"
  if (is.null(by)) {
    return(sample(colnames(obj), min(n_total, ncol(obj))))
  } else {
    tbl <- table(metas[[by]])
    frac <- n_total / sum(tbl)
    take <- pmax(1, round(as.numeric(tbl) * frac))
    cells <- unlist(mapply(function(level, m) {
      pool <- rownames(metas)[metas[[by]] == level]
      sample(pool, min(m, length(pool)))
    }, names(tbl), take, SIMPLIFY = FALSE), use.names = FALSE)
    return(unique(cells)[seq_len(min(length(unique(cells)), n_total))])
  }
}

## ---------- Metric helpers (all scalable) ----------

# Counts the number of connected components in the kNN graph (k=k_score)
comp_count <- function(emb2d, k = 15) {
  # Use exact kNN search for scoring metrics
  nn <- FNN::get.knn(emb2d, k = k)$nn.index
  
  # Create edges: for each cell (row i), link it to its k neighbors
  edges <- do.call(rbind, lapply(seq_len(nrow(emb2d)), function(i)
    cbind(i, nn[i, ])
  ))
  
  # Build undirected graph and count components
  g <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)
  length(igraph::components(g)$csize)
}

# Jaccard Overlap (Proxy for Trustworthiness/Manifold Preservation)
jaccard_overlap_highD_lowD <- function(highD, lowD, k = 15, n_eval = 5000) {
  n <- nrow(highD)
  idx <- sample.int(n, min(n_eval, n))
  
  # Find kNN in HighD space (e.g., 30-100 PCs)
  kn_hi <- FNN::get.knn(highD[idx, , drop=FALSE], k = k)$nn.index
  
  # Find kNN in LowD space (UMAP coordinates)
  kn_lo <- FNN::get.knn(lowD[idx, , drop=FALSE], k = k)$nn.index
  
  # Calculate mean Jaccard Index across sampled cells
  mean(sapply(seq_along(idx), function(i) {
    length(intersect(kn_hi[i, ], kn_lo[i, ])) / length(union(kn_hi[i, ], kn_lo[i, ]))
  }))
}

# Seed Stability (Jaccard)
seed_nn_jaccard <- function(embA, embB, k = 15, n_eval = 5000) {
  n <- nrow(embA)
  idx <- sample.int(n, min(n_eval, n))
  
  # Find kNN for UMAP run A
  knA <- FNN::get.knn(embA[idx, , drop=FALSE], k = k)$nn.index
  
  # Find kNN for UMAP run B
  knB <- FNN::get.knn(embB[idx, , drop=FALSE], k = k)$nn.index
  
  # Calculate mean Jaccard Index across sampled cells for 2D neighborhood overlap
  mean(sapply(seq_along(idx), function(i) {
    length(intersect(knA[i, ], knB[i, ])) / length(union(knA[i, ], knB[i, ]))
  }))
}

## ---------- Build subset ----------
sub_cells <- stratified_sample(obj, subset_n)
obj_sub   <- subset(obj, cells = sub_cells)
message(sprintf("Tuning on subset: %d cells; reduction: '%s' with %d PCs available",
                ncol(obj_sub), reduc, n_pcs_available))

## ---------- Heuristic for kmeans centers (if not provided) ----------
if (is.null(centers_guess)) {
  # Scales mildly with log(size); bounded for stability
  centers_guess <- max(25, min(60, round(log10(ncol(obj_sub)) * 12)))
}
message(sprintf("kmeans centers heuristic: %d", centers_guess))

## ---------- Grid search ----------
results <- list()
ctr <- 1

for (dims in dims_use_list) {
  # clip dims defensively
  dims <- dims[dims <= n_pcs_available]
  if (length(dims) < 5) {
    message(sprintf("Skipping dims=1:%d, fewer than 5 dimensions.", length(dims)))
    next
  }
  X_high_full <- Embeddings(obj_sub, reduc)[, dims, drop = FALSE]

  for (metric_umap in metric_grid) {
    for (k in k_grid) {
      for (md in min_dist_grid) {

        # Run two seeds; robust to uwot/annoy hiccups via tryCatch
        emb_list <- list()
        ok <- TRUE
        for (s in seeds) {
          set.seed(s)
          
          # Construct a unique reduction name for this parameter set and seed
          reduction_name_s <- paste0("umap.tune.k", k, ".md", md,
                                     ".d", length(dims), ".", metric_umap, ".s", s)
          
          emb_mat <- tryCatch({
            obj_tmp <- RunUMAP(
              obj_sub,
              reduction      = reduc,
              dims           = dims,
              metric         = metric_umap,
              n.neighbors    = k,
              min.dist       = md,
              init           = "pca",          # avoid spectral fallback for large data
              n.epochs       = n_epochs,
              n.trees        = annoy_trees,    # Annoy params
              search_k       = annoy_search_k,
              umap.method    = "uwot",
              n_threads      = force_threads,  # user-defined thread count
              reduction.name = reduction_name_s,
              verbose        = FALSE
            )
            Embeddings(obj_tmp, reduction_name_s)
          }, error = function(e) {
            message(sprintf("UMAP failed for dims=1:%d k=%d md=%.2f metric=%s seed=%d; skipping. Reason: %s",
                            length(dims), k, md, metric_umap, s, e$message))
            return(NULL)
          })
          if (is.null(emb_mat)) { ok <- FALSE; break }
          emb_list[[as.character(s)]] <- emb_mat
        }
        if (!ok) next

        # Use seed1 as reference
        emb1 <- emb_list[[as.character(seeds[1])]]
        
        # Guard against NA/Inf rows (should be rare with uwot defaults)
        keep1 <- which(rowSums(!is.finite(emb1)) == 0)
        if (length(keep1) < nrow(emb1)) {
          emb1 <- emb1[keep1, , drop = FALSE]
          X_high <- X_high_full[keep1, , drop = FALSE]
          # Ensure both seeds are subsetted equally
          emb2 <- emb_list[[as.character(seeds[2])]][keep1, , drop = FALSE]
          n_eval_pairs_adj <- min(n_eval_pairs, length(keep1))
        } else {
          X_high <- X_high_full
          emb2 <- emb_list[[as.character(seeds[2])]]
          n_eval_pairs_adj <- n_eval_pairs
        }

        # Metrics: components, HighD↔2D Jaccard, trust proxy
        components  <- comp_count(emb1, k = k_score)
        jaccard15   <- jaccard_overlap_highD_lowD(X_high, emb1, k = k_score, n_eval = n_eval_pairs_adj)
        trust15     <- jaccard15  # same proxy; keep a separate column for readability

        # Stability vs second seed
        seedJ <- seed_nn_jaccard(emb1, emb2, k = k_score, n_eval = n_eval_pairs_adj)

        # Clustering-based ARI (no dist())
        set.seed(1)
        # Use adjusted n-start for better stability
        nstart_km <- min(25, max(3, floor(sqrt(centers_guess))))
        km1 <- kmeans(emb1, centers = centers_guess, nstart = nstart_km, iter.max = 50)$cluster
        km2 <- kmeans(emb2, centers = centers_guess, nstart = nstart_km, iter.max = 50)$cluster
        ARI <- mclust::adjustedRandIndex(km1, km2)

        results[[ctr]] <- data.frame(
          dims           = length(dims),
          metric         = metric_umap,
          n_neighbors    = k,
          min_dist       = md,
          components     = as.integer(components),
          jaccard15      = jaccard15,
          trust15        = trust15,
          seed2D_jacc    = seedJ,
          ARI            = ARI,
          stringsAsFactors = FALSE
        )
        ctr <- ctr + 1

        message(sprintf("Tuned d=%-3d metric=%-9s k=%-3d md=%.2f -> comps=%-3d J15=%.3f ARI=%.3f J2D=%.3f",
                        length(dims), metric_umap, k, md, components, jaccard15, ARI, seedJ))
      }
    }
  }
}

if (length(results) == 0) stop("No successful UMAP runs in the grid. Check logs/params.")
res <- do.call(rbind, results)

## ---------- Composite score & ranking ----------
# The scoring prioritizes preserving the high-D manifold (jaccard/trust)
# and visualization stability (ARI/seedJ), while penalizing fragmentation.
res$score <- (0.35 * res$jaccard15) +
             (0.35 * res$trust15)  +
             (0.15 * pmin(1, res$ARI)) + # Cap ARI at 1
             (0.15 * res$seed2D_jacc) -
             0.05 * pmax(0, res$components - 5) # Penalize more than 5 disconnected parts

res <- res[order(-res$score), ]
rownames(res) <- NULL

## ---------- Save & show ----------
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- "/mnt/18T/chibao/gliomas/data/upstream/sc_sn/official/integrated_v5_optimized"
out_filename <- paste0("umap_grid_results_", stamp, ".csv")
out_csv <- file.path(output_dir, out_filename)
write.csv(res, out_csv, row.names = FALSE)
message(sprintf("Saved results: %s", out_csv))
print(head(res, 15))

## ---------- Apply best config to FULL OBJECT (optional) ----------
# To run the best configuration on your full Seurat object (not just the subset),
# uncomment the following block and ensure the 'best_obj' variable holds your
# full Seurat object (e.g., replace 'harmony_obj' with 'obj').

# best <- res[1, ]
# set.seed(1234)
# best_dims <- seq_len(best$dims)
# 
# # Replace 'harmony_obj' with the name of your full Seurat object (e.g., 'obj')
# best_obj <- obj 
# 
# best_obj <- RunUMAP(
#   best_obj,
#   reduction      = reduc,
#   dims           = best_dims,
#   metric         = best$metric,
#   n.neighbors    = best$n_neighbors,
#   min.dist       = best$min_dist,
#   init           = "pca",
#   n.epochs       = n_epochs,
#   n.trees        = annoy_trees,
#   search_k       = annoy_search_k,
#   umap.method    = "uwot",
#   reduction.name = sprintf("umap.optimized.k%d.p%d.md%.2f.%s",
#                            best$n_neighbors, best$dims, best$min_dist, best$metric),
#   verbose        = TRUE
# )
# message(sprintf("Final UMAP on full object complete: dims=1:%d, k=%d, min.dist=%.2f, metric=%s",
#                 best$dims, best$n_neighbors, best$min_dist, best$metric))

# You would then save the updated object:
# saveRDS(best_obj, "pbmcsca_with_optimized_umap.rds")
