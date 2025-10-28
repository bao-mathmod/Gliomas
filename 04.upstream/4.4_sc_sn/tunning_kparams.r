# #!/usr/bin/env Rscript
# suppressPackageStartupMessages({
#   library(Seurat)
#   library(FNN)       # kNN for metrics
#   library(igraph)    # connected components
#   library(mclust)    # adjustedRandIndex (no dist())
# })

# ## ===================== USER SETTINGS =====================
# obj             <- scvi_obj          # <-- your Seurat object
# reduc           <- "integrated.scvi"            # name of integrated reduction
# subset_n        <- 40000                # 20k–50k is fine
# # If you prefer to explicitly set dims, set dims_use_list <- list(1:30, 1:50, 1:70, 1:100)
# dims_use_list   <- NULL                 # auto-generate from available PCs if NULL
# k_grid          <- c(50, 100, 150, 200) # UMAP n.neighbors (more robust at atlas scale)
# min_dist_grid   <- c(0.15, 0.30, 0.50)  # UMAP min.dist
# metric_grid     <- c("euclidean", "cosine")  # try both, pick best
# n_epochs        <- 500                  # more epochs helps large graphs
# seeds           <- c(1234, 4321)       # for stability vs seed
# k_score         <- 15                   # k used in scoring metrics
# n_eval_pairs    <- 5000                 # pairs for Jaccard metrics
# annoy_trees     <- 75                   # Annoy index trees
# annoy_search_k  <- 6000                 # Annoy search_k (recall)
# force_threads   <- NULL                 # e.g., 20; set NULL to let uwot autodetect
# centers_guess   <- NULL                 # kmeans centers; NULL = auto heuristic
# set.seed(1234)
# ## =========================================================

# ## ---------- Sanity checks & auto dims ----------
# if (is.null(reduc) || !(reduc %in% Reductions(obj))) {
#   stop(sprintf("Reduction '%s' not found in object. Available: %s",
#                reduc, paste0(names(obj@reductions), collapse = ", ")))
# }
# Xr <- Embeddings(obj, reduc)
# n_pcs_available <- ncol(Xr)

# # Build dims_use_list automatically if not provided
# if (is.null(dims_use_list)) {
#   candidates <- c(30, 50, 70, 100)
#   candidates <- candidates[candidates <= n_pcs_available]
#   if (length(candidates) == 0) candidates <- min(30, n_pcs_available)
#   dims_use_list <- lapply(candidates, function(p) 1:p)
#   message("Auto dims_use_list: ", paste(sprintf("1:%d", sapply(dims_use_list, length)), collapse=", "))
# }

# ## ---------- Helper: stratified sample over batches if present ----------
# stratified_sample <- function(obj, n_total) {
#   metas <- obj@meta.data
#   by <- NULL
#   if ("orig.ident" %in% colnames(metas)) by <- "orig.ident"
#   else if ("batch" %in% colnames(metas)) by <- "batch"
#   if (is.null(by)) {
#     return(sample(colnames(obj), min(n_total, ncol(obj))))
#   } else {
#     tbl <- table(metas[[by]])
#     frac <- n_total / sum(tbl)
#     take <- pmax(1, round(as.numeric(tbl) * frac))
#     cells <- unlist(mapply(function(level, m) {
#       pool <- rownames(metas)[metas[[by]] == level]
#       sample(pool, min(m, length(pool)))
#     }, names(tbl), take, SIMPLIFY = FALSE), use.names = FALSE)
#     unique(cells)[seq_len(min(length(unique(cells)), n_total))]
#   }
# }

# ## ---------- Metric helpers (all scalable) ----------
# comp_count <- function(emb2d, k = 15) {
#   nn <- FNN::get.knn(emb2d, k = k)
#   edges <- do.call(rbind, lapply(seq_len(nrow(emb2d)), function(i)
#     cbind(i, nn$nn.index[i, ])
#   ))
#   g <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)
#   length(igraph::components(g)$csize)
# }

# jaccard_overlap_highD_lowD <- function(highD, lowD, k = 15, n_eval = 5000) {
#   n <- nrow(highD)
#   idx <- sample.int(n, min(n_eval, n))
#   kn_hi <- FNN::get.knn(highD[idx, , drop=FALSE], k = k)$nn.index
#   kn_lo <- FNN::get.knn(lowD[idx,  , drop=FALSE], k = k)$nn.index
#   mean(sapply(seq_along(idx), function(i) {
#     length(intersect(kn_hi[i, ], kn_lo[i, ])) / length(union(kn_hi[i, ], kn_lo[i, ]))
#   }))
# }

# seed_nn_jaccard <- function(embA, embB, k = 15, n_eval = 5000) {
#   n <- nrow(embA)
#   idx <- sample.int(n, min(n_eval, n))
#   knA <- FNN::get.knn(embA[idx, , drop=FALSE], k = k)$nn.index
#   knB <- FNN::get.knn(embB[idx, , drop=FALSE], k = k)$nn.index
#   mean(sapply(seq_along(idx), function(i) {
#     length(intersect(knA[i, ], knB[i, ])) / length(union(knA[i, ], knB[i, ]))
#   }))
# }

# ## ---------- Build subset ----------
# sub_cells <- stratified_sample(obj, subset_n)
# obj_sub   <- subset(obj, cells = sub_cells)
# message(sprintf("Tuning on subset: %d cells; reduction: '%s' with %d PCs available",
#                 ncol(obj_sub), reduc, n_pcs_available))

# ## ---------- Heuristic for kmeans centers (if not provided) ----------
# if (is.null(centers_guess)) {
#   # Scales mildly with size; bounded for stability
#   centers_guess <- max(25, min(60, round(log10(ncol(obj_sub)) * 12)))
# }
# message(sprintf("kmeans centers heuristic: %d", centers_guess))

# ## ---------- Grid search ----------
# results <- list()
# ctr <- 1

# for (dims in dims_use_list) {
#   # clip dims defensively in case of future changes
#   dims <- dims[dims <= n_pcs_available]
#   if (length(dims) < 5) next
#   X_high_full <- Embeddings(obj_sub, reduc)[, dims, drop = FALSE]

#   for (metric_umap in metric_grid) {
#     for (k in k_grid) {
#       for (md in min_dist_grid) {

#         # Run two seeds; robust to uwot/annoy hiccups via tryCatch
#         emb_list <- list()
#         ok <- TRUE
#         for (s in seeds) {
#           set.seed(s)
#           emb_mat <- tryCatch({
#             obj_tmp <- RunUMAP(
#               obj_sub,
#               reduction      = reduc,
#               dims           = dims,
#               metric         = metric_umap,
#               n.neighbors    = k,
#               min.dist       = md,
#               init           = "pca",         # avoid spectral fallback
#               n.epochs       = n_epochs,
#               n.trees        = annoy_trees,   # Annoy params
#               search_k       = annoy_search_k,
#               umap.method    = "uwot",
#               # pass threads if you want to force:
#               # n_threads    = force_threads,
#               reduction.name = paste0("umap.tune.k", k, ".md", md,
#                                       ".d", length(dims), ".", metric_umap, ".s", s),
#               verbose        = FALSE
#             )
#             Embeddings(obj_tmp, paste0("umap.tune.k", k, ".md", md,
#                                        ".d", length(dims), ".", metric_umap, ".s", s))
#           }, error = function(e) {
#             message(sprintf("UMAP failed for dims=1:%d k=%d md=%.2f metric=%s seed=%d; skipping. Reason: %s",
#                             length(dims), k, md, metric_umap, s, e$message))
#             return(NULL)
#           })
#           if (is.null(emb_mat)) { ok <- FALSE; break }
#           emb_list[[as.character(s)]] <- emb_mat
#         }
#         if (!ok) next

#         # Use seed1 as reference
#         emb1 <- emb_list[[as.character(seeds[1])]]
#         # Guard against NA/Inf rows
#         keep1 <- which(rowSums(!is.finite(emb1)) == 0)
#         if (length(keep1) < nrow(emb1)) emb1 <- emb1[keep1, , drop = FALSE]
#         X_high <- X_high_full[keep1, , drop = FALSE]

#         # Metrics: components, HighD↔2D Jaccard, trust proxy
#         components  <- comp_count(emb1, k = k_score)
#         jaccard15   <- jaccard_overlap_highD_lowD(X_high, emb1, k = k_score, n_eval = n_eval_pairs)
#         trust15     <- jaccard15  # same proxy; keep a separate column for readability

#         # Stability vs second seed
#         emb2 <- emb_list[[as.character(seeds[2])]][keep1, , drop = FALSE]
#         seedJ <- seed_nn_jaccard(emb1, emb2, k = k_score, n_eval = n_eval_pairs)

#         # Clustering-based ARI (no dist())
#         set.seed(1)
#         km1 <- kmeans(emb1, centers = centers_guess, nstart = 3)$cluster
#         km2 <- kmeans(emb2, centers = centers_guess, nstart = 3)$cluster
#         ARI <- mclust::adjustedRandIndex(km1, km2)

#         results[[ctr]] <- data.frame(
#           dims         = length(dims),
#           metric       = metric_umap,
#           n_neighbors  = k,
#           min_dist     = md,
#           components   = as.integer(components),
#           jaccard15    = jaccard15,
#           trust15      = trust15,
#           seed2D_jacc  = seedJ,
#           ARI          = ARI,
#           stringsAsFactors = FALSE
#         )
#         ctr <- ctr + 1

#         message(sprintf("Tuned d=%-3d metric=%-9s k=%-3d md=%.2f -> comps=%-3d J15=%.3f ARI=%.3f J2D=%.3f",
#                         length(dims), metric_umap, k, md, components, jaccard15, ARI, seedJ))
#       }
#     }
#   }
# }

# if (length(results) == 0) stop("No successful UMAP runs in the grid. Check logs/params.")
# res <- do.call(rbind, results)

# ## ---------- Composite score & ranking ----------
# # Weights: HighD↔2D Jaccard & trust (most important), then ARI & seed2D Jaccard; penalize fragmentation
# res$score <- (0.35 * res$jaccard15) +
#              (0.35 * res$trust15)   +
#              (0.15 * pmin(1, res$ARI)) +
#              (0.15 * res$seed2D_jacc) -
#              0.05 * pmax(0, res$components - 5)

# res <- res[order(-res$score), ]
# rownames(res) <- NULL

# ## ---------- Save & show ----------
# stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
# out_csv <- paste0("umap_grid_results_", stamp, ".csv")
# write.csv(res, out_csv, row.names = FALSE)
# message(sprintf("Saved results: %s", out_csv))
# print(head(res, 15))

# ## ---------- Apply best config to FULL OBJECT (optional) ----------
# # best <- res[1, ]
# # set.seed(1234)
# # best_dims <- seq_len(best$dims)
# # harmony_obj <- RunUMAP(
# #   harmony_obj,
# #   reduction      = reduc,
# #   dims           = best_dims,
# #   metric         = best$metric,
# #   n.neighbors    = best$n_neighbors,
# #   min.dist       = best$min_dist,
# #   init           = "pca",
# #   n.epochs       = n_epochs,
# #   n.trees        = annoy_trees,
# #   search_k       = annoy_search_k,
# #   umap.method    = "uwot",
# #   reduction.name = sprintf("umap.harmony.k%d.p%d.md%.2f.%s",
# #                            best$n_neighbors, best$dims, best$min_dist, best$metric),
# #   verbose        = TRUE
# # )
# # message(sprintf("Final UMAP on full object with dims=1:%d, k=%d, min.dist=%.2f, metric=%s",
# #                 best$dims, best$n_neighbors, best$min_dist, best$metric))

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(FNN)
  library(igraph)
})

# ----- Inputs -----
obj         <- scvi_obj                 # your Seurat with scVI reduction
reduc       <- "integrated.scvi"        # scVI latent space
L           <- ncol(Embeddings(obj, reduc))  # latent_dim, often 30
dims_use    <- 1:L
subset_n    <- 40000                    # tune on 30–50k cells
k_grid      <- c(60, 80, 100, 120, 150)
k_eval      <- 15                       # for metrics
n_eval_pairs<- 5000
set.seed(1234)

# ----- Stratified subset by batch/project -----
stratified_sample <- function(obj, n_total){
  md <- obj@meta.data
  by <- if ("orig.ident" %in% names(md)) "orig.ident" else if ("batch" %in% names(md)) "batch" else NULL
  if (is.null(by)) return(sample(colnames(obj), min(n_total, ncol(obj))))
  tab <- table(md[[by]]); frac <- n_total/sum(tab)
  take <- pmax(1, round(as.numeric(tab)*frac))
  cells <- unlist(mapply(function(level,m){
    pool <- rownames(md)[md[[by]]==level]; sample(pool, min(m, length(pool)))
  }, names(tab), take, SIMPLIFY = FALSE), use.names = FALSE)
  unique(cells)[seq_len(min(length(unique(cells)), n_total))]
}

sub_cells <- stratified_sample(obj, subset_n)
obj_sub   <- subset(obj, cells = sub_cells)
Z         <- Embeddings(obj_sub, reduc)[, dims_use, drop=FALSE]
stopifnot(all(is.finite(Z)))
message(sprintf("Tuning k.param on %d cells in %s (dims 1:%d)", ncol(obj_sub), reduc, length(dims_use)))

# ----- Metrics -----
comp_count <- function(emb2d, k = 15){  # connected components in 2D-like space (we'll use Z here)
  nn <- FNN::get.knn(emb2d, k = k)
  edges <- do.call(rbind, lapply(seq_len(nrow(emb2d)), function(i) cbind(i, nn$nn.index[i,])))
  g <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)
  length(igraph::components(g)$csize)
}

# HighD kNN overlap between a candidate k and a strong baseline (k_ref = max(k_grid))
jaccard_vs_ref <- function(Z, k, k_ref, n_eval = 5000){
  n <- nrow(Z); idx <- sample.int(n, min(n_eval, n))
  kn_k   <- FNN::get.knn(Z[idx,,drop=FALSE], k = k)$nn.index
  kn_ref <- FNN::get.knn(Z[idx,,drop=FALSE], k = k_ref)$nn.index
  mean(sapply(seq_along(idx), function(i){
    length(intersect(kn_k[i,], kn_ref[i,])) / length(union(kn_k[i,], kn_ref[i,]))
  }))
}

# ----- Evaluate grid -----
# k_ref <- max(k_grid)
# results <- lapply(k_grid, function(k){
#   # components in a derived 2D graph built from Z (just for connectivity diagnosis)
#   comps <- comp_count(Z, k = k_eval)      # stable small-k snapshot of fragmentation
#   # neighborhood stability vs a strong baseline k_ref
#   jac   <- jaccard_vs_ref(Z, k = k, k_ref = k_ref, n_eval = n_eval_pairs)
#   data.frame(k.param = k, components@k15 = comps, jaccard_vs_kref = jac)
# })
# res <- do.call(rbind, results)
k_ref <- max(k_grid)
results <- lapply(k_grid, function(k){
  # components in a derived 2D graph built from Z (just for connectivity diagnosis)
  comps <- comp_count(Z, k = k_eval)      # stable small-k snapshot of fragmentation
  # neighborhood stability vs a strong baseline k_ref
  jac   <- jaccard_vs_ref(Z, k = k, k_ref = k_ref, n_eval = n_eval_pairs)
  data.frame(k.param = k, components.k15 = comps, jaccard_vs_kref = jac) # Changed '@' to '.' for safety
})
res <- do.call(rbind, results)
# ----- Choose k: first k with (components <= 5) and decent overlap -----
res$ok_conn <- res$components.k15 <= 5
thr_jac <- quantile(res$jaccard_vs_kref, 0.5, na.rm = TRUE)  # median as a simple threshold
res$ok_jac <- res$jaccard_vs_kref >= thr_jac
res$score  <- (0.6*scale(res$jaccard_vs_kref)) - (0.4*scale(res$components.k15))
res <- res[order(-res$ok_conn, -res$ok_jac, -res$score), ]

print(res)

best_k <- res$k.param[1]
message(sprintf("Chosen k.param = %d (components@15=%d, Jaccard_vs_%d=%.3f)",
                best_k, res$components.k15[1], k_ref, res$jaccard_vs_kref[1]))

# # ----- Apply on FULL object, build SNN and run UMAP from the graph -----
# scvi_obj <- FindNeighbors(
#   scvi_obj,
#   reduction  = reduc,
#   dims       = dims_use,
#   k.param    = best_k,
#   graph.name = "scvi_snn"
# )

# scvi_obj <- RunUMAP(
#   scvi_obj,
#   graph          = "scvi_snn",
#   init           = "pca",     # avoids spectral init
#   min.dist       = 0.30,
#   n.epochs       = 500,
#   umap.method    = "uwot",
#   reduction.name = sprintf("umap.scvi.k%d", best_k),
#   verbose        = TRUE
# )
