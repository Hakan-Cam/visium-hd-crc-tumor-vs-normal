#!/usr/bin/env Rscript
# validate_pipeline_logic.R
#
# NOT part of the deliverable pipeline. This is a throwaway smoke test,
# run once in the sandbox that built this repo, against small synthetic
# data shaped like the real thing -- to catch real Seurat/dplyr API bugs
# in 02_visium_hd_analysis.Rmd's logic *before* handing the notebook off,
# since the real GEO data and a few packages (spacexr, fgsea, msigdbr,
# arrow, qs2, scatterpie) aren't reachable from this environment.
#
# It exercises: SCTransform -> SketchData -> cluster -> ProjectData,
# the two-sample merge + joint sketch-clustering extension, FindAllMarkers,
# FindMarkers, and the human-gene regex filters -- i.e. everything in the
# notebook that doesn't require spacexr/fgsea/msigdbr/arrow or a real
# spatial image. It does NOT validate RCTD, GSEA, or scatterpie code,
# which needs manual review instead (done separately).

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})
set.seed(123)

make_fake_sample <- function(n_cells = 3000, n_genes = 800, condition, seed) {
  set.seed(seed)
  gene_names <- c(
    paste0("GENE", seq_len(n_genes - 15)),
    paste0("MT-", c("ND1","ND2","CO1","CO2","ATP6","CO3","ND3","ND4","ND5","ND6","CYB")),  # 11
    paste0("AC", sample(100000:999999, 4), ".", sample(1:3, 4, replace = TRUE))              # 4 "clone-named" genes
  )
  stopifnot(length(gene_names) == n_genes)

  m <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.1,
                              rand.x = function(n) rpois(n, lambda = 3) + 1)
  rownames(m) <- gene_names
  colnames(m) <- paste0(condition, "_bin", seq_len(n_cells))

  obj <- CreateSeuratObject(counts = m, assay = "Spatial8um")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$condition <- condition
  obj$patient <- "TEST"
  obj
}

cat("== Building two fake samples ==\n")
tumor  <- make_fake_sample(condition = "tumor", seed = 1)
normal <- make_fake_sample(condition = "normal_adjacent", seed = 2)
cat("OK:", ncol(tumor), "+", ncol(normal), "fake bins\n\n")

process_fake <- function(obj, resolution = 0.5, n_sketch_max = 1500) {
  assay_name <- "Spatial8um"
  DefaultAssay(obj) <- assay_name
  obj <- subset(obj, subset = nCount_Spatial8um >= 1 & percent.mt < 30)

  # SketchData needs to run on the RAW-counts assay, named explicitly --
  # leaving `assay` at its default picks up whatever DefaultAssay()
  # currently is, and running it against an SCT assay (i.e. sketching
  # *after* SCTransform) breaks: LeverageScore's internal Cells() lookup
  # on an SCTAssay expects a "counts" model and gets "data" instead. So:
  # sketch first, from the plain counts assay, then normalize only the
  # ~50K-cell sketch (cheap) rather than the full object (expensive) --
  # which is also just a more efficient order of operations.
  #
  # LeverageScore also needs normalized data + a variable-feature set to
  # actually score cells -- caught on real data, not by this test
  # originally (see the "sketch actually shrank" check right below,
  # added after that): skipping this doesn't error, it just leaves
  # SketchData unable to subsample at all, silently keeping every cell.
  obj <- NormalizeData(obj, assay = assay_name, verbose = FALSE) %>%
    FindVariableFeatures(assay = assay_name, verbose = FALSE)

  n_before_sketch <- ncol(obj)
  obj <- SketchData(obj, assay = assay_name, ncells = min(n_sketch_max, ncol(obj)),
                     method = "LeverageScore", sketched.assay = "sketch")
  stopifnot(
    "SketchData did not actually subsample -- sketch has the same cell count as the source" =
      ncol(obj[["sketch"]]) < n_before_sketch
  )

  DefaultAssay(obj) <- "sketch"
  obj <- SCTransform(obj, assay = "sketch", new.assay.name = "sketch", vars.to.regress = "percent.mt",
                      return.only.var.genes = FALSE, verbose = FALSE) %>%
    RunPCA(npcs = 10, verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, k.param = 30, verbose = FALSE) %>%
    FindClusters(resolution = resolution, verbose = FALSE) %>%
    RunUMAP(dims = 1:10, return.model = TRUE, verbose = FALSE)

  # normalization.method = "SCT" here means ProjectData normalizes the
  # FULL (un-sketched) assay itself as part of projecting -- we don't
  # need to (and shouldn't) pre-run SCTransform on the full object.
  obj <- ProjectData(
    object = obj, assay = assay_name, sketched.assay = "sketch",
    full.reduction = "pca.full", sketched.reduction = "pca", umap.model = "umap",
    dims = 1:10, normalization.method = "SCT",
    refdata = list(seurat_clusters = "seurat_clusters")
  )
  DefaultAssay(obj) <- assay_name
  obj <- NormalizeData(obj, assay = assay_name, verbose = FALSE) %>%
    ScaleData(assay = assay_name, verbose = FALSE)
  obj
}

cat("== Per-sample: SCTransform -> SketchData -> cluster -> ProjectData ==\n")
tumor  <- process_fake(tumor)
normal <- process_fake(normal)
cat("OK. tumor clusters:", nlevels(tumor$seurat_clusters),
    " normal clusters:", nlevels(normal$seurat_clusters), "\n")
stopifnot("seurat_clusters" %in% colnames(tumor@meta.data))
stopifnot(!is.null(tumor@reductions$full.umap))
cat("full.umap reduction present: TRUE\n\n")

cat("== Human gene regex filters ==\n")
bad <- rownames(tumor)[str_detect(rownames(tumor), "^A[CL][0-9]{5,6}\\.[0-9]+$")]
cat("Clone-named genes matched:", length(bad), "(expect 4)\n")
stopifnot(length(bad) == 4)
mt_pct_range <- range(tumor$percent.mt)
cat("percent.mt range:", paste(round(mt_pct_range, 2), collapse = " - "), "(expect > 0)\n")
stopifnot(mt_pct_range[2] > 0)
cat("OK\n\n")

cat("== FindAllMarkers on sketch-projected clusters ==\n")
Idents(tumor) <- "seurat_clusters"
DefaultAssay(tumor) <- "Spatial8um"
genes_we_like <- setdiff(rownames(tumor), bad)
markers <- FindAllMarkers(tumor, only.pos = TRUE, features = genes_we_like,
                           min.pct = 0.1, logfc.threshold = 0.2, verbose = FALSE)
cat("OK:", nrow(markers), "marker rows returned, columns:", paste(colnames(markers), collapse=", "), "\n\n")

cat("== Merge + joint sketch-clustering (the un-official extension) ==\n")
# Merge the RAW per-sample objects (before their individual SCT/sketch
# processing above) so the joint SCTransform fit below is a single model
# across both conditions -- merging two objects that were each already
# SCTransform'd independently would combine two *different* SCT models
# under one assay name, which isn't statistically meaningful.
tumor_raw  <- make_fake_sample(condition = "tumor", seed = 1)
normal_raw <- make_fake_sample(condition = "normal_adjacent", seed = 2)
obj_list <- list(Tumor = tumor_raw, Normal = normal_raw)
spatial <- merge(x = obj_list[[1]], y = obj_list[2:length(obj_list)], add.cell.ids = names(obj_list))
DefaultAssay(spatial) <- "Spatial8um"
cat("Merged object:", ncol(spatial), "bins,", length(Layers(spatial, assay = "Spatial8um")), "Spatial8um layers (pre-join)\n")

# merge() keeps each original object's counts in its OWN layer
# ("counts.1", "counts.2", ...) rather than combining them -- a
# deliberate Seurat v5 design (cheap merges, join only when you need to).
# Left unjoined, SCTransform fits a SEPARATE model per layer instead of
# one pooled model across both conditions, and ProjectData then refuses
# to project against multiple reference SCT models. Join before sketching.
spatial <- JoinLayers(spatial, assay = "Spatial8um")
cat("Spatial8um layers after JoinLayers:", length(Layers(spatial, assay = "Spatial8um")), "(expect 1)\n")

spatial <- NormalizeData(spatial, assay = "Spatial8um", verbose = FALSE) %>%
  FindVariableFeatures(assay = "Spatial8um", verbose = FALSE)

n_sketch_joint <- min(1500, ncol(spatial))
n_before_joint_sketch <- ncol(spatial)
spatial <- SketchData(spatial, assay = "Spatial8um", ncells = n_sketch_joint,
                       method = "LeverageScore", sketched.assay = "sketch")
stopifnot(
  "Joint SketchData did not actually subsample" =
    ncol(spatial[["sketch"]]) < n_before_joint_sketch
)
DefaultAssay(spatial) <- "sketch"
spatial <- SCTransform(spatial, assay = "sketch", new.assay.name = "sketch", vars.to.regress = "percent.mt",
                        return.only.var.genes = FALSE, verbose = FALSE) %>%
  RunPCA(npcs = 10, verbose = FALSE) %>%
  FindNeighbors(dims = 1:10, k.param = 30, verbose = FALSE) %>%
  FindClusters(resolution = 0.6, verbose = FALSE) %>%
  RunUMAP(dims = 1:10, return.model = TRUE, verbose = FALSE)
spatial <- ProjectData(
  object = spatial, assay = "Spatial8um", sketched.assay = "sketch",
  full.reduction = "pca.full", sketched.reduction = "pca", umap.model = "umap",
  dims = 1:10, normalization.method = "SCT",
  refdata = list(seurat_clusters = "seurat_clusters")
)
DefaultAssay(spatial) <- "Spatial8um"
cat("OK. Joint clusters:", nlevels(spatial$seurat_clusters), "\n")
print(table(spatial$condition, spatial$seurat_clusters))

cat("\n== JoinLayers + renormalize ==\n")
spatial <- JoinLayers(spatial, assay = "Spatial8um")
spatial <- NormalizeData(spatial, assay = "Spatial8um", verbose = FALSE) %>%
  ScaleData(assay = "Spatial8um", features = rownames(spatial), verbose = FALSE)
cat("OK:", length(Layers(spatial, assay = "Spatial8um")), "layer(s) in Spatial8um after JoinLayers (expect 1)\n\n")

cat("== FindMarkers, tumor vs normal_adjacent (the DE step) ==\n")
DefaultAssay(spatial) <- "Spatial8um"
degs <- FindMarkers(spatial, ident.1 = "tumor", ident.2 = "normal_adjacent", group.by = "condition",
                     min.pct = 0.1, logfc.threshold = 0.1, assay = "Spatial8um", verbose = FALSE) %>%
  rownames_to_column("gene")
cat("OK:", nrow(degs), "DE rows, columns:", paste(colnames(degs), collapse=", "), "\n")
stopifnot(all(c("gene","avg_log2FC","p_val_adj") %in% colnames(degs)))

cat("\n== GSEA ranking construction (fgsea itself not installed here, just the dplyr chain) ==\n")
gsea_ranks <- degs %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene, avg_log2FC) %>% deframe()
cat("OK: ranking vector of length", length(gsea_ranks), "built from ALL", nrow(degs),
    "tested genes (not just p_val_adj<=0.05 subset)\n")
stopifnot(length(gsea_ranks) == nrow(degs))

cat("\n=== ALL CHECKS PASSED ===\n")
