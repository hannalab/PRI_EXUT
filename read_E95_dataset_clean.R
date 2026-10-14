## =============================================================================
## E9.5 embryo scRNA-seq: QC, CCA-anchor integration across samples, clustering
## R 4.3.1
## Written by Dr. Noa Novershtern 
## https://github.com/hannalab/PRI_EXUT
## =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(dplyr)
library(clustree)
library(svglite)
library(stringr)

set.seed(42)

# ── Paths ─────────────────────────────────────────────────────────────────────
WORKING_DIR <- "<enter your working directory here>"
DATA_DIR    <- "<path>/E75/outs/count/filtered_feature_bc_matrix/"
OUT_RDS      <- file.path(WORKING_DIR, "combined_E95_dataset.rds")

setwd(WORKING_DIR)

# ── Load 10x data ─────────────────────────────────────────────────────────────
# CellRanger aggr output; gem-group suffix encodes the source library (1-3)
E95.data   <- Read10X(data.dir = DATA_DIR)
E95.seurat <- CreateSeuratObject(
  counts      = E95.data,
  project     = "ExUt2",
  min.cells   = 3,
  min.features = 200
)

# ── Assign sample metadata ────────────────────────────────────────────────────
# gem-group index is embedded in the cell barcode as a trailing digit
groups <- as.integer(str_extract(colnames(E95.seurat), "\\d+$"))

# Capture condition label per cell
ct_map <- c(
  "1" = "InUt_E95_1e",
  "2" = "PRI_ExUt_E95_1e_d0_r1",
  "3" = "PRI_ExUt_E95_1e_d0_r2"
)

E95.seurat$sample <- unname(ct_map[as.character(groups)])

# Broad condition (in utero vs. ex utero)
type_map <- c("1" = "InUt", "2" = "PRI_ExUt", "3" = "PRI_ExUt")
E95.seurat$type <- unname(type_map[as.character(groups)])

cat("Cell counts per sample:\n"); print(table(E95.seurat$sample))
cat("Cell counts per type:\n");    print(table(E95.seurat$type))

# ── QC filtering ──────────────────────────────────────────────────────────────
E95.seurat[["percent.mt"]] <- PercentageFeatureSet(E95.seurat, pattern = "^mt-")

# Keep cells with 2,000–8,000 detected genes and < 5 % mitochondrial reads
E95.sub <- subset(
  E95.seurat,
  subset = nFeature_RNA > 2000 & nFeature_RNA < 8000 & percent.mt < 5
)
cat(sprintf("Cells after QC: %d (from %d)\n", ncol(E95.sub), ncol(E95.seurat)))

# ── Integration (CCA anchors) ─────────────────────────────────────────────────
dims <- 1:15

new_sub.list <- SplitObject(E95.sub, split.by = "sample")
new_sub.list <- lapply(new_sub.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

anchors  <- FindIntegrationAnchors(object.list = new_sub.list, dims = dims,
                                   anchor.features = 2000)
 combined <- IntegrateData(anchorset = anchors, dims = dims)

# ── Dimensionality reduction & clustering ─────────────────────────────────────
DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = dims, min.dist = 0.1)
combined <- FindNeighbors(combined, reduction = "pca", dims = dims)

# Sweep resolutions 0–1.2 to allow downstream comparison
combined <- FindClusters(combined, resolution = seq(0, 1.2, 0.1))

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(combined, file = OUT_RDS)
cat("Saved:", OUT_RDS, "\n")
