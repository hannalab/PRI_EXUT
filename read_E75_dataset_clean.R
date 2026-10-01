## =============================================================================
## E7.5 embryo scRNA-seq: QC, CCA-anchor integration across samples, clustering
## R 4.3.1
## Written by Dr. Noa Novershtern 
## https://github.com/hannalab/PRI_EXUT
## =============================================================================

library(Seurat)
library(ggplot2)

set.seed(42)

## ---- Paths ------------------------------------------------------------------
WORKING_DIR <- "<enter your working directory here>"
## Where to put the processed files from GEO
DATA_DIR    <- "<path>/E75/outs/count/filtered_feature_bc_matrix/"
setwd(WORKING_DIR)

## ---- Load Cell Ranger aggr output -------------------------------------------
E75.data   <- Read10X(data.dir = DATA_DIR)
E75.seurat <- CreateSeuratObject(counts = E75.data, project = "ExUt2",
                                 min.cells = 3, min.features = 200)

## ---- Sample annotation ------------------------------------------------------
## Cell Ranger aggr appends a library index to each barcode (e.g. "ACGT...-3").
## The index follows the order of libraries in the aggr CSV.
group <- sub(".*-", "", colnames(E75.seurat))

sample_names <- c("1" = "TFSEM_E75_1e",
                  "2" = "InUt_E75_1e_r1",
                  "3" = "PsI_ExUt_E75_1e_latestreak",
                  "4" = "PsI_ExUt_E75_1e_headfold",
                  "5" = "PrI_ExUt_E75_2e")

sample_types <- c("1" = "TFSEM",
                  "2" = "InUt",
                  "3" = "PsI_ExUt",
                  "4" = "PsI_ExUt",
                  "5" = "PrI_ExUt")

E75.seurat$sample <- unname(sample_names[group])
E75.seurat$type   <- unname(sample_types[group])

## ---- Quality control --------------------------------------------------------
## Keep cells with 2,000-8,000 detected genes and <10% mitochondrial reads.
E75.seurat[["percent.mt"]] <- PercentageFeatureSet(E75.seurat, pattern = "^mt-")
E75.sub <- subset(E75.seurat,
                  subset = nFeature_RNA > 2000 & nFeature_RNA < 8000 & percent.mt < 10)

## ---- Integration across samples (Seurat anchor-based) -----------------------
dims <- 1:15

obj.list <- SplitObject(E75.sub, split.by = "sample")
obj.list <- lapply(obj.list, function(x) {
  x <- NormalizeData(x)
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

anchors  <- FindIntegrationAnchors(object.list = obj.list, dims = dims,
                                   anchor.features = 2000)
combined <- IntegrateData(anchorset = anchors, dims = dims)
DefaultAssay(combined) <- "integrated"

## ---- Dimensionality reduction and clustering --------------------------------
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = dims, min.dist = 0.1)
combined <- FindNeighbors(combined, reduction = "pca", dims = dims)
## Clusters are computed over a range of resolutions; resolution 0.5 is used below.
combined <- FindClusters(combined, resolution = seq(0, 1.2, 0.1))

## ---- Save -------------------------------------------------------------------
saveRDS(combined, file = "combined_E75_dataset.rds")

sessionInfo()