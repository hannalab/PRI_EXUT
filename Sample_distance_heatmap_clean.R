# =============================================================================
## Sample distance in PCA space (ED Figure 7-b), can be run after 
## seurat_analysis_E45_datasets_clean.R
## Bedzhov (Chen et al, https://doi.org/10.1016/j.stem.2024.06.015)
## Pei (Cao et al, https://doi.org/10.1038/s41556-026-01971-3)
## Yilmaz (this paper)
## R 4.3.1
## Written by Dr. Noa Novershtern 
## https://github.com/hannalab/PRI_EXUT
# =============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(svglite)
library(pheatmap)
library(purrr)
library(gridExtra)
library(grid)

WORKING_DIR  <- "<Enter your path here>"
setwd(WORKING_DIR)
set.seed(42)

combined <- readRDS("data/combined_pei_bedzhov_yilmaz.rds")
combined_meta <- read.csv(
  "data/pei_yilmaz_bedzhov_meta.csv",
  row.names = 1,           # assumes first column is the cell barcode
  stringsAsFactors = FALSE
)

DimPlot(
  combined,
  reduction = "umap",
  group.by  = "seurat_clusters",
  label     = T,
  raster    = FALSE
)

combined$celltype2 <- case_when(
  combined$seurat_clusters %in% c("1", "6")           ~ "epiblast",
  combined$seurat_clusters %in% c("2", "3", "4", "8")  ~ "ectoderm",
  combined$seurat_clusters %in% c("0", "5", "7")       ~ "endoderm",
  TRUE ~ NA_character_
)

table(combined$celltype2, combined$seurat_clusters)  # sanity check

red_name <- "pca"     
n_dims   <- 1:5      

embed_full <- Embeddings(combined, reduction = red_name)[, n_dims]

# ---- Compute + plot sample distance matrix for one cell-type ----

# choose one cell type
#  ct = "epiblast"
#  ct = "ectoderm"
  ct = "endoderm"

  
  cells_ct <- colnames(combined)[combined$celltype2 == ct]
  embed_ct <- embed_full[cells_ct, , drop = FALSE]
  samples_ct <- combined$sample[cells_ct]
  
  # check sample sizes within this celltype (small samples -> noisy centroids)
  print(table(samples_ct))
  
  # centroid per sample
  centroids <- sapply(split(rownames(embed_ct), samples_ct), function(cells) {
    colMeans(embed_ct[cells, , drop = FALSE])
  })
  centroids <- t(centroids)   # samples x dims
  
  dist_matrix <- as.matrix(dist(centroids, method = "euclidean"))
  
  pheatmap(
    dist_matrix,
    main = paste0("Sample distance (", red_name, " space) — ", ct),
    color = colorRampPalette(c("firebrick", "white", "navy"))(100),
    display_numbers = T,
    breaks = seq(0, 30, length.out = 101),
    number_format = "%.1f",
    fontsize_number = 8
  )

