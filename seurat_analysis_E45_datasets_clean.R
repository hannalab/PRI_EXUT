# =============================================================================
# Integration of mouse peri-implantation scRNA-seq datasets
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

WORKING_DIR  <- "<Enter your path here>"
setwd(WORKING_DIR)
set.seed(42)

# ---- Load datasets ----
pei = readRDS("data/Pei_dataset.rds")
yilmaz = readRDS("data/E45_dataset.rds")
bedzhov = readRDS("data/Bedzhov_dataset.rds")
dims_use <- 1:5   # PCs used for UMAP and clustering

# ---- Use raw RNA counts; merge Seurat v5 layers into a single layer ----
DefaultAssay(pei) <- "RNA"
DefaultAssay(yilmaz) <- "RNA"
DefaultAssay(bedzhov) <- "RNA"

pei <- JoinLayers(pei, assay = "RNA")
yilmaz <- JoinLayers(yilmaz, assay = "RNA")
bedzhov <- JoinLayers(bedzhov, assay = "RNA")

# Label each cell with its source dataset
pei$dataset = "pei"
yilmaz$dataset = "yilmaz"
bedzhov$dataset = "bedzhov"

#---------------------------
# Put objects in a list
obj.list <- list(pei = pei, yilmaz = yilmaz,bedzhov = bedzhov)

# Normalize each dataset
obj.list <- lapply(obj.list, NormalizeData, verbose = FALSE)

# Find variable features
obj.list <- lapply(obj.list, FindVariableFeatures,
                   selection.method = "vst",
                   nfeatures = 2000,
                   verbose = FALSE)

# Select features for integration
features <- SelectIntegrationFeatures(
  object.list = obj.list,
  nfeatures = 2000
)

# Prepare objects for SCT-free integration
obj.list <- lapply(obj.list, ScaleData,features = features, verbose = FALSE)

obj.list <- lapply(obj.list, RunPCA, features = features, verbose = FALSE)

# Find anchors
integration_anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  anchor.features = features,
  reduction = "cca"      # use "cca" if datasets are very different
)

# Integrate
combined <- IntegrateData(
  anchorset = integration_anchors
)


# ---- Dimensionality reduction and clustering on the integrated assay ----
DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = dims_use,
                    return.model = TRUE, umap.method = "uwot")
combined <- FindNeighbors(combined, dims = dims_use)
combined <- FindClusters(combined, resolution = 0.1)


# ---- Build plotting table: UMAP coordinates + cell metadata ----
sc_df <- as.data.frame(Embeddings(combined, "umap"))
colnames(sc_df) <- c("UMAP_1","UMAP_2")
bg <- sc_df %>% dplyr::select(UMAP_1, UMAP_2)

sc_df$type <- "single-cell"
sc_df$sample <- combined$sample
sc_df$dataset <- combined$dataset
sc_df$celltype <- combined$celltype

# One colour per sample
my_colors = c("pink","purple","blue","brown3","orange","orange3","yellow","green","darkgreen","brown1")

# Figure 1-i 
ggplot(sc_df, aes(UMAP_1, UMAP_2)) +
  # grey background = every cell, repeated in each facet (no type_stage column)
  geom_point(data = bg, color = "grey90", size = 1) +
  # colored foreground = only the cells belonging to each facet
  geom_point(aes(color = sample), size = 1) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ dataset, ncol=3) +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

# Figure ED7-a
DimPlot(
  combined,
  reduction = "umap",
  group.by  = "seurat_clusters",
  label     = T,
  raster    = FALSE
)


# ---- Annotate clusters by lineage (based on marker expression) ----
combined$celltype2 <- case_when(
  combined$seurat_clusters %in% c("1", "6")           ~ "epiblast",
  combined$seurat_clusters %in% c("2", "3", "4", "8")  ~ "ectoderm",
  combined$seurat_clusters %in% c("0", "5", "7")       ~ "endoderm",
  TRUE ~ NA_character_
)

table(combined$celltype2, combined$seurat_clusters)  # sanity check

# ---- Save results ----
metadata <- combined@meta.data[c("sample","dataset","celltype2","seurat_clusters")]
saveRDS(combined,"data/combined_pei_bedzhov_yilmaz.rds")
write.csv(metadata,"data/pei_bedzhov_yilmaz_meta.csv")
