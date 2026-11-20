# =============================================================================
## Read mouse peri-implantation scRNA-seq 10x data, 
## time point: E4.5+24h,+48h,+72h. Generate Seurat object
## R 4.3.1
## Written by Dr. Noa Novershtern 
## https://github.com/hannalab/PRI_EXUT
# =============================================================================
library(cowplot)
library(dplyr)
library(patchwork)
library(Seurat)
library(Matrix)
library(ggplot2)
library(stringr)
library(Cairo)
library(svglite)

set.seed(42)

# ── Paths ─────────────────────────────────────────────────────────────────────
WORKING_DIR  <- "<Enter your path here>"
DATA_DIR     <- "<path>/AGG_24_48_72/outs/count/filtered_feature_bc_matrix/"
OUT_RDS      <- file.path(WORKING_DIR, "data/E45_dataset.rds")

setwd(WORKING_DIR)


# ── Load 10x data ─────────────────────────────────────────────────────────────
# CellRanger aggr output; gem-group suffix encodes the source library (1-3)
T24_48_72.data   <- Read10X(data.dir = DATA_DIR)
T24_48_72.seurat <- CreateSeuratObject(
  counts      = T24_48_72.data,
  project     = "ExUt2",
  min.cells   = 3,
  min.features = 200
)

# ── Assign sample metadata ────────────────────────────────────────────────────
# gem-group index is embedded in the cell barcode as a trailing digit
groups <- as.integer(str_extract(colnames(T24_48_72.seurat), "\\d+$"))

table(groups)
# Capture condition label per cell
ct_map <- c(
  "1" = "PRI_E45_24h",
  "2" = "PRI_E45_48h",
  "3" = "PRI_E45_72h"
)

ct <- ct_map[as.character(groups)]
names(ct) <- colnames(T24_48_72.seurat)
T24_48_72.seurat$sample <- ct

cat("Cell counts per capture:\n"); print(table(T24_48_72.seurat$sample))

# ── QC filtering ──────────────────────────────────────────────────────────────
T24_48_72.seurat[["percent.mt"]] <- PercentageFeatureSet(T24_48_72.seurat, pattern = "^mt-")

# Keep cells with 500–8,000 detected genes and < 10 % mitochondrial reads
T24_48_72.sub <- subset(
  T24_48_72.seurat,
  subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 10
)
cat(sprintf("Cells after QC: %d (from %d)\n", ncol(T24_48_72.sub), ncol(T24_48_72.seurat)))

# ── Integration (CCA anchors) ─────────────────────────────────────────────────
dims <- 1:15
new_sub.list <- SplitObject(T24_48_72.sub, split.by = "sample")
new_sub.list <- lapply(new_sub.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


anchors  <- FindIntegrationAnchors(object.list = new_sub.list, dims = dims,
                                   anchor.features = 2000,reduction = "cca")
combined <- IntegrateData(anchorset = anchors, dims = dims,k.weight = 30)

# ── Dimensionality reduction & clustering ─────────────────────────────────────
DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = dims, min.dist = 0.1)
combined <- FindNeighbors(combined, reduction = "pca", dims = dims)

# Sweep resolutions 0–1.2 to allow downstream comparison
combined <- FindClusters(combined, resolution = c(0.1,0.2,0.3))

combined$seurat_clusters = combined$integrated_snn_res.0.2


# ── Quick visualization ───────────────────────────────────────────────────────
# Adjust my_colors to match the number of clusters at this resolution
my_colors <- scales::hue_pal()(length(unique(combined$integrated_snn_res.0.2)))

DimPlot(
  combined,
  reduction = "umap",
  group.by = "integrated_snn_res.0.2",
  label     = FALSE,
  raster    = FALSE
)

cells <- combined$orig.ident == "ExUt2"
combined$orig.ident[cells] <- combined$seurat_clusters[cells]

#---------- annotations QC ---------------------
Early_PrE = c("Gata4","Gata6","Foxa2","Pdgfra")
Visceral_endoderm = c("Apoa1","Apoa2","Ttr","Afp","Hnf4a","Apoe","Cystm1","Emb","Dkk1","Amot","Spink1")
Parietal_endoderm = c("Lama1","Pth1r","Lamb1","Dab2","Sparc","Snai1")
epiblast = c("Otx2","Sox2","Zic2","Pou5f1") 
Uncommited_EPC =c("Gjb3","Krt19","Lgals1") 
Chorion_prog = c("Cdx2","Eomes","Tfap2c","Elf5","Esrrb","Bmp4","Nat8l") 
morula_markers=c("Pdzd3","Npl","Rab29","Oasl2","Sult5a1","Dnah2os","Ifi27")

dotplot_features = c(Early_PrE,Visceral_endoderm,Parietal_endoderm,Uncommited_EPC,Chorion_prog,epiblast,morula_markers)

DotPlot(combined,features = dotplot_features,cluster.idents = F,assay = "RNA",scale =T,scale.min =0,group.by = "seurat_clusters") + scale_y_discrete(expand = expansion(mult = c(0.1, 0.1))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size =8)) 

######## save annotation to object #####
cluster_to_celltype <- c(
  "0" = "endoderm",
  "1" = "ectoderm",
  "2" = "ectoderm",
  "3" = "ectoderm",
  "4" = "epiblast"
)
combined$celltype <- unname(cluster_to_celltype[as.character(combined$seurat_clusters)])
table(combined$seurat_clusters, combined$celltype, useNA = "ifany")

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(combined, file = OUT_RDS)
cat("Saved:", OUT_RDS, "\n")

