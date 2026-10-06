## =============================================================================
## E8.5 embryo scRNA-seq: QC, CCA-anchor integration across samples, clustering
## R 4.3.1
## Written by Dr. Noa Novershtern 
## https://github.com/hannalab/PRI_EXUT
## =============================================================================

library(Seurat)
set.seed(42)

# ---- Paths (edit to your environment) ---------------------------------------
DATA_DIR <- "<path>/E85/outs/count/filtered_feature_bc_matrix/"
OUTPUT_FILE <- "combined_E85_dataset.rds"
source("E85_colors.R")

# ---- Load data and QC filtering ---------------------------------------------
counts <- Read10X(data.dir = DATA_DIR)
E85.seurat <- CreateSeuratObject(counts = counts, project = "E85",
                          min.cells = 3, min.features = 200)
E85.seurat[["percent.mt"]] <- PercentageFeatureSet(E85.seurat, pattern = "^mt-")

# Keep cells with 2,000-8,000 detected genes and <10% mitochondrial reads
E85.seurat <- subset(E85.seurat, subset = nFeature_RNA > 2000 & nFeature_RNA < 8000 &
                percent.mt < 10)
rm(counts)

# ---- Sample annotation ------------------------------------------------------
# Cell Ranger aggr appends a library index to each barcode (e.g. "ACGT...-3").
# Map that index to the sample name (ct) and the sample group (type).
sample_names <- c(
  "1"  = "InUt_E85_1e_r1",  "2"  = "InUt_E85_1e_r2",  "3"  = "InUt_E85_1e_r4",
  "4"  = "InUt_E85_2e_r1",  "5"  = "InUt_E85_1e_r5",  "6"  = "InUt_E85_3e_r1",
  "7"  = "PsI_ExUt_E85_1e",  "8"  = "PsI_ExUt_E85_3e",
  "9"  = "TFSEM_E85_5e_r1", "10" = "TFSEM_E85_5e_r2",
  "11" = "TFSEM_E85_1e_r1", "12" = "TFSEM_E85_5e_r3",
  "13" = "SEM_E85_2e_r1",   "14" = "SEM_E85_1e_r1",   "15" = "SEM_E85_2e_r2",
  "16" = "SEM_E85_2e_r3",   "17" = "SEM_E85_2e_r4",
  "18" = "PrI_ExUt_E85_1e"
)
sample_groups <- c(
  setNames(rep("InUt",    6), 1:6),
  setNames(rep("PSI_ExUt", 2), 7:8),
  setNames(rep("TFSEM",   4), 9:12),
  setNames(rep("SEM",     5), 13:17),
  "18" = "PRI_ExUt"
)

lib_index <- sub(".*-", "", colnames(E85.seurat))
E85.seurat$ct   <- unname(sample_names[lib_index])
E85.seurat$type <- unname(sample_groups[lib_index])
stopifnot(!anyNA(E85.seurat$ct), !anyNA(E85.seurat$type))
print(table(E85.seurat$ct))

# ---- Integration (Seurat anchor-based, one batch per sample) ----------------
dims <- 1:15

obj_list <- SplitObject(E85.seurat, split.by = "ct")
obj_list <- lapply(obj_list, function(x) {
  x <- NormalizeData(x)
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

 anchors  <- FindIntegrationAnchors(object.list = obj_list, dims = dims,
                                   anchor.features = 2000)
combined <- IntegrateData(anchorset = anchors, dims = dims)
DefaultAssay(combined) <- "integrated"

# ---- Dimension reduction and clustering --------------------------------
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = dims, min.dist = 0.1)
combined <- FindNeighbors(combined, reduction = "pca", dims = dims)
combined <- FindClusters(combined, resolution = seq(0, 1.2, 0.1))

saveRDS(combined, file = OUTPUT_FILE)


sessionInfo()