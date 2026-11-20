# =============================================================================
## Plots heatmap for selected diapause gene list. Generates Fig 1j & ED Fig 8
## R 4.3.1
## Written by Dr. Noa Novershtern 
## https://github.com/hannalab/PRI_EXUT
# =============================================================================
library(Seurat)
library(pheatmap)
library(dplyr)
library(purrr)

WORKING_DIR  <- "<Enter path here>"
setwd(WORKING_DIR)
set.seed(42)

heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)
combined <- readRDS("data/combined_pei_bedzhov_yilmaz.rds")


## ---------------------------------------------------------------
## 1. Parse categorized gene list file
## ---------------------------------------------------------------
lines <- readLines("data/diapause_genes.list")
lines <- lines[nzchar(trimws(lines))]
lines

########## read file Alejandro_genes.list with multiple categories #######
gene_df <- purrr::map_dfr(lines, function(line) {
  parts <- strsplit(line, ":")[[1]]
  category <- trimws(parts[1])
  genes_raw <- trimws(parts[2])
  genes <- trimws(strsplit(genes_raw, ",")[[1]])
  data.frame(category = category, gene = genes, stringsAsFactors = FALSE)
})
gene_df$category <- factor(gene_df$category, levels = unique(gene_df$category))
gene_df <- gene_df %>% arrange(category)
gene_order_full <- gene_df$gene

gene_order_full

## ---------------------------------------------------------------
## 2. Subset by celltype -> average by sample ->
##    z-score per dataset -> clip -> drop zero-variance -> plot
## ---------------------------------------------------------------

  seurat_obj = combined
  gene_order   = gene_order_full
  gene_df      = gene_df
  celltype_use = "epiblast"
  title        = "Epiblast"

 celltype_col = "celltype2"
  sample_col   = "sample"
  dataset_col  = "dataset"

  sample_order <- c(
    "EDG75", "EDG95", "Rea", "E55",
    "E4.5", "E5.0", "E5.5",
    "PRI_E45_24h", "PRI_E45_48h", "PRI_E45_72h"
  )
  
  
    obj_sub <- subset(seurat_obj, subset = !!sym(celltype_col) == celltype_use)
  
  genes_use <- gene_order[gene_order %in% rownames(obj_sub@assays$RNA)]
  missing   <- setdiff(gene_order, genes_use)
  if (length(missing) > 0) {
    message("[", title, "] Missing ", length(missing), " genes: ",
            paste(missing, collapse = ", "), " out of ", length(gene_order), " genes.")
  }
  
  obj_sub <- JoinLayers(obj_sub, assay = "RNA")
  mat <- GetAssayData(obj_sub, assay = "RNA", layer = "data")[genes_use, , drop = FALSE]
  
  samples <- as.character(obj_sub@meta.data[[sample_col]])
  names(samples) <- colnames(obj_sub)
  samples <- samples[colnames(mat)]
  
  datasets <- as.character(obj_sub@meta.data[[dataset_col]])
  names(datasets) <- colnames(obj_sub)
  sample_to_dataset <- setNames(datasets, samples)[unique(samples)]
  
  # ---- average expression per sample ----
  mat_avg <- sapply(unique(samples), function(s) {
    cols <- names(samples)[samples == s]
    Matrix::rowMeans(mat[, cols, drop = FALSE])
  })
  rownames(mat_avg) <- rownames(mat)
  
  # ---- z-score per gene, SEPARATELY within each dataset's sample-columns ----
  mat_z <- mat_avg
  for (ds in unique(sample_to_dataset)) {
    cols_ds <- names(sample_to_dataset)[sample_to_dataset == ds]
    block <- mat_avg[, cols_ds, drop = FALSE]
    if (length(cols_ds) > 1) {
      mat_z[, cols_ds] <- t(scale(t(block)))
    } else {
      mat_z[, cols_ds] <- 0   # can't z-score a single sample
    }
  }
  mat_z[!is.finite(mat_z)] <- 0
  
  # ---- drop zero-variance genes (flat across all datasets) ----
  zero_var_genes <- rownames(mat_z)[apply(mat_z, 1, function(x) all(x == 0))]
  if (length(zero_var_genes) > 0) {
    message("[", title, "] Removing ", length(zero_var_genes), " zero-variance genes: ",
            paste(zero_var_genes, collapse = ", "))
  }
  mat_z <- mat_z[!rownames(mat_z) %in% zero_var_genes, , drop = FALSE]
  
  # ---- clip to +/-2 ----
  mat_z_clipped <- mat_z
  mat_z_clipped[mat_z_clipped >  2] <-  2
  mat_z_clipped[mat_z_clipped < -2] <- -2
  
  # ---- preserve gene order / category annotation (minus dropped genes) ----
  final_genes <- intersect(gene_order, rownames(mat_z_clipped))
  mat_z_clipped <- mat_z_clipped[final_genes, , drop = FALSE]
  
  gene_df_final <- gene_df %>% filter(gene %in% final_genes) %>%
    arrange(match(gene, final_genes))
  
  annotation_row <- data.frame(category = gene_df_final$category)
  rownames(annotation_row) <- gene_df_final$gene
  
  row_gaps <- cumsum(table(gene_df_final$category)[levels(gene_df_final$category)])
  row_gaps <- row_gaps[row_gaps > 0]
  row_gaps <- row_gaps[-length(row_gaps)]
  
  # ---- column order/annotation (group by dataset) ----
  col_order <- order(sample_to_dataset[colnames(mat_z_clipped)])
  mat_z_clipped <- mat_z_clipped[, col_order, drop = FALSE]
  ds_ordered <- sample_to_dataset[colnames(mat_z_clipped)]
  
  annotation_col <- data.frame(dataset = ds_ordered)
  rownames(annotation_col) <- colnames(mat_z_clipped)
  
  col_gaps <- cumsum(table(ds_ordered))
  col_gaps <- col_gaps[-length(col_gaps)]
  
  breaks <- seq(-2, 2, length.out = length(heat_colors) + 1)
  
  # ---- reorder matrix columns to match ----
  sample_order_present <- intersect(sample_order, colnames(mat_z_clipped))  # only samples actually in this celltype
  mat_z_clipped <- mat_z_clipped[, sample_order_present, drop = FALSE]
  
  ds_ordered <- sample_to_dataset[colnames(mat_z_clipped)]
  annotation_col <- data.frame(dataset = ds_ordered)
  rownames(annotation_col) <- colnames(mat_z_clipped)
  
  col_gaps <- cumsum(rle(as.character(ds_ordered))$lengths)
  col_gaps <- col_gaps[-length(col_gaps)]
  
  
  pheatmap(
    mat_z_clipped,
    cluster_rows   = FALSE,
    cluster_cols   = FALSE,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    gaps_col       = col_gaps,
    gaps_row       = row_gaps,
    show_colnames  = TRUE,
    show_rownames  = TRUE,
    fontsize_row   = 7,
    scale          = "none",
    color          = heat_colors,
    breaks         = breaks,
    border_color   = NA,
    main           = title
)
