# =============================================================================
# E8.5 single-cell RNA-seq: downstream analysis and figures
#
# https://github.com/hannalab/PRI_ExUt
# Written by Dr. Noa Novershtern
# Tested with R 4.3.1 and Seurat v5.
# =============================================================================
WORKING_DIR    <- "<enter your working directory here>"
setwd(WORKING_DIR)

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(svglite)
library(clustree)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
set.seed(42)
source("E85_colors.R")

# ---- Load data --------------------------------------------------------------
input_rds  <- "combined_E85_dataset.rds"
E85 <- readRDS(input_rds)


# ---- UMAPs ------------------------------------------------------------------
p <- DimPlot(E85, reduction = "umap", group.by = "seurat_clusters",
             pt.size = 1, cols = my_colors) + NoAxes()

ggsave("E85_dataset_umap_combined.svg", p, width = 6, height = 6)

p <- DimPlot(E85, reduction = "umap", group.by = "seurat_clusters",
             split.by = "type", pt.size = 1, cols = my_colors) + NoAxes()
ggsave("E85_dataset_umap_split.svg", p, width = 12, height = 6)



# --- PRI UMAP vs control -----------------------

type_colors <- c(SEM = "brown1", TFSEM = "pink", PRI_ExUt = "brown3",
                 PSI_ExUt = "grey", InUt = "grey")

plot_highlight <- function(obj, highlight_types, colors, file) {
  # Keep the full-dataset UMAP range so the plot aligns with other panels
  umap <- Embeddings(obj, "umap")
  x_lim <- range(umap[, 1])
  y_lim <- range(umap[, 2])
  
  cells <- colnames(obj)[obj$type %in% highlight_types]
  p <- DimPlot(obj, reduction = "umap", cells = cells, group.by = "type",
               order = "PRI_ExUt", raster = FALSE, pt.size = 0.1, cols = colors) +
    coord_cartesian(xlim = x_lim, ylim = y_lim)
  ggsave(file, p, width = 6, height = 6)
}

plot_highlight(E85, c("PRI_ExUt", "InUt"), type_colors, "E85_PRI_vs_InUt.svg")
plot_highlight(E85, c("PRI_ExUt", "PSI_ExUt"),type_colors,  "E85_PRI_vs_PSI.svg")
            
# ---- QC and clustering resolution -------------------------------------------
p <- clustree(E85@meta.data, prefix = "integrated_snn_res.")
ggsave("E85_clustree.svg", p, device = svglite, width = 10, height = 10)

# ---- Marker dot plot (main figure) ------------------------------------------

cluster_mapping <- read.csv("E85.cluster_order.list", sep = "\t",
                            header = TRUE, stringsAsFactors = FALSE)
order_map <- setNames(cluster_mapping$order, cluster_mapping$seurat_cluster)
E85$order <- unname(order_map[as.character(E85$seurat_clusters)])


dotplot_groups <- list(
  endoderm         = c("Apoa1", "Apoa2", "Ttr", "Afp"),
  exe_endoderm     = c("Lama1", "Pdgfra", "Pth1r"),
  hindgut          = c("Sox17", "Klf5", "Cyp26a1"),
  midgut           = c("Nepn", "Shh", "Onecut2"),
  foregut          = c("Tbx1", "Otx2", "Trh"),
  blood            = c("Hba-a2", "Cited4", "Blvrb", "Gypa", "Hbb-bs"),
  hpc              = c("Cd33", "Ptprc", "Runx1"),
  endothelium      = c("Cd34", "Cdh5", "Kdr"),
  brain            = c("Pax6", "Six3", "Sox2"),
  neural           = c("Pax3", "Fgf15", "Tfap2a"),
  surface_ectoderm = c("Foxg1", "Trp63", "Grhl2", "Grhl3", "Pax7"),
  cardiac          = c("Myh6", "Ttn", "Actc1"),
  uncommitted_epc  = c("Gjb3", "Lgals1"),
  chorion_prog     = c("Cdx2", "Eomes", "Tfap2c", "Elf5"),
  intermediate_chorion = c("Ndrg1", "Irx2", "Irx3"),
  chorion          = c("Gcm1", "Ascl2", "Plet1"),
  exe_mesoderm     = c("Bmp4", "Bnc2", "Col5a2", "Lum", "Postn", "Hoxc8", "Tdo2"),
  mesoderm         = c("T", "Isl1"),
  allantois        = c("Tbx4", "Vcam1", "Rspo3", "Pitx1"),
  nmps             = c("Hoxa9", "Wnt3a", "Fgf8", "Gdf11"),
  somites          = c("Foxc2", "Tbx18", "Meox1"),
  notochord        = c("Noto"),
  pgc              = c("Dppa3", "Prdm1", "Pou5f1")
)
dotplot_features <- unlist(dotplot_groups, use.names = FALSE)
vline_positions  <- cumsum(lengths(dotplot_groups))[-length(dotplot_groups)] + 0.5

p <- DotPlot(E85, features = dotplot_features, group.by = "order",
             assay = "RNA", scale = TRUE, scale.min = 0) +
  scale_y_discrete(expand = expansion(mult = c(0.04, 0.04))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  geom_vline(xintercept = vline_positions, color = "grey40", linewidth = 0.1)
p
ggsave("E85_dataset_dotplot.svg", p, width = 18, height = 6)

# ---- Cluster composition ----------------------------------------------------
# Fraction of cells per cluster within each group (type) 
# Clusters absent from a group/sample are kept as zero.

cluster_order <- c("6", "10", "18", "2", "1", "20", "19", "12", "8", "15", "21",
                   "14", "23", "5", "4", "13", "9", "3", "11", "22", "7", "17",
                   "16", "0")
cluster_colors <- sapply(cluster_order, function(cl) get(paste0("cluster", cl, "_col")))


cluster_fractions <- function(obj, by) {
  obj@meta.data %>%
    mutate(cluster = Idents(obj)) %>%
    count(cluster, .data[[by]]) %>%
    complete(cluster, .data[[by]], fill = list(n = 0)) %>%
    group_by(.data[[by]]) %>%
    mutate(fraction = n / sum(n),
           label = ifelse(fraction >= 0.01,
                          paste0(round(fraction * 100, 1), "%"), "")) %>%
    ungroup() %>%
    mutate(cluster = factor(cluster, levels = cluster_order))
}

composition_plot <- function(df, by, labels = TRUE) {
  p <- ggplot(df, aes(x = .data[[by]], y = fraction, fill = cluster)) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = cluster_colors) +
    labs(x = NULL, y = "Fraction of cells", fill = "Cluster") +
    theme_classic()
  if (labels) {
    p <- p + geom_text(aes(label = label), position = position_fill(vjust = 0.5),
                       size = 4, color = "white")
  }
}

# By group 
df_type <- cluster_fractions(E85, "type")

ggsave("E85_dataset_barplot.svg", composition_plot(df_type, "type"),
       width = 6, height = 6)


# ---- PGC analysis (cluster 8) -----------------------------------------------
pgc_genes <- c("Pou5f1", "Nanog", "Klf4", "Sox2", "Dppa3", "Prdm1",
               "Bmp4", "Prdm14", "Alpl", "Nanos3")
cluster8 <- subset(E85, subset = seurat_clusters == 8)
cluster8$PGC_mean <- rowMeans(FetchData(cluster8, vars = pgc_genes))
print(table(cluster8$type))

plots <- FeaturePlot(cluster8, features = "PGC_mean", split.by = "type",
                     combine = FALSE, pt.size = 2, order = TRUE,
                     cols = c("lightblue", "red"))
plots <- lapply(plots, function(p) {
  p + theme(panel.border = element_blank(), axis.line = element_blank(),
            axis.ticks = element_blank())
})

ggsave("PGC_annot.svg", wrap_plots(plots), width = 18, height = 6)



sessionInfo()