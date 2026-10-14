# =============================================================================
# E9.5 single-cell RNA-seq: figures
#
# Written by Dr. Noa Novershtern
# https://github.com/hannalab/Mouse_ExUt
#
# Tested with R 4.3.1 and Seurat v5.
# =============================================================================

WORKING_DIR    <- "<enter your working directory here>"
setwd(WORKING_DIR)

library(Seurat)
library(ggplot2)
library(dplyr)
library(svglite)

source("E95_colors.R")

# ---- Load data --------------------------------------------------------------
E95 <- readRDS("combined_E95_dataset.rds")

# Clusters at resolution 0.3 are used throughout
E95$seurat_clusters <- E95$integrated_snn_res.0.3
Idents(E95) <- E95$seurat_clusters

# Gene-level plots use the RNA assay; the integrated assay is only for UMAP
DefaultAssay(E95) <- "RNA"

# Plotting order of clusters (from cluster_order.list)
cluster_mapping <- read.csv("E95_order.list", sep = "\t",
                            header = TRUE, stringsAsFactors = FALSE)
order_map <- setNames(cluster_mapping$order, cluster_mapping$seurat_cluster)
E95$order <- unname(order_map[as.character(E95$seurat_clusters)])
stopifnot(!anyNA(E95$order))

# ---- UMAPs ------------------------------------------------------------------
p <- DimPlot(E95, reduction = "umap", group.by = "seurat_clusters",
             pt.size = 1, cols = my_colors) + NoAxes()

ggsave("E95_dataset_umap_combined.svg", p, width = 6, height = 6)

p <- DimPlot(E95, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
             pt.size = 1, split.by = "type", cols = my_colors) + NoAxes()

ggsave("E95_dataset_umap_split.svg", p, width = 12, height = 6)

p <- DimPlot(E95, reduction = "umap", group.by = "type", split.by = "type",
             raster = FALSE, cols = c("brown3", "black")) + NoAxes()

p <- DimPlot(E95, reduction = "umap", group.by = "type",
             raster = TRUE, cols = c("brown3", "black")) + NoAxes()

ggsave("E95_dataset_umap_by_type.svg", p, width = 12, height = 6)

# ---- Marker dot plot --------------------------------------------------------
# Marker sets are kept as a named list so the plotted order and the separator
# lines between groups stay in sync.
feature_groups <- list(
  endoderm        = c("Apoa1", "Apoa2", "Ttr", "Afp"),
  exe_endoderm    = c("Lama1", "Pdgfra", "Pth1r"),
  hindgut         = c("Sox17", "Klf5", "Cyp26a1"),
  hpc             = c("Cd33", "Ptprc", "Runx1"),
  mep             = c("Tspan32", "Gata1", "Snca", "Arhgdib", "Tal1"),
  blood           = c("Hba-a2", "Cited4", "Blvrb", "Gypa", "Hbb-bs"),
  endothelium     = c("Cd34", "Cdh5", "Kdr"),
  brain           = c("Pax6", "Six3", "Sox2"),
  neural          = c("Pax3", "Fgf15", "Otx2"),
  dorsal_hindbrain_progenitors = c("Msx3", "Fabp7"),
  dorsal_midbrain_neurons      = c("Btbd17", "Pou3f2", "Dll1", "Trp53i11"),
  hindbrain_floor_plate        = c("Foxa2", "Igfbp5"),
  optic_vesicle   = c("Wnt8b", "Cntnap2"),
  nmps            = c("Hoxa9", "Hoxc10", "Cdx4"),
  cardiac         = c("Myh6", "Ttn", "Actc1"),
  chorion         = c("Gcm1", "Ascl2", "Plet1"),
  intermediate_chorion = c("Ndrg1", "Irx2", "Irx3"),
  exe_mesoderm    = c("Bmp4", "Bnc2", "Col5a2", "Lum", "Postn", "Hoxc8", "Tdo2"),
  surface_ectoderm = c("Foxg1", "Trp63", "Grhl2", "Grhl3", "Pax7"),
  epidermis       = c("Fermt1", "Wnt4", "Krt7", "Npnt", "Wnt6", "Sfn"),
  endotome        = c("Prrx2", "Twist1", "Fst", "Foxs1", "Ebf1", "Cdh11"),
  sclerotome      = c("Tcf15", "Foxd1", "Arg1", "Col26a1")
)
dotplot_features <- unlist(feature_groups, use.names = FALSE)
vline_positions  <- cumsum(lengths(feature_groups))[-length(feature_groups)] + 0.5

p <- DotPlot(E95, features = dotplot_features, group.by = "order",
             assay = "RNA", scale = TRUE, scale.min = 0) +
  scale_y_discrete(expand = expansion(mult = c(0.04, 0.04))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  geom_vline(xintercept = vline_positions, color = "grey40", linewidth = 0.1)
ggsave("E95_dataset_dotplot.svg", p, width = 22, height = 6)


# ---- Cluster composition by group ------------------------------------------
# Fraction of cells per cluster within each group; clusters absent from a
# group are kept as zero.
cluster_order <- c("3", "2", "20", "5", "15", "8", "6", "10", "16", "12", "18",
                   "17", "4", "1", "7", "21", "13", "0", "14", "11", "19", "9")
cluster_colors <- sapply(cluster_order, function(cl) get(paste0("cluster", cl, "_col")))

df <- E95@meta.data %>%
  mutate(cluster = Idents(E95)) %>%
  count(cluster, type) %>%
  tidyr::complete(cluster, type, fill = list(n = 0)) %>%
  group_by(type) %>%
  mutate(fraction = n / sum(n),
         label = ifelse(fraction >= 0.01,
                        paste0(round(fraction * 100, 1), "%"), "")) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = cluster_order))

p <- ggplot(df, aes(x = type, y = fraction, fill = cluster)) +
  geom_col(position = "fill") +
  geom_text(aes(label = label), position = position_fill(vjust = 0.5),
            size = 4, color = "white") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = cluster_colors) +
  labs(x = NULL, y = "Fraction of cells", fill = "Cluster") +
  theme_classic()
ggsave("E95_dataset_barplot.svg", p, width = 6, height = 6)


sessionInfo()