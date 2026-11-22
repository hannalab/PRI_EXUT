# =============================================================================
# E7.5: per-cluster comparison of mean gene expression, PRI_ExUt 
# vs InUt/PSI_ExUt (ED Fig 10)
#
# For each Seurat cluster (integrated_snn_res.0.5), the mean normalised
# expression of every gene is computed separately in PRI_ExUt and InUt/PSI_ExUt
# cells.Genes are plotted as PRI_ExUt vs InUt (one panel per cluster), and the
# Pearson correlation (R) between the two conditions is reported per cluster.
#
# R 4.3.1
# Written by Dr. Noa Novershtern 
# https://github.com/hannalab/PRI_EXUT
# =============================================================================
WORKING_DIR = "<Enter path here>" 
setwd(WORKING_DIR)
library(Seurat)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
source("E75_colors.R")

# ---- Load data; use resolution-0.5 clusters ----
E75.sub<-readRDS("data/combined_E75_dataset.rds")
E75.sub$seurat_clusters <- E75.sub$integrated_snn_res.0.5
Idents(E75.sub) <- E75.sub$integrated_snn_res.0.5
DimPlot(E75.sub, reduction = "umap", group.by="seurat_clusters",label=F,pt.size = 1,cols = my_colors) + NoAxes() 

# Join layers if integrated object
E75.sub <- JoinLayers(E75.sub, assay = "RNA")

# Pull normalized expression + metadata
# (cells x genes data frame of log-normalised RNA expression)
expr_df <- as.data.frame(t(as.matrix(
  LayerData(E75.sub, assay = "RNA", layer = "data")
)))


table(E75.sub$seurat_clusters)
expr_df$type <- E75.sub$type
expr_df$seurat_clusters <- E75.sub$seurat_clusters

# Filter to ExUt and InUt only
expr_df <- expr_df[expr_df$type %in% c("PRI_ExUt", "PSI_ExUt"), ]
#expr_df <- expr_df[expr_df$type %in% c("PRI_ExUt", "InUt"), ]
# Compute mean expression per gene, per ct, per type
genes <- setdiff(colnames(expr_df), c("seurat_clusters", "type"))

mean_expr <- expr_df %>%
  group_by(seurat_clusters, type) %>%
  summarise(across(all_of(genes), mean), .groups = "drop")

# Reshape to wide: one row per gene per cluster, columns = ExUt and InUt means
plot_data <- mean_expr %>%
  tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "mean_exp") %>%
  tidyr::pivot_wider(names_from = type, values_from = mean_exp)

# ---- Pearson correlation between conditions, per cluster ----
r_labels <- plot_data %>%
  group_by(seurat_clusters) %>%
  summarise(
      r = cor(PRI_ExUt, PSI_ExUt, use = "complete.obs"),
    #r = cor(PRI_ExUt, InUt, use = "complete.obs"),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("R = ", round(r, 3)))


# ---- Scatter plot: one panel per cluster, background shaded by R ----
# Get ct categories for faceting
#ggplot(plot_data, aes(x = PRI_ExUt, y = InUt)) +
ggplot(plot_data, aes(x = PRI_ExUt, y = PSI_ExUt)) +
  
  # background rectangle per facet, colored by R
  geom_rect(
    data = r_labels,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = r),
    inherit.aes = FALSE,
    alpha = 0.4
  ) +
  geom_point(alpha = 0.5, size = 0.5, color = "grey25", shape=20) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +   # y = x reference
  geom_text(
    data = r_labels,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.5,
    inherit.aes = FALSE,
    size = 3.5, fontface = "bold"
  ) +
  facet_wrap(~ seurat_clusters, scales = "free") +
  scale_fill_gradient(low = "white", high = "steelblue", name = "R",limits = c(0.9, 1)) +
  labs(
    x     = "Mean expression (PRI_ExUt)",
    y     = "Mean expression (PSI_ExUt)"
    #y     = "Mean expression (InUt)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold")
  )
