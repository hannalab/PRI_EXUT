# =============================================================================
## Compare ModuleScore() in 3 time points for pei's (in-utero) 
## and yilmaz (ex-utero) epiblast. Generates ED Fig 7c
## R 4.3.1
## Written by Dr. Noa Novershtern 
## https://github.com/hannalab/PRI_EXUT
# =============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(svglite)
library(patchwork)

WORKING_DIR <- "<Enter path here>"
setwd(WORKING_DIR)
set.seed(42)
dims_use <- 1:5

# ---- Settings ----
SCORE_MAX     <- 1.5   # scores are rescaled to 0 – SCORE_MAX
N_TOP_MARKERS <- 1000
DATASET_COLS  <- c(pei = "blue", yilmaz = "#2E8B57")
STAGE_LABELS  <- c("1 (E4.5 / 24h)", "2 (E5.0 / 48h)", "3 (E5.5 / 72h)")

# ---- Load data ----
combined <- readRDS("data/combined_pei_bedzhov_yilmaz.rds")
DefaultAssay(combined) <- "RNA"
combined <- JoinLayers(combined)

table(combined$sample)
combined$dataset <- case_when(
  combined$sample %in% c("E4.5", "E5.0", "E5.5")                      ~ "pei",
  combined$sample %in% c("E55", "EDG75", "EDG95", "Rea")              ~ "bedzhov",
  combined$sample %in% c("PRI_E45_24h", "PRI_E45_48h", "PRI_E45_72h") ~ "yilmaz"
)

table(combined$celltype2,combined$dataset)

# Visualize umap only with in-utero based (pei's) annotations
DimPlot(
  combined,
  reduction = "umap",
  group.by = "level3_celltypes",
  label     = FALSE,
  raster    = FALSE,
  pt.size   = 0.1,
  order=T
)

# =====================================================================
# Helper functions
# =====================================================================

# Min-max scale to 0 – SCORE_MAX (instead of 0 – 1)
rescale_score <- function(score, new_max = SCORE_MAX) {
  rng <- range(score, na.rm = TRUE)
  (score - rng[1]) / (rng[2] - rng[1]) * new_max
}

# Marker genes for given cluster(s) -> module score -> rescaled column
add_cluster_score <- function(obj, clusters, score_name) {
  markers <- FindMarkers(
    obj, ident.1 = clusters, group.by = "seurat_clusters",
    only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25
  )
  top_markers <- markers %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = N_TOP_MARKERS) %>%
    rownames()
  
  obj <- AddModuleScore(obj, features = list(top_markers), name = score_name, seed = 42)
  obj[[paste0(score_name, "_norm")]] <- rescale_score(obj[[paste0(score_name, "1")]][, 1])
  obj
}

# Map each dataset's 3 timepoints onto a shared ordinal stage (1/2/3)
# so pei E4.5/E5.0/E5.5 and yilmaz 24h/48h/72h overlap on the same x positions
get_stage_points <- function(obj, celltypes, score_col, label) {
  obj@meta.data %>%
    filter(dataset %in% c("pei", "yilmaz"), celltype2 %in% celltypes) %>%
    mutate(
      stage = case_when(
        sample %in% c("E4.5", "PRI_E45_24h") ~ 1,
        sample %in% c("E5.0", "PRI_E45_48h") ~ 2,
        sample %in% c("E5.5", "PRI_E45_72h") ~ 3,
        TRUE ~ NA_real_
      ),
      cell_type_plot = label,
      score = .data[[score_col]]
    ) %>%
    filter(!is.na(stage))
}

# Jittered cells + per-dataset mean line
plot_progression <- function(points, means, title, y_lab, line_width = 2, pt_alpha = 0.05) {
  ggplot() +
    geom_jitter(
      data = points,
      aes(x = stage, y = score, color = dataset),
      width = 0.15, height = 0, alpha = pt_alpha, size = 0.8, show.legend = FALSE
    ) +
    geom_line(
      data = means,
      aes(x = stage, y = mean_score, color = dataset, group = dataset),
      linewidth = line_width
    ) +
    geom_point(
      data = means,
      aes(x = stage, y = mean_score, color = dataset),
      size = 3
    ) +
    scale_color_manual(values = DATASET_COLS) +
    scale_x_continuous(breaks = 1:3, labels = STAGE_LABELS) +
    scale_y_continuous(breaks = seq(0, 1.4, by = 0.2)) +
    coord_cartesian(ylim = c(0, 1.4)) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

# =====================================================================
# Scores (0 – 1.2)
# =====================================================================
combined <- add_cluster_score(combined, "1",          "Cluster1Score")  # epiblast
summary(combined$Cluster1Score_norm)  # sanity check: should span 0 – 1.2

# =====================================================================
# Epiblast plot
# =====================================================================
epiblast_points <- get_stage_points(combined, "epiblast",           "Cluster1Score_norm", "epiblast")

group_means <- function(points) {
  points %>%
    group_by(dataset, stage, sample) %>%
    summarise(mean_score = mean(score), n_cells = n(), .groups = "drop")
}

epiblast_plot <- plot_progression(epiblast_points, group_means(epiblast_points),
                                  "epiblast", "E5.5 epiblast similarity score")

epiblast_plot


