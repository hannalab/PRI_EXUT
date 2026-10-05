## =============================================================================
## E7.5 dataset: UMAPs, cluster resolution tree, marker dot plot, and
## cluster-composition bar plots.
## Written by Dr. Noa Novershtern
## https://github.com/hannalab/PRI_EXUT
## R 4.3.1
## =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(clustree)
library(svglite)
set.seed(42)

WORKING_DIR <- "<enter your working directory here>"
setwd(WORKING_DIR)

## Defines my_colors 
source("E75_colors.R")

## ---- Load data --------------------------------------------------------------
E75.sub <- readRDS("combined_E75_dataset.rds")

## Clustering at resolution 0.5 is used throughout
E75.sub$seurat_clusters <- E75.sub$integrated_snn_res.0.5
Idents(E75.sub) <- "integrated_snn_res.0.5"

## ---- UMAPs ------------------------------------------------------------------
p <- DimPlot(E75.sub, reduction = "umap", group.by = "seurat_clusters",
             pt.size = 1, cols = my_colors)
ggsave("E75_dataset_umap.svg", p, device = svglite, width = 6, height = 6)


p <- DimPlot(E75.sub, reduction = "umap", group.by = "integrated_snn_res.0.5",
             split.by = "type", pt.size = 1, cols = my_colors)

ggsave("E75_dataset_umap_split.svg", p, device = svglite, width = 6, height = 6)

## ---- Cluster stability across resolutions -----------------------------------
p <- clustree(E75.sub@meta.data, prefix = "integrated_snn_res.")
ggsave("E75_clustree.svg", p, device = svglite, width = 10, height = 10)


## ---- Marker dot plot --------------------------------------------------------
## Marker genes per cell type; comments give the matching cluster(s) at res 0.5
feature_groups <- list(
  Epiblast             = c("Otx2", "Sox2", "Zic2", "Zic3", "Pou5f1"),           
  Nascent_mesoderm     = c("Mixl1", "T", "Fgf8", "Mesp1"),                      
  Notochord            = c("Noto"),                                             
  Blood_prog           = c("Runx1", "Gata2", "Lmo2", "Tal1"),                   
  Erythroids           = c("Hbb-y", "Flt3", "Spi1"),                            
  Visceral_end         = c("Ttr", "Hnf4a", "Lrp2", "Afp", "Apoa1"),             
  Parietal_end         = c("Lama1", "Lamb1", "Dab2", "Sparc", "Snai1"),         
  DE                   = c("Cer1", "Sox17", "Hesx1", "Foxa1", "Lhx1", "Sox7", "Trh"), 
  Early_PrE            = c("Gata4", "Gata6", "Foxa2", "Pdgfra"),                
  Uncommitted_EPC      = c("Gjb3", "Krt19", "Lgals1"),                          
  Chorion_prog         = c("Cdx2", "Eomes", "Tfap2c", "Elf5", "Esrrb", "Bmp4", "Nat8l"),
  Chorion              = c("Gcm1", "Ascl2", "Plet1"),                           
  Intermediate_Chorion = c("Fgfr2", "Ndrg1", "Irx2", "Irx3"),                   
  TGC_prog             = c("Adm", "Fosl1", "Hand1", "Trpm5", "Maged2", "Prl5a1"),
  Amnion               = c("Postn", "Hoxc8", "Tdo2", "Vim"),                    
  SpT_Gly              = c("Dlx3", "Ncam1", "Pcdh12", "Tpbpa")                  
)
dotplot_features <- unlist(feature_groups, use.names = FALSE)

## Vertical separators between marker groups
vline_positions <- cumsum(lengths(feature_groups))[-length(feature_groups)] + 0.5

## Map clusters to their display order on the dot plot (order.list, tab-separated:
## columns seurat_cluster, order)
cluster_mapping <- read.delim("E75_order.list", stringsAsFactors = FALSE)
mapping_vector  <- setNames(cluster_mapping$order, cluster_mapping$seurat_cluster)
E75.sub$order   <- unname(mapping_vector[as.character(E75.sub$integrated_snn_res.0.5)])

p <- DotPlot(E75.sub, features = dotplot_features, group.by = "order",
             assay = "RNA", scale = TRUE, scale.min = 0) +
  scale_y_discrete(expand = expansion(mult = c(0.04, 0.04))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  geom_vline(xintercept = vline_positions, color = "grey40", linewidth = 0.1)
ggsave("E75_dataset_dotplot.svg", p, device = svglite, width = 18, height = 6)

## ---- Cluster composition per sample type ------------------------------------
## Clusters are stacked in this order (bottom to top within each bar)
cluster_order <- c("5", "9", "1", "8", "6", "14", "4", "11", "7",
                   "15", "13", "2", "10", "12", "3", "16", "0")

## Look up each cluster's color (clusterN_col, defined in E75_colors.R) in that order
my_colors_barplot <- setNames(
  unlist(mget(paste0("cluster", cluster_order, "_col"))),
  cluster_order
)

## Fraction of cells in each cluster within each sample type.
## Segments under 1% are left unlabeled.
df <- E75.sub@meta.data %>%
  mutate(cluster = factor(Idents(E75.sub), levels = cluster_order)) %>%
  count(cluster, type) %>%
  group_by(type) %>%
  mutate(fraction = n / sum(n),
         label = ifelse(fraction >= 0.01, paste0(round(fraction * 100, 1), "%"), "")) %>%
  ungroup()

## Stacked bar plot with percentage labels centered in each segment
p <- ggplot(df, aes(x = type, y = fraction, fill = cluster)) +
  geom_col(position = "fill") +
  geom_text(aes(label = label), position = position_fill(vjust = 0.5),
            size = 4, color = "white") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = my_colors_barplot) +
  labs(x = "Type", y = "Fraction of cells", fill = "Cluster") +
  theme_classic()

ggsave("E75_barplot.svg", p, device = svglite, width = 4, height = 6)



sessionInfo()


