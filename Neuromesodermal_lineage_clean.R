# =============================================================================
# Neuro-mesodermal lineage (ED Fig. d-e)
#
# Subsets the combined all-time-point dataset to epiblast, neural and
# mesodermal cell types, recomputes PCA/UMAP/clusters on the subset, assigns
# the condition label (InUt, SEM, PRI_ExUt, PSI_ExUt) from the sample name,
# and plots the UMAP coloured by cell type, faceted by condition x stage.
#
# R 4.3.1
# Written by Dr. Noa Novershtern 
# https://github.com/hannalab/PRI_EXUT
# =============================================================================
library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(data.table)
library(tidyr)


WORKING_DIR  <- "<Enter path here>"
setwd(WORKING_DIR)
set.seed(1234)

# ---- Load data ----
all.sub <- readRDS("combined_all_dataset.rds")

# Colours per cell type (alphabetical order of cell-type names)
neuro_meso_colors = c(
  "#FFCC99",#"Cardiomyocytes",
  "#D58700",#"Dorsal_Hindbrain_Prog",
  "#C2A14D",#"Dorsal_Midbrain_Neurons",
  "red4",#"Endotome",
  "yellow2", #Epiblast
  "coral3",#"Hindbrain_Floor_Plate",
  "#330000",#"Migratory_Neural_Crest",
  "#B983FF",#"Nascent_Meso",
  "#8B4513",#"Neural_Crest",
  "orange",#"Neural_Tube",
  "bisque3", #"Neuromesodermal_Prog",
  "#D8BFD8",#"Somites",
  "#C8A2B8",#"Somitic_Mesoderm",
  "#E76BF3"#"Tailbud_Mesoderm"
)


# ---- Subset to neural / mesodermal cell types ----
dims <- 1:15   # PCs used for UMAP and clustering
neuro_meso.sub = subset(all.sub,ct %in% c("Epiblast",
                                          "Cardiomyocytes",
                                          "Dorsal_Hindbrain_Prog",
                                          "Dorsal_Midbrain_Neurons",
                                          "Endotome",
                                          "Hindbrain_Floor_Plate",
                                          "Migratory_Neural_Crest",
                                          "Nascent_Meso",
                                          "Neural_Crest",
                                          "Neural_Tube",
                                          "Neuromesodermal_Prog",
                                          "Somites",
                                          "Somitic_Mesoderm",
                                          "Tailbud_Mesoderm"
))
# ---- Re-embed the subset ----
neuro_meso.sub <- RunPCA(neuro_meso.sub, npcs = 30, verbose = FALSE)
neuro_meso.sub <- RunUMAP(neuro_meso.sub, reduction = "pca", dims = dims, min.dist = 0.1)
neuro_meso.sub <- FindNeighbors(neuro_meso.sub, reduction = "pca", dims = dims)
neuro_meso.sub <- FindClusters(neuro_meso.sub, resolution = 0.2)

# ---- Assign condition (type) from sample name ----
s <- neuro_meso.sub$sample

table(neuro_meso.sub$type)
type <- dplyr::case_when(
  str_detect(s, "^InUt_")              ~ "InUt",
  str_detect(s, "^SEM_")               ~ "SEM",
  str_detect(s, "^TFSEM_")             ~ "SEM",
  str_detect(s, "^ExUt_.*_d6")         ~ "PSI_ExUt",
  str_detect(s, "^ExUt_.*_d0")         ~ "PRI_ExUt",
  TRUE                                 ~ NA_character_   # catch anything unmatched
)



neuro_meso.sub$type <- type
neuro_meso.sub$type_stage = paste(neuro_meso.sub$type,neuro_meso.sub$stage,sep="_")

# ---- UMAP coloured by cell type ----
DimPlot(neuro_meso.sub, reduction = "umap", group.by="ct",label=T,pt.size = 0.3,cols=neuro_meso_colors,label.size = 3,label.box = T) + NoAxes() 


# ---- UMAP faceted by condition x stage ----
# Pull UMAP coordinates + metadata into one data frame
umap_df <- as.data.frame(Embeddings(neuro_meso.sub, reduction = "umap"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$ct         <- neuro_meso.sub$ct
umap_df$type_stage <- neuro_meso.sub$type_stage

# Background layer: ALL cells, with type_stage dropped so it's shared across facets
bg <- umap_df %>% dplyr::select(UMAP_1, UMAP_2)

ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
  # grey background = every cell, repeated in each facet (no type_stage column)
  geom_point(data = bg, color = "grey90", size = 0.3) +
  # colored foreground = only the cells belonging to each facet
  geom_point(aes(color = ct), size = 0.3) +
  scale_color_manual(values = neuro_meso_colors) +
  facet_wrap(~ type_stage, ncol=3) +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))


saveRDS(neuro_meso.sub,"neuro_meso_sub.rds")
#neuro_meso.sub <- readRDS("neuro_meso_sub.rds")
