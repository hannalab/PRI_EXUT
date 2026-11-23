# Ex utero development of mouse embryos from fertilization to mid-gestation

Analysis code for the single-cell RNA-seq data in:

> **Ex utero development of mouse embryos from fertilization to mid-gestation**, Gulben Gurhan, Alperen Yilmaz, Dmitry Lokshtanov, Mohammad Ghosheh, Francesco Roncato, Juana De La O, Lihua Wang, Xiushuai Yang, Castle Raley, Inbal Serfaty, Nizar Khatib, Reli Rachel Kakun, Nassrala Zawahri, Suhail Elias, Sergey Viukov, Itay Maza, Bernardo Oldak, Noa Novershtern,  Alejandro Aguilera-Castrejon,  Jacob H. Hanna, bioRxiv (2026).
> https://www.biorxiv.org/content/10.64898/2026.11.17.767265v1

*This README is temporary and will be updated on publication.*

## Requirements

R 4.3.1, Seurat v5, and the packages: Matrix, dplyr, tidyr, purrr, stringr, ggplot2, cowplot, patchwork, gridExtra, pheatmap, clustree, svglite, Cairo.

## Data

Raw data: GEO accession GSE349365.
External datasets: Bedzhov (Chen et al., https://doi.org/10.1016/j.stem.2024.06.015) and Pei (Cao et al., https://doi.org/10.1038/s41556-026-01971-3).

Edit the `<path>` placeholders and `WORKING_DIR` at the top of each script before running. Processed objects are read from / written to `data/`.

## Workflow

1. `read_*_dataset_clean.R` builds a Seurat object per time point from the Cell Ranger output.
2. `seurat_analysis_*_clean.R` produces the UMAP, dot plot and composition figures.
3. The remaining scripts produce the comparison figures and can run once the objects from steps 1–2 exist.

## Files

**Building the datasets**

| File | Description |
|---|---|
| `read_E45_dataset_clean.R` | Reads the E4.5 +24h/+48h/+72h 10x data, applies QC, integrates samples, clusters and saves the annotated Seurat object. |
| `read_E75_dataset_clean.R` | Reads the E7.5 10x data, applies QC, integrates samples, clusters and saves the Seurat object. |
| `read_E85_dataset_clean.R` | Reads the E8.5 10x data, applies QC, integrates samples, clusters and saves the Seurat object. |
| `read_E95_dataset_clean.R` | Reads the E9.5 10x data, applies QC, integrates samples, clusters and saves the Seurat object. |

**Analysis and figures**

| File | Description |
|---|---|
| `seurat_analysis_E45_datasets_clean.R` | Integrates our E4.5 data with the Chen et al and Cao et al datasets and annotates lineages (Fig. 1i, ED Fig. 7a). |
| `Sample_distance_heatmap_clean.R` | Heatmap of sample distances in PCA space for the integrated E4.5 datasets (ED Fig. 7b). |
| `Epiblast_progression_clean.R` | Compares epiblast module scores across three time points in utero (Cao et al) and ex utero (this study) (ED Fig. 7c). |
| `Epiblast_heatmap_clean.R` | Heatmap of selected diapause genes across samples of the integrated E4.5 datasets (Fig. 1j, ED Fig. 8). |
| `seurat_analysis_E75_clean.R` | E7.5 UMAPs, clustering-resolution tree, marker dot plot and cluster-composition bar plots. (Fig. 2e-f, ED Fig. 9a-b)|
| `seurat_analysis_E85_clean.R` | E8.5 UMAPs, clustering-resolution tree, marker dot plot and cluster-composition bar plots. (Fig. 3g-h, ED Fig. 12a-d)|
| `seurat_analysis_E95_clean.R` | E9.5 UMAPs, marker dot plot and cluster-composition bar plots. (Fig. 4e-f, ED Fig. 13a,c) |
| `E75_scatterplot_clean.R` | Per-cluster scatter plots and Pearson R of mean gene expression between PRI_ExUt and control at E7.5 (ED Fig. 10). |
| `E85_scatterplot_clean.R` | Per-cluster scatter plots and Pearson R of mean gene expression between PRI_ExUt and control at E8.5 (ED Fig. 10). |
| `E95_scatterplot_clean.R` | Per-cluster scatter plots and Pearson R of mean gene expression between PRI_ExUt and InUt at E9.5 (ED Fig. 10). |
| `Neuromesodermal_lineage_clean.R` | Subsets the combined all-time-point dataset and plot UMAP of Neuro-mesodermal cells only. (ED Fig. 13d-e) | 

**Helper files**

| File | Description |
|---|---|
| `E75_colors.R`, `E85_colors.R`, `E95_colors.R` | Cluster colour palettes for each time point. |
| `E75_order.list`, `E85_order.list`, `E95_order.list` | Display order of clusters in the dot plots (tab-separated: cluster, order). |

## Contact

Code written by Dr. Noa Novershtern, Hanna Lab. Questions: noa.novershtern@weizmann.ac.il 
