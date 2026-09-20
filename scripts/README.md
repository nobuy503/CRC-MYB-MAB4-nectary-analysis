# Analysis scripts

The analysis scripts used in the manuscript will be organized here in execution order. Each final script will document its required inputs, software packages, parameters and generated outputs.

## Source scripts received

- `source/01_supplementary_figure_1_flower_snRNAseq.R`: exploratory and figure-generation code for the flower snRNA-seq analysis associated with Supplementary Figure 1, using the flower portion of NCBI GEO GSE226097.
- `source/02_supplementary_figure_1_silique_snRNAseq.R`: corresponding exploratory and figure-generation code for the silique snRNA-seq analysis associated with Supplementary Figure 1, using the silique portion of NCBI GEO GSE226097.
- `source/03_figure_1_integrated_umap_featureplots.R`: UMAP visualizations from a previously integrated flower/silique Seurat object. The script supports Fig. 1b, c, e, f and h–j. It does not contain the integration workflow, Fig. 1d heatmap code or Fig. 1g SWEET9 plot.
- `source/04_cluster_marker_heatmaps_and_tables.R`: per-cluster marker identification using Seurat `FindAllMarkers` (Wilcoxon test), selection of the top markers per cluster, cluster-marker heatmaps and output tables associated with Supplementary Data 1. This is not the Flower-versus-Silique DEG heatmap shown in Fig. 1d.

The source scripts contain exploratory commands as well as code used for figure preparation. The final reproducible versions will include explicit input instructions, streamlined plotting steps and named output files after the relevant figure panels and object-preparation workflow are confirmed.
