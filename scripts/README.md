# Analysis scripts

The analysis scripts used in the manuscript will be organized here in execution order. Each final script will document its required inputs, software packages, parameters and generated outputs.

## Source scripts received

- `source/01_supplementary_figure_1_flower_snRNAseq.R`: exploratory and figure-generation code for the flower snRNA-seq analysis associated with Supplementary Figure 1, using the flower portion of NCBI GEO GSE226097.
- `source/02_supplementary_figure_1_silique_snRNAseq.R`: corresponding exploratory and figure-generation code for the silique snRNA-seq analysis associated with Supplementary Figure 1, using the silique portion of NCBI GEO GSE226097.
- `source/03_figure_1_integrated_umap_featureplots.R`: UMAP visualizations from a previously integrated flower/silique Seurat object. The script supports Fig. 1b, c, e, f and h–j. It does not contain the integration workflow, Fig. 1d heatmap code or Fig. 1g SWEET9 plot.
- `source/04_cluster_marker_heatmaps_and_tables.R`: per-cluster marker identification using Seurat `FindAllMarkers` (Wilcoxon test), selection of the top markers per cluster, cluster-marker heatmaps and output tables associated with Supplementary Data 1. This is not the Flower-versus-Silique DEG heatmap shown in Fig. 1d.
- `source/05_flower_silique_deg_and_cluster_marker_tables.R`: generates the cluster-specific marker results corresponding to Supplementary Data 1 and the top 100 Flower-versus-Silique DEG table corresponding to Supplementary Data 2. It supplies the gene list used for Fig. 1d but not the heatmap-drawing code. It contains no Supplementary Data 3 analysis.
- `source/06_legacy_flower_silique_integration_and_visualization.R`: earlier Seurat anchor-based integration of the Flower and Silique datasets, followed by PCA, UMAP, clustering and marker-gene visualization, including SWEET9 (AT2G39060). The obsolete combined-dataset pseudotime section from the supplied file was excluded because the final manuscript uses a Flower-only pseudotime analysis.
- `source/07_figure_1_flower_only_pseudotime_monocle3.R`: Flower-only Monocle3 trajectory analysis supporting Fig. 1q–t. It constructs a Monocle3 object from the public flower snRNA-seq Seurat object, learns a trajectory, assigns pseudotime from root node `Y_111`, projects CRC and AT1G65970 expression, and fits GAM-smoothed expression curves for the six genes shown in Fig. 1t. Fig. 1p cluster labels are inherited from the input object's metadata rather than explicitly selected in this script.

- `source/08_figure_1_cluster15_candidate_gene_exploration.R`: exploratory workflow used to identify cluster 15-enriched candidate genes and inspect them with violin and dot plots. This documents the selection process underlying Fig. 1k–n, but is not the exact final four-gene plotting script; AT1G55670 and AT1G55330 are not explicitly hard-coded in the source and may have entered dynamically through ranked marker lists.

- `source/09_figure_1_go_enrichment_dotplots.R`: plots previously generated agriGO v2.0 SEA results for clusters 0, 15 and 19; the cluster 15 plot is associated with Fig. 1o. This script visualizes imported results and does not itself perform GO enrichment. The supplied source maps raw P values, so an FDR column must be used explicitly if the final panel is described as showing FDR-adjusted P values.

The source scripts contain exploratory commands as well as code used for figure preparation. The final reproducible versions will include explicit input instructions, streamlined plotting steps and named output files after the relevant figure panels and object-preparation workflow are confirmed.
