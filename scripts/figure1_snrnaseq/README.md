# Figure 1: public snRNA-seq analysis

This folder contains the public flower and silique single-nucleus RNA-seq
analyses used for Fig. 1 and the associated supplementary analyses. The public
dataset was obtained from NCBI GEO accession GSE226097.

## Final figure scripts

- `final/12_figure_1d_flower_silique_deg_heatmap.R`: Flower-versus-Silique
  Wilcoxon analysis and the top-100 DEG heatmap in Fig. 1d.
- `final/13_figure_1k_n_fixed_gene_violin_plots.R`: violin plots for
  At1g65970, MYB57, At1g55670 and At1g55330 in Fig. 1k-n.
- `final/14_figure_1p_t_flower_pseudotime_final.R`: Flower-only UMAP,
  Monocle3 trajectory, expression maps and pseudotime curves in Fig. 1p-t.

The final scripts expect processed objects under `data/processed/` and write
outputs to `results/figure1/`.

## Source and exploratory scripts

- `source/01_supplementary_figure_1_flower_snRNAseq.R`: flower snRNA-seq
  exploration and Supplementary Fig. 1 visualizations.
- `source/02_supplementary_figure_1_silique_snRNAseq.R`: silique snRNA-seq
  exploration and Supplementary Fig. 1 visualizations.
- `source/03_figure_1_integrated_umap_featureplots.R`: integrated
  Flower/Silique UMAP and marker FeaturePlots supporting Fig. 1b, c and e-j.
- `source/04_cluster_marker_heatmaps_and_tables.R`: cluster marker analysis
  associated with Supplementary Data 1.
- `source/05_flower_silique_deg_and_cluster_marker_tables.R`: cluster-marker
  and Flower-versus-Silique DEG tables associated with Supplementary Data 1-2.
- `source/06_legacy_flower_silique_integration_and_visualization.R`: earlier
  anchor-based integration and marker visualization.
- `source/07_figure_1_flower_only_pseudotime_monocle3.R`: original
  Flower-only Monocle3 analysis.
- `source/08_figure_1_cluster15_candidate_gene_exploration.R`: candidate-gene
  exploration underlying Fig. 1k-n.
- `source/09_figure_1_go_enrichment_dotplots.R`: visualization of imported
  agriGO v2.0 SEA results underlying Fig. 1o.
- `source/10_cluster_marker_heatmap_exploration.R`: alternative cluster-marker
  heatmaps.
- `source/11_cluster15_reclustering_and_heatmap_exploration.R`: exploratory
  reclustering and heatmaps for integrated cluster 15.
