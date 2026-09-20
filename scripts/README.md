# Analysis scripts

The analysis scripts used in the manuscript will be organized here in execution order. Each final script will document its required inputs, software packages, parameters and generated outputs.

## Source scripts received

- `source/01_supplementary_figure_1_flower_snRNAseq.R`: exploratory and figure-generation code for the flower snRNA-seq analysis associated with Supplementary Figure 1, using the flower portion of NCBI GEO GSE226097.

The source script currently assumes that the Seurat object `GSE226097_flower_230221` is already present in the R environment. The final reproducible version will include explicit input loading and output generation after the relevant figure panels and object-preparation workflow are confirmed.
