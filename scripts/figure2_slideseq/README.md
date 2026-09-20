# Figure 2: Slide-seq spatial transcriptomics

## Vendor preprocessing

Raw paired-end FASTQ files were processed using the Curio Seeker Pipeline on
the Curio Bioscience Bioinformatics Cloud Portal powered by Latch. When
multiple sequencing runs were available for a library, FASTQ files
corresponding to the same read direction were concatenated. Read 1 was trimmed
to 55 nt using the FASTQ Sequence Trimmer workflow, while Read 1 and Read 2
remained paired. Each library was processed independently against the
Arabidopsis thaliana TAIR10 reference genome.

The vendor pipeline generated bead-level count matrices, matched-bead spatial
coordinates, quality-control metrics and quality-metric spatial maps, Seurat
objects and AnnData (h5ad) files. The quality-metric spatial maps shown in
Fig. 2b and Fig. 2l were generated directly by the Curio Seeker Pipeline and
were not produced using custom R code.
No custom code was used for the initial Curio Seeker processing. The analysis
procedure was based on the Curio Latch Portal User Guide (December 2024).
The exact pipeline version and execution identifiers should be transcribed
from each sample Report.html file when available.

Libraries analyzed:
- WT_3
- WT_5
- crc_2

Library-level sequencing and barcode-recovery metrics are reported in
Supplementary Data 5.

## Exploratory source scripts

- `exploratory/01_wt3_seurat_clustering_and_marker_plots.R`: WT_3 Seurat
  normalization, PCA, clustering, UMAP, spatial cluster mapping and marker-gene
  visualization underlying Fig. 2c, d and f-k.
- `exploratory/02_wt5_seurat_clustering_and_marker_plots.R`: the corresponding
  WT_5 workflow, used to inspect marker localization and the reproducibility of
  the clustering pattern in an independently prepared wild-type library.

Both wild-type scripts use 30 principal components for PCA, dimensions 1-10
for neighbor finding and UMAP, and a Seurat clustering resolution of 0.5.
These scripts are retained as exploratory, library-specific analyses rather
than presented as a single combined-sample workflow.

## Still required

- crc_2 clustering and spatial marker maps.
- WT-versus-crc differential-expression analysis for Fig. 2m.
- Gene-set overlap analysis for Fig. 2n.
- Final WT-versus-crc dot plot for Fig. 2o.
- crc spatial marker plots for Fig. 2p-s.
- A consolidated final script with explicit input and output filenames.
