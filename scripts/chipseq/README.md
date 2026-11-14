# ChIP-seq analysis

This directory contains paired-end ChIP-seq workflows contributed for the CRC–MYB–MAB4 nectary analysis repository. The scripts are preserved as supplied; their commands and parameters have not been revised.

## Scripts

- `chipseq_pipeline_narrow.sh`: narrow-peak workflow, suitable for transcription factor ChIP-seq. It uses MACS2 2.2.6 with a matched Input BAM as control (`-c`), `-f BAMPE`, `-q 0.01`, and `--nomodel`.
- `chipseq_pipeline_broad.sh`: broad-region workflow, suitable for broad histone-mark enrichment. It uses SICER2 2.0 (`sicer`) with `-s tair10 -w 200 -rt 1 -f 150 -egf 0.74 -fdr 0.01 -g 600 -e 1000 -cpu 1`. As originally performed, SICER2 peak calling does **not** use the Input BAM; Input is included in the signal-profile analysis.

## Processing steps

1. Trimmomatic: paired-end adapter/quality trimming.
2. Bowtie2: alignment to a user-provided genome index.
3. SAMtools: BAM conversion, sorting and indexing.
4. Extraction of uniquely mapped reads, based on the original SAM filtering commands.
5. Picard MarkDuplicates: duplicate removal.
6. deepTools bamCoverage: RPKM-normalized BigWig tracks (5-bp bins).
7. MACS2 (narrow) or SICER2 (broad): enrichment region detection.
8. BEDTools: region merging; featureCounts: read counts over merged regions.
9. deepTools computeMatrix/plotHeatmap: scaled gene-body heatmap with 3-kb flanks.

## Inputs and setup

The scripts expect paired-end FASTQ files at `../rawdata/<sample>_1.fastq.gz` and `../rawdata/<sample>_2.fastq.gz`. Before running, set `INPUT`, `IP_SAMPLES`, `RAWDATA_DIR`, `TRIMMOMATIC_JAR`, `ADAPTER_FA`, `PICARD_JAR`, `BT2_INDEX`, and `GENE_BED` in the chosen script. For the narrow workflow, verify `GSIZE` for the reference genome. Run from a working directory with write permissions.

For narrow analysis, requirements are Trimmomatic 0.38, Bowtie2 2.3.4.2, SAMtools 1.9, Picard 2.18.11, deepTools 3.2.1, MACS2 2.2.6, BEDTools 2.27.1, and featureCounts 1.6.3. For broad analysis, substitute SICER2 2.0 for MACS2.

```bash
bash scripts/chipseq/chipseq_pipeline_narrow.sh
# OR:
bash scripts/chipseq/chipseq_pipeline_broad.sh
```

Run only the relevant workflow in its intended working directory. Outputs are written under `Trimmomatic/`, `Bowtie2_map/`, `SamtoBam/`, `Picard_dedup/`, `BigWig/`, the peak-calling directory, `featureCounts/`, and `ChIPseqfig/`.

## Reproducibility notes

These scripts document original analysis commands and are **not** validated as a portable, turnkey pipeline. In particular, the original `grep -v 'XS:'` alignment filter, SAM-to-BAM processing, duplicate handling, BED-to-SAF conversion, and paired-end counting options should be checked against the original data and environment before reuse. A histone modification does not always require broad-peak calling (for example, H3K4me3 may show narrow enrichment). Choose the analysis mode based on signal characteristics and experimental goals.

No raw sequencing data or large BAM/BigWig outputs are included here.
