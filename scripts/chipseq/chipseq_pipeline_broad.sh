#!/usr/bin/env bash
# ChIP-seq 解析（ペアエンド / broad peak）
#
# 使用ツール:
#   Trimmomatic 0.38 / Bowtie2 2.3.4.2 / SAMtools 1.9 / Picard 2.18.11
#   deepTools 3.2.1 / SICER2 2.0 / BEDTools 2.27.1 / featureCounts 1.6.3

set -euo pipefail

#--- 設定 ---------------------------------------------------------------------
INPUT="input_sample"                  # Inputのサンプル名（手順9のプロファイル図で使用）
IP_SAMPLES=(ip_sample1 ip_sample2)    # IPサンプル名
RAWDATA_DIR="../rawdata"              # <サンプル名>_1.fastq.gz / _2.fastq.gz

TRIMMOMATIC_JAR="trimmomatic-0.38.jar"
ADAPTER_FA="adapter.fa"
PICARD_JAR="picard.jar"
BT2_INDEX="genome_index"              # Bowtie2インデックスのプレフィックス
GENE_BED="genes.bed"                  # 遺伝子BED

ALL_SAMPLES=("${INPUT}" "${IP_SAMPLES[@]}")

mkdir -p Trimmomatic Bowtie2_map SamtoBam Picard_dedup BigWig \
         SICER2_peak_call featureCounts ChIPseqfig

#--- 1〜6: 全サンプル共通の前処理 -----------------------------------------------
for s in "${ALL_SAMPLES[@]}"; do
  # 1. トリミング
  java -jar ${TRIMMOMATIC_JAR} PE -threads 2 \
    ${RAWDATA_DIR}/${s}_1.fastq.gz ${RAWDATA_DIR}/${s}_2.fastq.gz \
    ./Trimmomatic/${s}_R1.trimmed.fastq.gz ./Trimmomatic/${s}_R1_unpaired.fastq.gz \
    ./Trimmomatic/${s}_R2.trimmed.fastq.gz ./Trimmomatic/${s}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTER_FA}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36

  # 2. マッピング
  bowtie2 -x ${BT2_INDEX} \
    -1 ./Trimmomatic/${s}_R1.trimmed.fastq.gz \
    -2 ./Trimmomatic/${s}_R2.trimmed.fastq.gz \
    -S ./Bowtie2_map/${s}.sam

  # 3. SAM -> ソート済みBAM
  samtools view -bS ./Bowtie2_map/${s}.sam > ./SamtoBam/${s}.bam
  samtools sort ./SamtoBam/${s}.bam -o ./SamtoBam/${s}.sorted.bam
  samtools index ./SamtoBam/${s}.sorted.bam

  # 4. ユニークマッピングリードの抽出
  grep -v 'XS:' ./Bowtie2_map/${s}.sam \
    | awk '{if($1~/^@/ || $5>0) print $0}' > ./Bowtie2_map/${s}.uniq.sam
  samtools view -bS ./Bowtie2_map/${s}.uniq.sam > ./SamtoBam/${s}.uniq.bam
  samtools sort ./SamtoBam/${s}.uniq.bam -o ./SamtoBam/${s}.sorted_uniq.bam
  samtools index ./SamtoBam/${s}.sorted_uniq.bam

  # 5. 重複リード除去
  java -XX:ParallelGCThreads=4 -jar ${PICARD_JAR} MarkDuplicates \
    TMP_DIR=./Picard_dedup SORTING_COLLECTION_SIZE_RATIO=0.15 \
    ASSUME_SORTED=true REMOVE_DUPLICATES=true CREATE_INDEX=true \
    INPUT=./SamtoBam/${s}.sorted_uniq.bam \
    O=./Picard_dedup/${s}.sorted_dedup.bam \
    METRICS_FILE=./Picard_dedup/${s}_markduplicates_metrics.txt
  mv ./Picard_dedup/${s}.sorted_dedup.bai ./Picard_dedup/${s}.sorted_dedup.bam.bai

  # 6. BigWig 作成
  bamCoverage -b ./Picard_dedup/${s}.sorted_dedup.bam -o ./BigWig/${s}.bw \
    --exactScaling --normalizeUsing RPKM -bs 5
done

#--- 7. ピークコール（SICER2）------------
for ip in "${IP_SAMPLES[@]}"; do
  mkdir -p ./SICER2_peak_call/${ip}
  sicer -t ./Picard_dedup/${ip}.sorted_dedup.bam \
    -s tair10 -w 200 -rt 1 -f 150 -egf 0.74 -fdr 0.01 \
    -o ./SICER2_peak_call/${ip} -g 600 -e 1000 -cpu 1
  cp ./SICER2_peak_call/${ip}/${ip}.sorted_dedup-W200-G600.scoreisland \
     ./SICER2_peak_call/${ip}/${ip}.sorted_dedup-W200-G600.scoreisland.bed
done

#--- 8. ピークの統合とリードカウント --------------------------------------------
# 各IPサンプルのscoreisland.bedを連結
cat ./SICER2_peak_call/*/*.sorted_dedup-W200-G600.scoreisland.bed > all_broadPeak.bed

bedtools sort  -i all_broadPeak.bed > all_broadPeak.sorted.bed
bedtools merge -i all_broadPeak.sorted.bed > all_broadPeak.sorted_merged.bed
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}
     {print $1"_"$2+1"_"$3, $1, $2+1, $3, "."}' \
  all_broadPeak.sorted_merged.bed > broadPeak.sort.saf

for s in "${ALL_SAMPLES[@]}"; do
  featureCounts -F SAF -a broadPeak.sort.saf \
    -o ./featureCounts/${s}.peakcount_withcommand.txt \
    ./Picard_dedup/${s}.sorted_dedup.bam
  sed '1d' ./featureCounts/${s}.peakcount_withcommand.txt > ./featureCounts/${s}.peakcount.txt
done

#--- 9. 遺伝子周辺のシグナルプロファイル ----------------------------------------
BW_FILES=(./BigWig/${INPUT}.bw)
for ip in "${IP_SAMPLES[@]}"; do
  BW_FILES+=(./BigWig/${ip}.bw)
done

computeMatrix scale-regions \
  -S "${BW_FILES[@]}" \
  -R ${GENE_BED} \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  -o ./ChIPseqfig/scale-regions_matrix.mat.gz -p 4

plotHeatmap -m ./ChIPseqfig/scale-regions_matrix.mat.gz \
  -out ./ChIPseqfig/all_scale-regions_Heatmap.png --colorMap RdBu_r
