# Exploratory selection of cluster 15-enriched genes for Fig. 1k-n
#
# Input:
#   data/processed/integrated_data.rds
#
# Purpose:
# - Annotate and inspect clusters in the integrated Flower/Silique dataset.
# - Identify genes enriched in cluster 15 using Seurat FindMarkers.
# - Compare cluster 15 with all other clusters and specifically with cluster 19.
# - Evaluate expression-frequency specificity across clusters.
# - Draw exploratory violin and dot plots used to select candidate genes.
#
# Provenance note:
# This is an exploratory source script, not the exact final plotting script for
# Fig. 1k-n. The final four genes (AT1G65970, AT3G01530, AT1G55670 and
# AT1G55330) are not all explicitly hard-coded together here. In particular,
# AT1G55670 and AT1G55330 do not occur by name in this file; they could have
# entered the displayed top-gene sets dynamically from the marker results.
# The script is retained to document the candidate-selection process.

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ---- データ読み込み ----
seurat_integrated <- readRDS(
  "results/figure1_cluster15_exploration/integrated_data.rds"
)

# ---- クラスターを確認 ----
Idents(seurat_integrated) <- "seurat_clusters"

cluster_ids <- levels(Idents(seurat_integrated))
cluster_ids

table(Idents(seurat_integrated))

# ---- 各クラスターのマーカー遺伝子を抽出 ----
DefaultAssay(seurat_integrated) <- "RNA"

markers_all <- FindAllMarkers(
  object = seurat_integrated,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25
)

dim(markers_all)
colnames(markers_all)
head(markers_all)

warnings()

Layers(seurat_integrated[["RNA"]])


seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])

Layers(seurat_integrated[["RNA"]])

table(seurat_integrated$seurat_clusters)

Idents(seurat_integrated) <- "seurat_clusters"

levels(Idents(seurat_integrated))

markers_all <- FindAllMarkers(
  seurat_integrated,
  assay = "RNA",
  only.pos = TRUE
)

dim(markers_all)

colnames(markers_all)

top20_markers <- markers_all %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()

dim(top20_markers)

DotPlot(
  seurat_integrated,
  features = c(
    "AT4G21750",  # ATML1：表皮
    "AT4G04890",  # PDF2：表皮
    "AT3G24140",  # FAMA：孔辺細胞
    "AT1G22690",  # GC1：孔辺細胞
    "AT5G07230",  # A9：タペート
    "AT2G16910",  # AMS：タペート
    "AT5G45880",  # LAT52：花粉
    "AT3G12110",  # ACT11：花粉
    "AT1G79430",  # APL：師部
    "AT1G22710",  # SUC2：師部
    "AT4G32880",  # ATHB8：維管束
    "AT1G69180",  # CRC：雌しべ・蜜腺
    "AT2G39060"   # SWEET9：蜜腺
  ),
  group.by = "seurat_clusters"
) +
  RotatedAxis()

DotPlot(
  seurat_integrated,
  features = c(
    "AT1G01280",  # CYP703
    "AT1G02050",  # PSKA
    "AT1G69500",  # CYP704B1
    "AT2G16910",  # AMS
    "AT3G11980",  # FAR2
    "AT4G34850",  # PSKB
    "AT4G35420",  # TKPR1
    "AT5G07230"   # A9
  ),
  group.by = "seurat_clusters"
) +
  RotatedAxis()


cluster_sample_table <- table(
  Cluster = seurat_integrated$seurat_clusters,
  Sample  = seurat_integrated$orig.ident
)

cluster_sample_percent <- round(
  prop.table(cluster_sample_table, margin = 1) * 100,
  1
)

cluster_sample_percent

write.csv(
  as.data.frame.matrix(cluster_sample_percent),
  "results/figure1_cluster15_exploration/cluster_sample_percent.csv"
)

nectary_markers <- c(
  "AT1G69180",  # CRC
  "AT2G39060",  # SWEET9
  "AT3G25810",  # TPS24
  "AT2G30650",
  "AT4G12530",
  "AT5G24270",
  "AT2G36190",
  "AT2G26580",  # YAB5
  "AT1G23300",
  "AT5G38120",
  "AT1G03770",
  "AT4G31820"   # MAB4
)

DotPlot(
  seurat_integrated,
  features = nectary_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

flower_markers <- c(
  # Pollen
  "AT1G19890",
  "AT1G64570",
  "AT5G16020",
  "AT5G55490",
  
  # Ovule
  "AT1G23420",
  
  # Petal
  "AT2G43680",
  "AT3G54340",
  "AT5G20240",
  "AT5G25980",
  "AT5G65590",
  
  # Stamen / Anther
  "AT2G21870",
  "AT3G12630",
  "AT3G47350",
  "AT3G60970",
  
  # Carpel / Stigma
  "AT2G16940",
  "AT3G04620",
  "AT4G18960",
  "AT5G49460",
  
  # Sepal
  "AT1G13930",
  "AT1G52510",
  "AT2G17033",
  "AT5G48620"
)

DotPlot(
  seurat_integrated,
  features = flower_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()


silique_markers <- c(
  # Seed coat
  "AT1G61720",  # BAN
  "AT3G59030",  # TT12
  "AT5G23260",  # TT16
  "AT1G79840",  # GL2
  
  # Embryo
  "AT1G21970",  # LEC1
  "AT3G24650",  # ABI3
  "AT5G59340",  # WOX2
  "AT5G45980",  # WOX8
  
  # Endosperm
  "AT5G60440",  # AGL62
  "AT2G35670",  # FIS2
  "AT1G02580",  # MEA
  "AT1G65330",  # PHE1
  
  # Fruit wall / valve
  "AT5G60910",  # FUL
  "AT3G58780",  # SHP1
  "AT2G42830",  # SHP2
  "AT4G00120"   # IND
)

DotPlot(
  seurat_integrated,
  features = silique_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()


surface_markers <- c(
  # Epidermis
  "AT4G21750",  # ATML1
  "AT4G04890",  # PDF2
  "AT1G79840",  # GL2
  "AT5G61590",  # TTG1
  
  # Guard cell
  "AT3G24140",  # FAMA
  "AT1G22690",  # GC1
  "AT1G12860",  # SCRM
  "AT3G26744",  # ICE1/SCRM2
  
  # Trichome
  "AT3G27920",  # GL1
  "AT1G64690",  # BLT
  "AT2G30432",  # TRY
  "AT2G46410"   # CPC
)

Rplot14_Surface <- DotPlot(
  seurat_integrated,
  features = surface_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

Rplot14_Surface


vascular_markers <- c(
  # Xylem
  "AT1G71930",  # VND7
  "AT4G35350",  # XCP1
  "AT5G17420",  # IRX3
  "AT5G54690",  # IRX8
  
  # Phloem / Companion cell
  "AT1G79430",  # APL
  "AT1G22710",  # SUC2
  "AT3G01680",  # SEOR1
  "AT3G12730",  # SEOR2
  
  # Procambium / Developing vasculature
  "AT4G32880",  # ATHB8
  "AT3G25710",  # TMO5
  "AT2G27230",  # LHW
  "AT5G61480"   # PXY
)

Rplot15_Vasculature <- DotPlot(
  seurat_integrated,
  features = vascular_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

Rplot15_Vasculature

meristem_markers <- c(
  # Stem cell / organizing center
  "AT2G17950",  # WUS
  "AT1G75820",  # CLV1
  "AT1G65380",  # CLV2
  "AT2G27250",  # CLV3
  
  # Meristem / organ initiation
  "AT1G62360",  # STM
  "AT5G03790",  # ATH1
  "AT5G61850",  # LFY
  "AT4G37750",  # ANT
  
  # Dividing cells
  "AT5G43077",  # CYCB1;1
  "AT1G20610",  # CYCB2;3
  "AT2G26760",  # CYCB1;4
  "AT2G28740"   # HISTONE H4
)

Rplot16_Meristem_Dividing <- DotPlot(
  seurat_integrated,
  features = meristem_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

Rplot16_Meristem_Dividing


photosynthetic_markers <- c(
  # Photosynthetic mesophyll-like cells
  "AT1G67090",  # RBCS1A
  "AT5G38420",  # RBCS2B
  "AT1G29910",  # CAB3
  "AT1G29920",  # LHCB1.3
  "AT3G54890",  # LHCA1
  "AT1G61520",  # LHCA3
  "ATCG00490",  # RBCL
  "ATCG00020",  # PSBA
  
  # Bundle sheath / vascular-associated ground tissue
  "AT3G54220",  # SCR
  "AT5G41920",  # SCL23
  "AT5G07700",  # MYB76
  "AT4G32880"   # ATHB8
)

Rplot17_Photosynthetic_GroundTissue <- DotPlot(
  seurat_integrated,
  features = photosynthetic_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

Rplot17_Photosynthetic_GroundTissue



fruit_markers <- c(
  # Valve / fruit wall
  "AT5G60910",  # FUL
  "AT1G69180",  # CRC
  "AT1G70510",  # KNAT2
  "AT1G23380",  # KNAT6
  
  # Valve margin / dehiscence zone
  "AT3G58780",  # SHP1
  "AT2G42830",  # SHP2
  "AT4G00120",  # IND
  "AT5G67110",  # ALC
  
  # Replum
  "AT5G02030",  # RPL
  "AT4G08150",  # BP/KNAT1
  "AT1G62990",  # KNAT7
  "AT1G65620"   # AS1
)

Rplot18_FruitWall_ValveMargin <- DotPlot(
  seurat_integrated,
  features = fruit_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

Rplot18_FruitWall_ValveMargin


seed_maturation_markers <- c(
  # Embryo maturation
  "AT3G24650",  # ABI3
  "AT3G26790",  # FUS3
  "AT1G21970",  # LEC1
  "AT1G28300",  # LEC2
  
  # Seed storage proteins
  "AT4G27140",  # SEED STORAGE ALBUMIN 1
  "AT4G28520",  # CRU3
  "AT5G44120",  # CRA1
  "AT1G03880",  # CRU2
  
  # Oleosin / lipid storage
  "AT4G25140",  # OLE1
  "AT5G40420",  # OLE2
  "AT3G01570",  # OLE4
  "AT2G25890"   # OLE5
)

Rplot19_SeedMaturation <- DotPlot(
  seurat_integrated,
  features = seed_maturation_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

Rplot19_SeedMaturation



seedcoat_subtype_markers <- c(
  # Mucilage secretory cells / outer seed coat
  "AT1G79840",  # GL2
  "AT1G53500",  # MUM4 / RHM2
  "AT5G63800",  # MUM2
  "AT2G47670",  # PMEI6
  
  # Seed coat differentiation
  "AT2G37260",  # TTG2
  "AT5G35550",  # TT2
  "AT5G41315",  # TTG1
  "AT3G50990",  # PER36
  
  # Endothelium / proanthocyanidin-producing cells
  "AT1G61720",  # BAN
  "AT3G59030",  # TT12
  "AT5G23260",  # TT16
  "AT1G17260"   # AHA10
)

Rplot20_SeedCoatSubtypes <- DotPlot(
  seurat_integrated,
  features = seedcoat_subtype_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

Rplot20_SeedCoatSubtypes


unresolved_clusters <- c(1, 2, 4, 5, 8, 10, 19, 20)

unresolved_markers <- FindAllMarkers(
  seurat_integrated,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25
)

top_unresolved_markers <- unresolved_markers %>%
  filter(cluster %in% unresolved_clusters) %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 15,
    with_ties = FALSE
  ) %>%
  ungroup()

write.csv(
  top_unresolved_markers,
  "results/figure1_cluster15_exploration/Rplot21_UnresolvedClusterMarkers.csv",
  row.names = FALSE
)

Rplot21_UnresolvedClusters <- DoHeatmap(
  seurat_integrated,
  features = unique(top_unresolved_markers$gene),
  group.by = "seurat_clusters",
  cells = WhichCells(
    seurat_integrated,
    idents = unresolved_clusters
  ),
  assay = "RNA"
) +
  NoLegend()

Rplot21_UnresolvedClusters


# 使用するクラスター情報を統一
Idents(seurat_integrated) <- "seurat_clusters"

unresolved_clusters <- c("1", "2", "4", "5", "8", "10", "19", "20")

# 各クラスターのマーカーを再計算
unresolved_markers <- FindAllMarkers(
  seurat_integrated,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.20,
  logfc.threshold = 0.25
)

top_unresolved_markers <- unresolved_markers %>%
  filter(cluster %in% unresolved_clusters) %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 15,
    with_ties = FALSE
  ) %>%
  ungroup()

# CSV出力
write.csv(
  top_unresolved_markers,
  "results/figure1_cluster15_exploration/Rplot21_UnresolvedClusterMarkers.csv",
  row.names = FALSE
)

# 未分類クラスターだけを抽出
unresolved_object <- subset(
  seurat_integrated,
  idents = unresolved_clusters
)

Rplot21_UnresolvedClusters <- DoHeatmap(
  unresolved_object,
  features = unique(top_unresolved_markers$gene),
  group.by = "seurat_clusters",
  assay = "RNA"https://c.p02.c4a.im/images/item/7239255/fa591c6fa806251e4f40333d745e9d3f27a02768aefde9f53e3ab49d8734f7d0?d=250x250
) +
  NoLegend()

Rplot21_UnresolvedClusters


write.csv(
  top_unresolved_markers,
  "results/figure1_cluster15_exploration/Rplot22_UnresolvedClusterMarkers.csv",
  row.names = FALSE
)

validation_markers <- c(
  # Microspore / developing pollen
  "AT5G22260",  # MS1
  "AT1G01280",  # CYP703A2
  "AT1G02050",  # LAP6
  "AT4G14080",  # MIKC*
  
  # Mature pollen / pollen tube
  "AT5G17480",  # PC1
  "AT2G04750",  # FIM5
  "AT5G41310",  # pollen-expressed kinesin
  "AT3G62230",  # ROH1
  
  # Chalazal endosperm
  "AT2G44240",  # chalazal cyst marker
  "AT4G13380",  # chalazal cyst marker
  "AT1G48910",  # NPF4.5
  "AT1G02580",  # MEA
  
  # Seed coat suberin / lipid transfer
  "AT2G48130",  # LTPG15
  "AT2G48140",  # EDA4
  "AT5G58860",  # CYP86A1
  "AT5G41040",  # ASFT
  
  # Xylem / lignified cells
  "AT4G37990",  # CAD8
  "AT4G35350",  # XCP1
  "AT5G17420",  # IRX3
  "AT5G54690"   # IRX8
)

Rplot23_UnresolvedValidation <- DotPlot(
  seurat_integrated,
  features = validation_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

Rplot23_UnresolvedValidation



nectary_core_markers <- c(
  "AT4G12530", # LTP family protein
  "AT3G25810", # TPS24
  "AT2G30650", # CoA-thioester hydrolase-related
  "AT1G69180", # CRC
  "AT2G39060", # SWEET9（旧annotation: MtN3 family）
  "AT5G24270", # calcium sensor-related
  "AT2G36190", # beta-fructosidase
  "AT1G77110", # auxin transport-related
  "AT2G26580", # YABBY5
  "AT5G44630", # terpene synthase
  "AT5G38120", # 4CL family
  "AT1G23300", # MATE transporter
  "AT2G22680", # RING-finger protein
  "AT1G70270",
  "AT5G60760",
  "AT4G34680",
  "AT1G03770",
  "AT1G23200", # pectinesterase
  "AT5G44620", # cytochrome P450
  "AT2G42830"  # AGL5/SHP2
)

# 実際にobject内に存在する遺伝子だけを使用
nectary_core_present <- intersect(
  nectary_core_markers,
  rownames(seurat_integrated)
)

length(nectary_core_present)


nectary_avg <- AverageExpression(
  seurat_integrated,
  assays = "RNA",
  features = nectary_core_present,
  group.by = "seurat_clusters",
  slot = "data"
)$RNA

# 遺伝子ごとにZ-score化
nectary_avg_z <- t(scale(t(as.matrix(nectary_avg))))
nectary_avg_z[nectary_avg_z > 2]  <- 2
nectary_avg_z[nectary_avg_z < -2] <- -2

pheatmap::pheatmap(
  nectary_avg_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  main = "Nectary-enriched gene expression",
  fontsize_row = 8
)



seurat_integrated <- AddModuleScore(
  seurat_integrated,
  features = list(nectary_core_present),
  assay = "RNA",
  name = "NectaryCore"
)

Rplot24_NectaryScore_Violin <- VlnPlot(
  seurat_integrated,
  features = "NectaryCore1",
  group.by = "seurat_clusters",
  pt.size = 0
)

Rplot24_NectaryScore_UMAP <- FeaturePlot(
  seurat_integrated,
  features = "NectaryCore1",
  reduction = "umap"
)

Rplot24_NectaryScore_Violin
Rplot24_NectaryScore_UMAP

nectary_score_summary <- seurat_integrated@meta.data %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(
    mean_nectary_score = mean(NectaryCore1, na.rm = TRUE),
    median_nectary_score = median(NectaryCore1, na.rm = TRUE),
    n_cells = dplyr::n()
  ) %>%
  dplyr::arrange(desc(mean_nectary_score))

nectary_score_summary


nectary_avg <- AverageExpression(
  seurat_integrated,
  assays = "RNA",
  features = nectary_core_present,
  group.by = "seurat_clusters",
  slot = "data"
)$RNA

nectary_avg_z <- t(scale(t(as.matrix(nectary_avg))))
nectary_avg_z[is.na(nectary_avg_z)] <- 0
nectary_avg_z[nectary_avg_z > 2] <- 2
nectary_avg_z[nectary_avg_z < -2] <- -2

Rplot25_NectaryHeatmap <- pheatmap::pheatmap(
  nectary_avg_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = NA,
  main = "Nectary-enriched genes across clusters",
  fontsize_row = 8
)

Rplot25_NectaryHeatmap


nectary_key_markers <- c(
  "AT1G69180", # CRC
  "AT2G39060", # SWEET9
  "AT3G25810", # TPS24
  "AT4G12530",
  "AT5G44630",
  "AT5G44620",
  "AT1G77110"
)

Rplot26_NectaryKeyMarkers <- DotPlot(
  seurat_integrated,
  features = nectary_key_markers,
  group.by = "seurat_clusters",
  assay = "RNA"
) +
  RotatedAxis()

Rplot26_NectaryKeyMarkers


Rplot27_FloralMYBs <- DotPlot(
  seurat_integrated,
  features = c(
    "AT3G27810", # MYB21
    "AT5G40350", # MYB24
    "AT3G01530"  # MYB57
  ),
  group.by = "seurat_clusters",
  assay = "RNA"
) +
  RotatedAxis()

Rplot27_FloralMYBs


#violin plotへ

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

# データの場所
outdir <- path.expand("results/figure1_cluster15_exploration")

# Seurat objectを読み込む
obj <- readRDS(file.path(outdir, "data/processed/integrated_data.rds"))


DefaultAssay(obj) <- "RNA"

# 分割されているRNA layersを結合
obj <- JoinLayers(
  object = obj,
  assay = "RNA"
)

Idents(obj) <- "seurat_clusters"

# 結合を確認
Layers(obj[["RNA"]])



markers15_all <- FindMarkers(
  object = obj,
  ident.1 = "15",
  ident.2 = NULL,
  assay = "RNA",
  test.use = "wilcox",
  only.pos = FALSE,
  min.pct = 0,
  logfc.threshold = 0
) |>
  tibble::rownames_to_column("gene") |>
  dplyr::mutate(
    pct_diff  = pct.1 - pct.2,
    pct_ratio = (pct.1 + 0.001) / (pct.2 + 0.001)
  ) |>
  dplyr::arrange(desc(avg_log2FC))


nrow(markers15_all)
head(markers15_all, 20)

markers15_all |>
  summarise(
    total       = n(),
    pct1_ge_001 = sum(pct.1 >= 0.01, na.rm = TRUE),
    pct1_ge_002 = sum(pct.1 >= 0.02, na.rm = TRUE),
    pct1_ge_005 = sum(pct.1 >= 0.05, na.rm = TRUE),
    pct1_ge_010 = sum(pct.1 >= 0.10, na.rm = TRUE)
  )

markers15_all |>
  filter(
    pct.1 >= 0.05,
    avg_log2FC > 0,
    !is.na(p_val_adj)
  ) |>
  summarise(
    total       = n(),
    pct2_le_001 = sum(pct.2 <= 0.01),
    pct2_le_002 = sum(pct.2 <= 0.02),
    pct2_le_005 = sum(pct.2 <= 0.05),
    diff_ge_005 = sum(pct_diff >= 0.05),
    diff_ge_010 = sum(pct_diff >= 0.10),
    diff_ge_020 = sum(pct_diff >= 0.20)
  )


markers15_candidates <- markers15_all |>
  filter(
    pct.1 >= 0.05,
    pct_diff >= 0.05,
    avg_log2FC > 0,
    !is.na(p_val_adj)
  ) |>
  arrange(
    desc(pct_diff),
    desc(avg_log2FC),
    desc(pct.1)
  )

write.csv(
  markers15_candidates,
  file.path(outdir, "Cluster15_candidates_pctdiff005.csv"),
  row.names = FALSE
)

nrow(markers15_candidates)
head(markers15_candidates, 30)





suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# 今表示した上位30遺伝子
genes_page <- head(markers15_candidates$gene, 30)

# RNA/data layerから発現量を取得
expr <- FetchData(
  obj,
  vars = genes_page,
  assay = "RNA",
  layer = "data"
)

df <- cbind(
  expr,
  seurat_clusters = Idents(obj)
)

df_long <- df |>
  mutate(
    seurat_clusters = factor(
      seurat_clusters,
      levels = as.character(0:21)
    )
  ) |>
  pivot_longer(
    cols = all_of(genes_page),
    names_to = "gene",
    values_to = "expr"
  ) |>
  drop_na(expr)

# DimPlotと同じクラスタ色を取得
p_dim <- DimPlot(
  obj,
  group.by = "seurat_clusters"
)

gb <- ggplot_build(p_dim)
scale_col <- gb$plot$scales$get_scales("colour")

levs_dim <- scale_col$get_limits()
cols_dim <- scale_col$palette(length(levs_dim))
names(cols_dim) <- levs_dim

# バイオリンプロット
p30 <- ggplot(
  df_long,
  aes(
    x = seurat_clusters,
    y = expr,
    fill = seurat_clusters
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    width = 1.2,
    adjust = 4,
    colour = NA
  ) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    width = 1.2,
    adjust = 4,
    fill = NA,
    colour = "black",
    linewidth = 0.25
  ) +
  facet_wrap(
    ~gene,
    ncol = 3,
    scales = "free_y"
  ) +
  scale_fill_manual(values = cols_dim) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    axis.text.x = element_text(
      size = 5,
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.title = element_text(size = 10)
  ) +
  xlab("Cluster") +
  ylab("Log-normalized expression")

p30


markers15_vs19 <- FindMarkers(
  object = obj,
  ident.1 = "15",
  ident.2 = "19",
  assay = "RNA",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0,
  logfc.threshold = 0
) |>
  tibble::rownames_to_column("gene") |>
  dplyr::mutate(
    pct_diff = pct.1 - pct.2
  ) |>
  dplyr::arrange(
    dplyr::desc(pct_diff),
    dplyr::desc(avg_log2FC)
  )

nrow(markers15_vs19)
head(markers15_vs19, 30)


markers15_combined <- markers15_candidates |>
  dplyr::rename(
    avg_log2FC_vs_all = avg_log2FC,
    pct15_vs_all      = pct.1,
    pct_other         = pct.2,
    pct_diff_vs_all   = pct_diff,
    p_val_adj_vs_all  = p_val_adj
  ) |>
  dplyr::left_join(
    markers15_vs19 |>
      dplyr::select(
        gene,
        avg_log2FC,
        pct.1,
        pct.2,
        pct_diff,
        p_val_adj
      ) |>
      dplyr::rename(
        avg_log2FC_vs19 = avg_log2FC,
        pct15_vs19      = pct.1,
        pct19           = pct.2,
        pct_diff_vs19   = pct_diff,
        p_val_adj_vs19  = p_val_adj
      ),
    by = "gene"
  )

nrow(markers15_combined)

head(
  markers15_combined |>
    dplyr::arrange(
      dplyr::desc(pct_diff_vs19),
      dplyr::desc(pct_diff_vs_all)
    ),
  30
)


library(Matrix)

candidate_genes <- markers15_combined$gene

# RNA/data layerを取得
rna_data <- LayerData(
  obj,
  assay = "RNA",
  layer = "data"
)

# 候補遺伝子のみ
rna_candidate <- rna_data[
  candidate_genes,
  ,
  drop = FALSE
]

cluster_id <- as.character(Idents(obj))
cluster_levels <- as.character(0:21)

# 各遺伝子について、各クラスターで発現する細胞の割合を計算
pct_by_cluster <- sapply(
  cluster_levels,
  function(cl) {
    Matrix::rowMeans(
      rna_candidate[
        ,
        cluster_id == cl,
        drop = FALSE
      ] > 0
    )
  }
)

colnames(pct_by_cluster) <- paste0("pct_cluster", cluster_levels)

pct_by_cluster <- data.frame(
  gene = rownames(pct_by_cluster),
  pct_by_cluster,
  row.names = NULL,
  check.names = FALSE
)

dim(pct_by_cluster)

head(pct_by_cluster[, 1:8], 10)



pct_specificity <- pct_by_cluster |>
  dplyr::rowwise() |>
  dplyr::mutate(
    # Cluster 15での検出率
    pct_cluster15_value = pct_cluster15,
    
    # Cluster 15以外で最も高い検出率
    max_pct_other = max(
      dplyr::c_across(
        -c(gene, pct_cluster15)
      ),
      na.rm = TRUE
    ),
    
    # Cluster 15と最大競合クラスターとの差
    pct_diff_vs_max_other =
      pct_cluster15_value - max_pct_other
  ) |>
  dplyr::ungroup()

# 最大競合クラスター名を追加
other_cluster_cols <- setdiff(
  paste0("pct_cluster", 0:21),
  "pct_cluster15"
)

pct_specificity$max_other_cluster <- apply(
  pct_specificity[, other_cluster_cols, drop = FALSE],
  1,
  function(x) {
    sub(
      "pct_cluster",
      "",
      other_cluster_cols[which.max(x)]
    )
  }
)

# Cluster 15との差が大きい順
pct_specificity <- pct_specificity |>
  dplyr::arrange(
    dplyr::desc(pct_diff_vs_max_other),
    dplyr::desc(pct_cluster15_value)
  )

head(
  pct_specificity |>
    dplyr::select(
      gene,
      pct_cluster15_value,
      max_other_cluster,
      max_pct_other,
      pct_diff_vs_max_other
    ),
  30
)



# Cluster 15以外の列を明示
other_cluster_cols <- setdiff(
  paste0("pct_cluster", 0:21),
  "pct_cluster15"
)

# Cluster 15以外で最大の検出率
other_pct_matrix <- as.matrix(
  pct_by_cluster[, other_cluster_cols, drop = FALSE]
)

max_col_index <- max.col(
  other_pct_matrix,
  ties.method = "first"
)

pct_specificity <- pct_by_cluster |>
  dplyr::mutate(
    pct_cluster15_value = pct_cluster15,
    max_other_cluster = sub(
      "pct_cluster",
      "",
      other_cluster_cols[max_col_index]
    ),
    max_pct_other = other_pct_matrix[
      cbind(seq_len(nrow(other_pct_matrix)), max_col_index)
    ],
    pct_diff_vs_max_other =
      pct_cluster15_value - max_pct_other
  ) |>
  dplyr::arrange(
    dplyr::desc(pct_diff_vs_max_other),
    dplyr::desc(pct_cluster15_value)
  )

head(
  pct_specificity |>
    dplyr::select(
      gene,
      pct_cluster15_value,
      max_other_cluster,
      max_pct_other,
      pct_diff_vs_max_other
    ),
  30
)



pct_specificity |>
  dplyr::summarise(
    total = dplyr::n(),
    cluster15_highest = sum(
      pct_diff_vs_max_other > 0,
      na.rm = TRUE
    ),
    diff_ge_002 = sum(
      pct_diff_vs_max_other >= 0.02,
      na.rm = TRUE
    ),
    diff_ge_005 = sum(
      pct_diff_vs_max_other >= 0.05,
      na.rm = TRUE
    ),
    diff_ge_010 = sum(
      pct_diff_vs_max_other >= 0.10,
      na.rm = TRUE
    )
  )




markers15_specific_82 |>
  dplyr::select(
    gene,
    pct_cluster15_value,
    max_other_cluster,
    max_pct_other,
    pct_diff_vs_max_other,
    avg_log2FC_vs_all,
    pct_diff_vs_all,
    avg_log2FC_vs19,
    pct_diff_vs19
  ) |>
  tibble::as_tibble() |>
  print(n = 82, width = Inf)


print(
  markers15_specific_82 |>
    dplyr::select(
      gene,
      pct_cluster15_value,
      max_other_cluster,
      max_pct_other,
      pct_diff_vs_max_other,
      avg_log2FC_vs_all,
      pct_diff_vs_all,
      avg_log2FC_vs19,
      pct_diff_vs19
    ),
  row.names = FALSE
)



genes82 <- markers15_specific_82$gene

p82 <- DotPlot(
  obj,
  features = genes82,
  group.by = "seurat_clusters",
  assay = "RNA",
  dot.scale = 5
) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 7
    ),
    axis.text.y = element_text(size = 6)
  ) +
  labs(
    x = "Cluster",
    y = "Gene"
  )

p82


ggsave(
  "Cluster15_candidate_82genes_DotPlot.pdf",
  plot = p82,
  width = 10,
  height = 18
)


genes20 <- markers15_specific_82 |>
  dplyr::slice_head(n = 20) |>
  dplyr::pull(gene)

p20 <- DotPlot(
  obj,
  features = genes20,
  group.by = "seurat_clusters",
  assay = "RNA",
  dot.scale = 7
) +
  RotatedAxis() +
  labs(
    x = "Gene",
    y = "Cluster"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 9
    ),
    axis.text.y = element_text(size = 9)
  )

p20


nectary_gene_list <- unique(toupper(nectary_gene_list))

nectary_genes_present <- intersect(
  nectary_gene_list,
  rownames(obj[["RNA"]])
)

length(nectary_gene_list)
length(nectary_genes_present)