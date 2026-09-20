# Flower-versus-silique DEGs and cluster-marker tables
#
# This source script supports:
#   - Supplementary Data 1: cluster-specific marker genes identified with
#     Seurat FindAllMarkers using a Wilcoxon rank-sum test.
#   - Supplementary Data 2: the top 100 differentially expressed genes between
#     flower and silique datasets, identified with Seurat FindMarkers using a
#     Wilcoxon rank-sum test.
#
# The internal S1/S2 labels in the supplied script predate the final manuscript
# numbering and should not be used to infer the published file numbering.
# No analysis corresponding to Supplementary Data 3 is included.
#
# The script selects the genes used in the Fig. 1d heatmap but does not contain
# the heatmap-drawing code.
#
# P values reported as zero by Seurat because of numerical underflow are
# replaced only for tabular display with one tenth of the smallest non-zero
# value (or 1e-300 if all values are zero). This replacement does not affect
# DEG filtering or ranking.
#
# The original local Desktop directory was replaced with repository-relative
# paths. The processed input is expected at integrated_data.rds; output files
# are written relative to the working directory.

## ----------------------------------------------------------
## Flower vs Silique DEGs (Seurat::FindMarkers, Wilcoxon)
## ・p_val, p_val_adj が 0 になったものは
##   「表示用に」小さな値に置き換えて CSV に出力
## ----------------------------------------------------------

library(Seurat)
library(dplyr)
library(readr)

## 0. パス設定 ------------------------------------------------
base_dir <- "."   # 自分のパスに必要なら変更
integrated_rds <- file.path(base_dir, "integrated_data.rds")

## 1. オブジェクト読み込み ----------------------------------
integrated <- readRDS(integrated_rds)

## 2. RNA アッセイのレイヤーを確認して、複数あれば JoinLayers ----
Layers(integrated[["RNA"]])

if (length(Layers(integrated[["RNA"]])) > 1) {
  integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])
}

Layers(integrated[["RNA"]])  # "counts" など1つだけになっていればOK

## 3. 条件（Flower / Silique）を Idents に設定 ----------------
table(integrated$orig.ident)  # "Flower" / "Silique" があるか確認

Idents(integrated) <- integrated$orig.ident
table(Idents(integrated))

## 4. Flower vs Silique の差次的発現解析（Wilcoxon） ----------
deg_flower_vs_silique <- FindMarkers(
  object          = integrated,
  ident.1         = "Flower",
  ident.2         = "Silique",
  logfc.threshold = 0,      # しきい値ナシで全部テスト
  min.pct         = 0.1,    # どちらかの群で10%以上の細胞で発現
  test.use        = "wilcox"
)

deg_flower_vs_silique <- as.data.frame(deg_flower_vs_silique)
deg_flower_vs_silique$gene <- rownames(deg_flower_vs_silique)

## 5. p_val / p_val_adj が 0 のものを「小さな値」に置き換え ----
##   ・Seurat の出力で 0 になったものは、
##     観測された最小の非ゼロ p の 1/10 に置き換える（表示用）
##   ・すべて 0 の場合の保険として 1e-300 を使う

nonzero_p <- deg_flower_vs_silique$p_val[
  deg_flower_vs_silique$p_val > 0
]

if (length(nonzero_p) == 0) {
  eps_p <- 1e-300
} else {
  eps_p <- min(nonzero_p, na.rm = TRUE) / 10
}

nonzero_p_adj <- deg_flower_vs_silique$p_val_adj[
  deg_flower_vs_silique$p_val_adj > 0
]

if (length(nonzero_p_adj) == 0) {
  eps_p_adj <- 1e-300
} else {
  eps_p_adj <- min(nonzero_p_adj, na.rm = TRUE) / 10
}

deg_flower_vs_silique <- deg_flower_vs_silique %>%
  mutate(
    p_val     = ifelse(p_val     == 0, eps_p,     p_val),
    p_val_adj = ifelse(p_val_adj == 0, eps_p_adj, p_val_adj)
  )

## 6. Flower+Silique 全細胞に対する平均カウント（baseMean_raw） ----
counts <- GetAssayData(integrated, assay = "RNA", layer = "counts")

cells_use  <- colnames(integrated)[integrated$orig.ident %in% c("Flower", "Silique")]
counts_sub <- counts[, cells_use, drop = FALSE]

gene_mean <- rowMeans(counts_sub)
# gene_mean の名前は gene ID（AT○○…）になっているはず

## 7. FDR < 0.05 ＆ |avg_log2FC| 上位100個を抽出 ---------------
deg_top100 <- deg_flower_vs_silique %>%
  filter(!is.na(p_val_adj), p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 100) %>%
  mutate(
    baseMean_raw = gene_mean[gene]  # gene_mean から平均カウントを引っ張る
  ) %>%
  select(
    gene,
    avg_log2FC,
    baseMean_raw,
    p_val,
    p_val_adj,
    pct.1,
    pct.2,
    everything()
  )

## 8. CSV 出力 -------------------------------------------------
out_path <- file.path(
  base_dir,
  "Supplementary_Table_S1_top100_DEG_flower_vs_silique_SeuratWilcox.csv"
)

write_csv(deg_top100, out_path)
message("Done. Output: ", out_path)









## ----------------------------------------------------------
## (b) Table S2: クラスター特異的マーカー
##  - FindAllMarkers (Wilcoxon)
##  - p_val, p_val_adj の 0 を小さな値に置き換え（表示用）
##  - 列順：gene, avg_log2FC, baseMean_raw, p_val, p_val_adj, pct.1, pct.2, ...
## ----------------------------------------------------------

library(Seurat)
library(dplyr)
library(readr)

## 0. パス設定 ------------------------------------------------
base_dir <- "."
integrated_rds <- file.path(base_dir, "integrated_data.rds")

## すでに integrated があればこの行はスキップでOK
# integrated <- readRDS(integrated_rds)

## 1. クラスターIDを Idents に設定 ----------------------------
# メタデータの列名チェック（必要なら）
colnames(integrated@meta.data)[1:10]

# seurat_clusters を Idents に使う
Idents(integrated) <- integrated$seurat_clusters
table(Idents(integrated))  # クラスター番号の分布を確認

## 2. クラスター特異的マーカーを FindAllMarkers で取得 --------
markers_all <- FindAllMarkers(
  object          = integrated,
  only.pos        = TRUE,      # クラスターで高発現の遺伝子のみ
  min.pct         = 0.25,      # 少なくとも一方の群で 25% 以上の細胞が発現
  logfc.threshold = 0.25,      # |log2FC| >= 0.25
  test.use        = "wilcox"
)

markers_all <- as.data.frame(markers_all)

## gene 列の確認（なければ rownames から作る）----------------
if (!"gene" %in% colnames(markers_all)) {
  markers_all$gene <- rownames(markers_all)
}

## 3. p_val / p_val_adj が 0 のものを小さな値に置き換え --------
##    ・Seurat の出力で 0 になったものは、
##      観測された最小の非ゼロ p の 1/10 に置き換える（表示用）
##    ・すべて 0 の場合の保険として 1e-300 を使う

nonzero_p <- markers_all$p_val[markers_all$p_val > 0]
if (length(nonzero_p) == 0) {
  eps_p <- 1e-300
} else {
  eps_p <- min(nonzero_p, na.rm = TRUE) / 10
}

nonzero_p_adj <- markers_all$p_val_adj[markers_all$p_val_adj > 0]
if (length(nonzero_p_adj) == 0) {
  eps_p_adj <- 1e-300
} else {
  eps_p_adj <- min(nonzero_p_adj, na.rm = TRUE) / 10
}

markers_all <- markers_all %>%
  mutate(
    p_val     = ifelse(p_val     == 0, eps_p,     p_val),
    p_val_adj = ifelse(p_val_adj == 0, eps_p_adj, p_val_adj)
  )

## 4. 全細胞に対する平均カウント（baseMean_raw）を計算 ---------
##    （クラスター関係なく、integrated 内の全細胞で平均）
counts <- GetAssayData(integrated, assay = "RNA", layer = "counts")

cells_use  <- colnames(integrated)
counts_sub <- counts[, cells_use, drop = FALSE]

gene_mean <- rowMeans(counts_sub)  # 名前が gene ID になっているはず

## 5. 各クラスターで Top 50 マーカーを抽出 --------------------
markers_top50 <- markers_all %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%  # avg_log2FC 大きい順
  slice_head(n = 50) %>%
  ungroup() %>%
  mutate(
    baseMean_raw = gene_mean[gene]  # 平均カウントを追加
  ) %>%
  select(
    gene,
    avg_log2FC,
    baseMean_raw,
    p_val,
    p_val_adj,
    pct.1,
    pct.2,
    everything()
  )

## 6. Supplementary Table S2 として保存 ------------------------
out_path_s2 <- file.path(
  base_dir,
  "Supplementary_Table_S2_top50_markers_per_cluster.csv"
)

write_csv(markers_top50, out_path_s2)
message("Done. Output S2: ", out_path_s2)



## 各クラスター Top100 マーカーだけ作る ----
markers_top100 <- markers_all %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 100) %>%
  ungroup() %>%
  mutate(
    baseMean_raw = gene_mean[gene]
  ) %>%
  select(
    gene,
    avg_log2FC,
    baseMean_raw,
    p_val,
    p_val_adj,
    pct.1,
    pct.2,
    everything()
  )

## CSV 出力 ----
out_path_s2_100 <- file.path(
  base_dir,
  "Supplementary_Table_S2_top100_markers_per_cluster.csv"
)

write_csv(markers_top100, out_path_s2_100)
message("Done. Output S2 Top100: ", out_path_s2_100)




## ネクタリー用：Top100 ＆ p=0 を小さい値に置き換え ----

# まず markers_top100 全体で p_val / p_val_adj の 0 を処理
nonzero_p <- markers_top100$p_val[markers_top100$p_val > 0]
eps_p <- if (length(nonzero_p) == 0) 1e-300 else min(nonzero_p, na.rm = TRUE) / 10

nonzero_p_adj <- markers_top100$p_val_adj[markers_top100$p_val_adj > 0]
eps_p_adj <- if (length(nonzero_p_adj) == 0) 1e-300 else min(nonzero_p_adj, na.rm = TRUE) / 10

markers_top100_clean <- markers_top100 %>%
  mutate(
    p_val     = ifelse(p_val     == 0, eps_p,     p_val),
    p_val_adj = ifelse(p_val_adj == 0, eps_p_adj, p_val_adj),
    # baseMean_raw がまだ無ければ gene_mean から追加（あればそのまま）
    baseMean_raw = if ("baseMean_raw" %in% colnames(.)) baseMean_raw else gene_mean[gene]
  ) %>%
  select(
    gene,
    avg_log2FC,
    baseMean_raw,
    p_val,
    p_val_adj,
    pct.1,
    pct.2,
    everything()
  )

## ネクタリーと見なすクラスターIDを指定（必要なものだけ）----
nectary_ids <- c("0", "15", "19")  # 好きなIDに書き換え

for (cid in nectary_ids) {
  tmp <- markers_top100_clean %>% filter(cluster == cid)
  out_path <- file.path(
    base_dir,
    paste0("Nectary_cluster", cid, "_top100_markers.csv")
  )
  write_csv(tmp, out_path)
  message("Done. Nectary markers (cluster ", cid, "): ", out_path)
}
