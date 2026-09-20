# Cluster-marker identification, heatmaps and output tables
#
# This script analyzes a previously integrated flower/silique Seurat object.
# It identifies positive markers for each Seurat cluster using Wilcoxon tests,
# selects the top markers per cluster and generates cluster-marker heatmaps and
# output tables. It is associated with the cluster-specific marker analysis
# reported in Supplementary Data 1.
#
# Important: this is not the script for Fig. 1d, which shows the top 100
# differentially expressed genes between the flower and silique datasets.
#
# Local Desktop paths in the supplied source were replaced with repository-
# relative data/processed/ and results/ paths. Analysis parameters were retained.

library(Seurat)
library(ggplot2)

obj <- readRDS("data/processed/integrated_data.rds")

# 由来（Flower/Silique）で色分け
obj$dataset <- obj$orig.ident
obj$dataset <- gsub("^Flower.*","Flower", obj$dataset)
obj$dataset <- gsub("^Silique.*","Silique", obj$dataset)

p1 <- DimPlot(obj, reduction = "umap", group.by = "dataset", pt.size = 0.4) +
  labs(title = "UMAP by dataset (Flower vs Silique)")

p2 <- DimPlot(obj, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 0.4) +
  NoLegend() + labs(title = "UMAP by clusters (res=0.5)")

print(p1); print(p2)



library(Seurat)
library(dplyr)

obj <- readRDS("data/processed/integrated_data.rds")

# 0) DimPlotと同じクラスタをIdentsに
stopifnot("seurat_clusters" %in% colnames(obj@meta.data))
Idents(obj) <- factor(obj$seurat_clusters,
                      levels = sort(unique(as.character(obj$seurat_clusters))))

# 1) RNAアッセイに切替
DefaultAssay(obj) <- "RNA"

# 2) layer状況を診断
cat("DefaultAssay:", DefaultAssay(obj), "\n")
print("RNA layers BEFORE join:"); print(Layers(obj[["RNA"]]))
print(str(obj[["RNA"]]@layers, max.level=1))

# 3) 必要なら正規化データを作る（data layer）
if (!"data" %in% Layers(obj[["RNA"]])) {
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
}

# 4) v5対策：Assayに対して JoinLayers を実行（これが肝）
if ("JoinLayers" %in% ls("package:SeuratObject")) {
  obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])
} else if ("JoinLayers" %in% ls("package:Seurat")) {
  # 万一こちらにエクスポートされている環境なら
  obj <- Seurat::JoinLayers(obj, assay = "RNA")
}

print("RNA layers AFTER join:"); print(Layers(obj[["RNA"]]))
# data が単一 layer として存在しているか再確認
stopifnot("data" %in% Layers(obj[["RNA"]]))

# 5) FindAllMarkers を再実行（軽量化オプションつき）
markers <- tryCatch({
  FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    assay = "RNA",
    test.use = "wilcox",
    layer = "data",                 # v5 で有効
    max.cells.per.ident = 4000,
    verbose = TRUE
  )
}, error = function(e) { message("FindAllMarkers error: ", conditionMessage(e)); NULL })

# 6) まだ空なら閾値を緩めて再試行
if (is.null(markers) || nrow(markers) == 0) {
  message("No markers found or NULL. Retrying with relaxed thresholds...")
  markers <- tryCatch({
    FindAllMarkers(
      obj,
      only.pos = TRUE,
      min.pct = 0.05,
      logfc.threshold = 0,
      assay = "RNA",
      test.use = "wilcox",
      layer = "data",
      max.cells.per.ident = 4000,
      verbose = TRUE
    )
  }, error = function(e) { message("Retry error: ", conditionMessage(e)); NULL })
}

# 7) それでも空なら：クラスタ毎にフォールバック
if (is.null(markers) || nrow(markers) == 0) {
  message("Fallback: per-cluster FindMarkers loop...")
  clusters <- levels(Idents(obj))
  markers_list <- lapply(clusters, function(cl) {
    m <- tryCatch({
      FindMarkers(
        obj, ident.1 = cl, ident.2 = NULL,
        min.pct = 0.1, logfc.threshold = 0.25,
        assay = "RNA", test.use = "wilcox", layer = "data"
      )
    }, error = function(e) NULL)
    if (is.null(m) || nrow(m) == 0) return(NULL)
    m$gene <- rownames(m); m$cluster <- cl; m
  })
  markers <- bind_rows(markers_list)
}

# 8) 列名標準化（v5では 'group' になることがある）
if (is.data.frame(markers) && nrow(markers) > 0) {
  cn <- colnames(markers)
  if (!"cluster" %in% cn) {
    if ("group" %in% cn) markers <- dplyr::rename(markers, cluster = group)
    if ("ident" %in% cn) markers <- dplyr::rename(markers, cluster = ident)
  }
}
stopifnot(is.data.frame(markers), nrow(markers) > 0, "cluster" %in% colnames(markers))

# 9) 上位100/クラスタ抽出
score_col <- dplyr::case_when(
  "avg_log2FC" %in% names(markers) ~ "avg_log2FC",
  "avg_log2fc" %in% names(markers) ~ "avg_log2fc",
  "avg_logFC"  %in% names(markers) ~ "avg_logFC",
  "log2FC"     %in% names(markers) ~ "log2FC",
  TRUE ~ NA_character_
)
stopifnot(!is.na(score_col))

top100_by_cluster <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = .data[[score_col]], n = 100, with_ties = FALSE) %>%
  ungroup()

features_for_heatmap <- unique(top100_by_cluster$gene)
length(features_for_heatmap)




library(Seurat)
library(dplyr)
library(ggplot2)

# 1) クラスタID固定（既にやっていればスキップ可）
Idents(obj) <- factor(obj$seurat_clusters,
                      levels = sort(unique(as.character(obj$seurat_clusters))))
DefaultAssay(obj) <- "RNA"

# 2) スケール行列を取得して“存在して非NA”な行列に絞る
S <- GetAssayData(obj, assay = "RNA", layer = "scale.data")
genes_ok <- intersect(features_for_heatmap, rownames(S))
cells_all <- colnames(obj)
# まずはクラスタごとに手動ダウンサンプル
set.seed(123)
N <- 300  # 重ければさらに小さく（200など）調整
cells_use <- unlist(lapply(levels(Idents(obj)), function(cl){
  v <- WhichCells(obj, idents = cl)
  if (length(v) > N) sample(v, N) else v
}), use.names = FALSE)
cells_ok <- intersect(cells_all, cells_use)

# 3) 完全NAの遺伝子をドロップ → その集合だけ再スケール
S_sub <- S[genes_ok, cells_ok, drop = FALSE]
keep_genes <- rowSums(!is.na(S_sub)) > 0
genes_ok <- rownames(S_sub)[keep_genes]

# 4) クリアして対象だけスケール（警告を消したい場合）
obj[["RNA"]]@layers$scale.data <- NULL
obj <- ScaleData(obj, features = genes_ok, verbose = FALSE)



png("results/heatmap_top100_markers_per_cluster_downsampled.png",
    width = 4200, height = 3200, res = 300)
p_hm <- DoHeatmap(
  obj,
  features = genes_ok,
  group.by = "seurat_clusters",   # ← DimPlotと一致
  assay = "RNA",
  slot = "scale.data",
  cells = cells_ok,
  raster = TRUE,
  draw.lines = FALSE,
  disp.min = -2,   # クリップ範囲を明示して警告を抑える
  disp.max =  2
) +
  labs(title = "Top100 marker genes per cluster (0–21), downsampled") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(face = "italic", size = 6))
print(p_hm)
dev.off()



library(dplyr)

# --- 0) 前提：markers がある（無ければ前の FindAllMarkers/フォールバックで作成）

# --- 1) 列名の正規化（v5だと group/avg_log2fc など）
m <- markers
if (!"cluster" %in% names(m)) {
  if ("group" %in% names(m)) m <- rename(m, cluster = group)
  if ("ident" %in% names(m)) m <- rename(m, cluster = ident)
}
# FC列を一本化
fc_col <- dplyr::case_when(
  "avg_log2FC" %in% names(m) ~ "avg_log2FC",
  "avg_log2fc" %in% names(m) ~ "avg_log2fc",
  "avg_logFC"  %in% names(m) ~ "avg_logFC",
  "log2FC"     %in% names(m) ~ "log2FC",
  TRUE ~ NA_character_
)
stopifnot(!is.na(fc_col))
m <- m %>%
  mutate(cluster = as.character(cluster)) %>%
  rename(logFC = !!fc_col) %>%
  # よく使う列だけ前に出す（存在すれば）
  relocate(cluster, gene, logFC, p_val_adj, pct.1, pct.2, .before = 1)

# --- 2) フルのマーカー表を書き出し（全クラスタ・しきい値ゆるめ後の全行）
out_all <- m %>%
  arrange(as.numeric(cluster), desc(logFC))
write.csv(out_all,
          file = "results/markers_all_clusters_full.csv",
          row.names = FALSE)
# 圧縮版（任意）
# write.csv(out_all, gzfile("results/markers_all_clusters_full.csv.gz"), row.names = FALSE)

# --- 3) 各クラスタ上位100（重複OK：各クラスタで100ずつ、合計≈ 22*100）
top100_by_cluster <- m %>%
  group_by(cluster) %>%
  slice_max(order_by = logFC, n = 100, with_ties = FALSE) %>%
  mutate(rank_in_cluster = row_number()) %>%
  ungroup() %>%
  arrange(as.numeric(cluster), rank_in_cluster)

write.csv(top100_by_cluster,
          file = "results/markers_top100_per_cluster.csv",
          row.names = FALSE)

# --- 4) 遺伝子名だけの一覧も（ヒートマップ用 features_for_heatmap と一致）
features_for_heatmap <- unique(top100_by_cluster$gene)
write.table(features_for_heatmap,
            file = "results/features_for_heatmap_top100_per_cluster.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# --- 5) 1クラスタ=1行で遺伝子名をカンマ連結した「見やすい表」も（査読メモ向け）
by_cluster_compact <- top100_by_cluster %>%
  group_by(cluster) %>%
  summarise(genes_csv = paste(gene, collapse = ","), .groups = "drop") %>%
  arrange(as.numeric(cluster))
write.csv(by_cluster_compact,
          file = "results/markers_top100_per_cluster_compact.csv",
          row.names = FALSE)

# --- 6) ついでに確認用の件数表
table_counts <- top100_by_cluster %>% count(cluster, name = "n_top")
print(table_counts)  # すべて 100 になっていればOK
write.csv(table_counts,
          file="results/markers_top100_counts_per_cluster.csv",
          row.names = FALSE)

# install.packages("openxlsx") が必要
library(openxlsx)
wb <- createWorkbook()
# 全体シート
addWorksheet(wb, "ALL_full")
writeData(wb, "ALL_full", out_all)
# 上位100まとめシート
addWorksheet(wb, "top100_all")
writeData(wb, "top100_all", top100_by_cluster)
# クラスタ別シート
for (cl in sort(unique(top100_by_cluster$cluster))) {
  addWorksheet(wb, paste0("cl", cl))
  writeData(wb, paste0("cl", cl),
            top100_by_cluster %>% filter(cluster == cl))
}
saveWorkbook(wb, "results/markers_top100_per_cluster.xlsx", overwrite = TRUE)



library(Seurat)
library(dplyr)
library(ggplot2)

# 前提: obj, markers があり、Idents(obj) は seurat_clusters と一致
Idents(obj) <- factor(obj$seurat_clusters,
                      levels = sort(unique(as.character(obj$seurat_clusters))))
DefaultAssay(obj) <- "RNA"

# fold-change 列を特定
score_col <- c("avg_log2FC","avg_log2fc","avg_logFC","log2FC")
score_col <- score_col[score_col %in% names(markers)][1]

# 各クラスタ上位30（重複遺伝子は1回だけ）
top30 <- markers %>%
  mutate(cluster = as.character(ifelse("cluster" %in% names(.), cluster,
                                       ifelse("group" %in% names(.), group, ident)))) %>%
  group_by(cluster) %>%
  slice_max(order_by = .data[[score_col]], n = 30, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(as.numeric(cluster), desc(.data[[score_col]])) %>%
  distinct(gene, .keep_all = TRUE)

features30 <- top30$gene

# 各クラスタ 最大200細胞にダウンサンプル
set.seed(123)
N <- 200
cells_use <- unlist(lapply(levels(Idents(obj)), function(cl){
  v <- WhichCells(obj, idents = cl)
  if (length(v) > N) sample(v, N) else v
}), use.names = FALSE)

# 対象遺伝子だけクリーンに再スケール
obj[["RNA"]]@layers$scale.data <- NULL
obj <- ScaleData(obj, features = features30, verbose = FALSE)

# 列順：クラスタ順のまま
cells_ordered <- unlist(lapply(levels(Idents(obj)), function(cl){
  intersect(WhichCells(obj, idents = cl), cells_use)
}), use.names = FALSE)

# 出力（PNG推奨）
png("results/heatmap_top30_markers_per_cluster_downsampled.png",
    width = 3600, height = 2800, res = 300)
p_hm30 <- DoHeatmap(
  obj,
  features = features30,
  group.by = "seurat_clusters",
  assay = "RNA", slot = "scale.data",
  cells = cells_ordered,
  raster = TRUE, draw.lines = FALSE,
  disp.min = -1.8, disp.max = 1.8
) +
  labs(title = "Top 30 marker genes per cluster (0–21), downsampled") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(face = "italic", size = 6))
print(p_hm30)
dev.off()



#ブロックっぽいヒートマップ
library(Seurat)
library(dplyr)
library(pheatmap)

# 前提: obj, top100_by_cluster があり、Idents(obj) は seurat_clusters と一致
stopifnot("seurat_clusters" %in% colnames(obj@meta.data))
cl_levels <- sort(unique(as.character(obj$seurat_clusters)))
Idents(obj) <- factor(as.character(obj$seurat_clusters), levels = cl_levels)
DefaultAssay(obj) <- "RNA"

features100 <- unique(top100_by_cluster$gene)
score_col <- "logFC"   # あなたの列名

# 遺伝子のoriginとrank（FC降順）を決める
gene_origin <- top100_by_cluster %>%
  mutate(cluster = as.character(cluster)) %>%
  group_by(cluster) %>%
  arrange(desc(.data[[score_col]]), .by_group = TRUE) %>%
  mutate(rank_in_origin = row_number()) %>%
  ungroup() %>%
  group_by(gene) %>% slice(1) %>% ungroup() %>%
  select(gene, origin = cluster, rank_in_origin)

# 平均発現（明示的に group.by = "seurat_clusters"）
avg <- AverageExpression(
  obj, features = features100, assays = "RNA", slot = "data",
  group.by = "seurat_clusters"
)
mat <- as.matrix(avg$RNA)  # genes x groups (列名にプレフィックスが付くことあり)

# 列名の数値部分（末尾の数字）を抽出して0–21の順に並べ替え
col_raw <- colnames(mat)
cl_num  <- suppressWarnings(as.integer(sub(".*?(\\d+)\\s*$", "\\1", col_raw)))
# 存在チェック
stopifnot(all(sort(unique(cl_num)) %in% as.integer(cl_levels)))
ord <- order(cl_num, na.last = NA)
mat <- mat[, ord, drop = FALSE]
cl_num <- cl_num[ord]

# 行方向Zスコア
mat_z <- t(scale(t(mat))); mat_z[is.na(mat_z)] <- 0
colnames(mat_z) <- as.character(cl_num)

# 行の並び（origin→rank）。gene_origin と照合して並び替え＆ブロック境界を算出
go <- gene_origin[match(rownames(mat_z), gene_origin$gene), ]
go$origin_num <- as.integer(as.character(go$origin))
row_ord <- order(go$origin_num, go$rank_in_origin, na.last = TRUE)
mat_z_ord <- mat_z[row_ord, , drop = FALSE]

origin_ord <- go$origin_num[row_ord]
gap_idx <- cumsum(table(factor(origin_ord, levels = sort(unique(origin_ord)))))
gap_idx <- gap_idx[gap_idx < nrow(mat_z_ord)]

# 平均ヒートマップ（ブロック感を確認）
pheatmap(
  mat_z_ord,
  filename = "results/heatmap_top100_clusterAverages_ordered.pdf",
  cluster_rows = FALSE, cluster_cols = FALSE,
  gaps_row = gap_idx,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 6,
  color = colorRampPalette(c("#313695","#f7f7f7","#a50026"))(201)
)


library(Seurat)
library(ggplot2)

features_ordered <- rownames(mat_z_ord)   # 平均版で確定した行順（ブロック順）

# クラスタ→遺伝子集合（Top100）を作成（originに基づく）
genes_by_cluster <- lapply(sort(unique(go$origin_num)), function(k){
  features <- gene_origin$gene[as.integer(gene_origin$origin) == k]
  unique(features)
})
names(genes_by_cluster) <- as.character(sort(unique(go$origin_num)))

# モジュールスコアを付与（順に MS1..MS22 になる）
obj <- AddModuleScore(obj, features = genes_by_cluster, name = "MS", assay = "RNA")
ms_cols <- grep("^MS", colnames(obj@meta.data), value = TRUE)
stopifnot(length(ms_cols) >= length(genes_by_cluster))

# 各クラスタ内で該当MSスコア降順に列順を決める
cells_ordered <- unlist(lapply(seq_along(genes_by_cluster), function(i){
  cl <- names(genes_by_cluster)[i]
  ms <- ms_cols[i]
  v  <- WhichCells(obj, idents = cl)
  v[order(obj@meta.data[v, ms], decreasing = TRUE)]
}), use.names = FALSE)

# 重いので各クラスタ先頭N細胞だけに制限（可視性重視）
N <- 150   # きつければ 120/100 に下げる
cells_ordered <- unlist(lapply(names(genes_by_cluster), function(cl){
  vv <- cells_ordered[cells_ordered %in% WhichCells(obj, idents = cl)]
  if (length(vv) > N) vv[seq_len(N)] else vv
}), use.names = FALSE)

# 対象遺伝子だけ再スケール（クリーンに）
obj[["RNA"]]@layers$scale.data <- NULL
obj <- ScaleData(obj, features = features_ordered, verbose = FALSE)

# ファイルへ直書き（重い描画を安定化）
png("results/heatmap_top100_percluster_ordered_percell.png",
    width = 4200, height = 3200, res = 300)
p_hm <- DoHeatmap(
  obj,
  features = features_ordered,
  group.by = "seurat_clusters",   # DimPlot と完全対応
  assay = "RNA", slot = "scale.data",
  cells = cells_ordered,
  raster = TRUE, draw.lines = FALSE,
  disp.min = -1.8, disp.max = 1.8
) +
  labs(title = "Top 100 per cluster — cells ordered by module score (0→21)") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(face = "italic", size = 6))
print(p_hm)
dev.off()

