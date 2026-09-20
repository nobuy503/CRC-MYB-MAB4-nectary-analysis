# Exploratory reclustering and heatmaps within integrated cluster 15
#
# Input:
#   data/processed/integrated_data.rds
#
# Analyses included:
# - Subsetting integrated Seurat cluster 15.
# - Reclustering and UMAP of cluster 15 cells.
# - Annotation by Flower or Silique origin.
# - Feature plots for nectary-associated genes.
# - Filtering genes by the fraction of cluster 15 cells expressing them.
# - Per-cell heatmaps of highly expressed genes with subcluster and tissue
#   annotations.
#
# Figure relationship:
# This is an exploratory analysis of heterogeneity within cluster 15. It is not
# the Flower-versus-Silique top-100 DEG heatmap shown in Fig. 1d, and it does
# not directly reproduce the final Flower-only pseudotime panels in Fig. 1p-t.
#
# Provenance note:
# The file contains repeated exploratory blocks and alternative thresholds.
# It is retained as analysis provenance rather than a streamlined final script.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

BASEDIR <- "results/cluster15_reclustering_exploration"
FIGDIR  <- file.path(BASEDIR, "figs")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

INTEGRATED_RDS <- "data/processed/integrated_data.rds"

## 0) オブジェクト読み込み（すでに obj があればそれを使ってもOK）
if (!exists("obj")) {
  stopifnot(file.exists(INTEGRATED_RDS))
  obj <- readRDS(INTEGRATED_RDS)
}

## 1) tissue 情報（._tissue）を確認・補完
if (!"._tissue" %in% colnames(obj@meta.data)) {
  tissue_guess <- NA_character_
  if ("orig.ident" %in% colnames(obj@meta.data)) {
    oi <- tolower(obj$orig.ident)
    tissue_guess <- ifelse(grepl("flower|flw|flwr", oi), "flower",
                           ifelse(grepl("silique|sil", oi), "silique", NA))
  }
  obj$._tissue <- tissue_guess
}

cat("Cells per tissue:\n")
print(table(obj$._tissue, useNA = "ifany"))

## 2) クラスタ15だけサブセット

stopifnot("seurat_clusters" %in% colnames(obj@meta.data))
Idents(obj) <- obj$seurat_clusters
cat("Original clusters:\n"); print(table(Idents(obj)))

stopifnot("15" %in% levels(Idents(obj)))
obj15 <- subset(obj, idents = "15")
cat("\nCluster 15 cells per tissue:\n")
print(table(obj15$._tissue, useNA = "ifany"))

## 3) PCA/UMAP に使うアッセイを決める（integrated → SCT → RNA の優先順位）

assay_graph <- if ("integrated" %in% Assays(obj15)) {
  "integrated"
} else if ("SCT" %in% Assays(obj15)) {
  "SCT"
} else {
  "RNA"
}
cat("\nUsing assay for PCA/UMAP:", assay_graph, "\n")
DefaultAssay(obj15) <- assay_graph

## 4) VariableFeatures の保証と前処理

if (assay_graph %in% c("integrated", "SCT")) {
  if (length(VariableFeatures(obj15)) == 0) {
    VariableFeatures(obj15) <- SelectIntegrationFeatures(list(obj15), nfeatures = 3000)
  }
} else {  # RNA の場合
  need_norm <- tryCatch({
    m <- GetAssayData(obj15, assay = "RNA", layer = "data")
    ncol(m) == 0
  }, error = function(e) TRUE)
  if (need_norm) {
    obj15 <- NormalizeData(obj15, verbose = FALSE)
  }
  obj15 <- FindVariableFeatures(obj15, selection.method = "vst",
                                nfeatures = 3000, verbose = FALSE)
  obj15 <- ScaleData(obj15, features = VariableFeatures(obj15), verbose = FALSE)
}

## 5) PCA → Neighbors → Clusters → UMAP

set.seed(1)
N_PCS <- 30

obj15 <- RunPCA(
  obj15,
  assay    = assay_graph,
  features = VariableFeatures(obj15),
  npcs     = N_PCS,
  verbose  = FALSE
)
obj15 <- FindNeighbors(obj15, dims = 1:N_PCS, verbose = FALSE)
obj15 <- FindClusters(obj15, resolution = 0.4, verbose = FALSE)
obj15 <- RunUMAP(obj15, dims = 1:N_PCS, verbose = FALSE)

cat("\nSubclusters inside cluster 15 (new seurat_clusters):\n")
print(table(obj15$seurat_clusters))

## 6) UMAP プロット（クラスタ別 / tissue 別）

p1 <- DimPlot(
  obj15,
  reduction = "umap",
  group.by  = "seurat_clusters",
  label = TRUE, repel = TRUE
) + ggtitle("Cluster 15 — reclustered (integrated)")

p2 <- DimPlot(
  obj15,
  reduction = "umap",
  group.by  = "._tissue"
) + ggtitle("Cluster 15 — tissue (flower vs silique)")

print(p1)
print(p2)


p1 <- DimPlot(
  obj15,
  reduction = "umap",
  group.by  = "seurat_clusters",
  label = FALSE, repel = TRUE
) + ggtitle("Cluster 15 — reclustered (integrated)")

print(p1)

ggsave(file.path(FIGDIR, "cluster15_integrated_UMAP_clusters.pdf"),
       p1, width = 6, height = 5)
ggsave(file.path(FIGDIR, "cluster15_integrated_UMAP_tissue.pdf"),
       p2, width = 6, height = 5)

## 7) 次の解析用に保存
OUT_RDS <- file.path(BASEDIR, "cluster15_integrated_reclustered.rds")
saveRDS(obj15, OUT_RDS)
cat("\nSaved object to:", normalizePath(OUT_RDS), "\n")







#1. 準備：cluster15オブジェクトを読み込み & RNA を可視化用にセット

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

BASEDIR   <- "results/cluster15_reclustering_exploration"
OBJ15_RDS <- file.path(BASEDIR, "cluster15_integrated_reclustered.rds")

stopifnot(file.exists(OBJ15_RDS))
obj15 <- readRDS(OBJ15_RDS)

## 可視化には RNA を使う
if ("RNA" %in% Assays(obj15)) {
  DefaultAssay(obj15) <- "RNA"
  
  # data layer が空なら正規化だけしておく
  need_norm <- tryCatch({
    m <- GetAssayData(obj15, assay = "RNA", layer = "data")
    ncol(m) == 0
  }, error = function(e) {
    m2 <- tryCatch(GetAssayData(obj15, assay = "RNA", slot = "data"),
                   error = function(e2) NULL)
    is.null(m2) || ncol(m2) == 0
  })
  
  if (need_norm) {
    message("Normalizing RNA assay for visualization ...")
    obj15 <- NormalizeData(obj15, verbose = FALSE)
  }
} else {
  stop("RNA assay が obj15 にありません。")
}



plot_umap_gene <- function(obj, gene) {
  rn  <- rownames(obj)
  hit <- grep(paste0("^", gene, "(\\.\\d+)?$"), rn, value = TRUE)
  
  if (length(hit) == 0) {
    stop("Gene not found in RNA assay: ", gene)
  }
  gene_use <- hit[1]
  cat("Using feature:", gene_use, "\n")
  
  p <- FeaturePlot(
    obj,
    features  = gene_use,
    reduction = "umap",
    order     = TRUE,
    pt.size   = 0.6,
    # 微量発現も見たいのでカットオフはかけない
    min.cutoff = NA
  ) + ggtitle(paste0(gene_use, " expression (UMAP)"))
  
  print(p)
  
  # 画像保存したければ（任意）
  out_png <- file.path(BASEDIR, paste0(gene_use, "_UMAP.png"))
  ggsave(out_png, p, width = 5, height = 4, dpi = 300)
  cat("Saved:", out_png, "\n")
}








# CRC
plot_umap_gene(obj15, "AT1G69180")

# 他の候補
plot_umap_gene(obj15, "AT4G31820")
plot_umap_gene(obj15, "AT2G39060")
plot_umap_gene(obj15, "AT3G01530")
plot_umap_gene(obj15, "AT3G27810")
plot_umap_gene(obj15, "AT1G65970")




suppressPackageStartupMessages({
  library(Seurat)
})

BASEDIR   <- "results/cluster15_reclustering_exploration"
OBJ15_RDS <- file.path(BASEDIR, "cluster15_integrated_reclustered.rds")

obj15 <- readRDS(OBJ15_RDS)

# RNAアッセイを使う
stopifnot("RNA" %in% Assays(obj15))
rna <- obj15[["RNA"]]

genes_all <- rownames(rna)
length(genes_all)          # 何遺伝子あるか確認

# コンソールに少しだけ表示
head(genes_all, 20)

# CSVに保存（1列だけ）
out_all <- file.path(BASEDIR, "cluster15_genes_all_RNA.csv")
write.csv(data.frame(gene = genes_all),
          file = out_all, row.names = FALSE)
cat("Saved all RNA genes to:", out_all, "\n")











suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

BASEDIR   <- "results/cluster15_reclustering_exploration"
OBJ15_RDS <- file.path(BASEDIR, "cluster15_integrated_reclustered.rds")

obj15 <- readRDS(OBJ15_RDS)
DefaultAssay(obj15) <- "RNA"

rna_assay <- obj15[["RNA"]]

layer_names <- Layers(rna_assay)
cat("Layers in RNA assay:\n")
print(layer_names)

## ---- Silique / Flower をそれぞれ取得 ----
mat_s <- NULL
mat_f <- NULL

if ("data.Silique" %in% layer_names) {
  mat_s <- LayerData(rna_assay, "data.Silique")
  cat("data.Silique dim:", paste(dim(mat_s), collapse=" x "), "\n")
}
if ("data.Flower" %in% layer_names) {
  mat_f <- LayerData(rna_assay, "data.Flower")
  cat("data.Flower dim:", paste(dim(mat_f), collapse=" x "), "\n")
}

if (is.null(mat_s) && is.null(mat_f)) {
  stop("data.Silique / data.Flower が見つかりません。")
}

## ---- 各レイヤーごとに「1細胞でも >0 の遺伝子」を取得 ----
expr_s <- character(0)
expr_f <- character(0)

if (!is.null(mat_s)) {
  expr_s <- rownames(mat_s)[Matrix::rowSums(mat_s > 0) > 0]
  cat("Silique: expressed genes:", length(expr_s), "\n")
}
if (!is.null(mat_f)) {
  expr_f <- rownames(mat_f)[Matrix::rowSums(mat_f > 0) > 0]
  cat("Flower : expressed genes:", length(expr_f), "\n")
}

## ---- union をとって「クラスタ15で発現している遺伝子」の集合に ----
expr_genes <- sort(unique(c(expr_s, expr_f)))
cat("Total expressed genes in cluster 15 (Silique ∪ Flower):", length(expr_genes), "\n")

head(expr_genes, 20)

## ---- CSV に保存 ----
out_expr <- file.path(BASEDIR, "cluster15_genes_expressed_RNA.csv")
write.csv(data.frame(gene = expr_genes),
          file = out_expr, row.names = FALSE)
cat("Saved expressed RNA genes to:", out_expr, "\n")





suppressPackageStartupMessages({
  library(Matrix)
})

# さっきと同じく data.Silique / data.Flower を使う
mat_s <- LayerData(rna_assay, "data.Silique")
mat_f <- LayerData(rna_assay, "data.Flower")

# 細胞ごとの発現 >0 をカウント
ncell_s <- Matrix::rowSums(mat_s > 0)
ncell_f <- Matrix::rowSums(mat_f > 0)
ncell_tot <- ncell_s + ncell_f

# クラスタ15全体の細胞数
n_cells_total <- ncol(mat_s) + ncol(mat_f)
n_cells_total

# 例：5%以上の細胞で発現している遺伝子
min_prop <- 0.05  # 5%
min_cells <- ceiling(min_prop * n_cells_total)

genes_5pct <- names(ncell_tot)[ncell_tot >= min_cells]
length(genes_5pct)
head(genes_5pct, 20)

# CSVに保存
out_5pct <- file.path(BASEDIR, "cluster15_genes_expressed_5pct.csv")
write.csv(data.frame(gene = genes_5pct),
          file = out_5pct, row.names = FALSE)
cat("Saved genes expressed in ≥5% of cluster15 cells to:", out_5pct, "\n")






genes_10cells <- names(ncell_tot)[ncell_tot >= 10]
length(genes_10cells)
head(genes_10cells, 20)

out_10 <- file.path(BASEDIR, "cluster15_genes_expressed_10cells.csv")
write.csv(data.frame(gene = genes_10cells),
          file = out_10, row.names = FALSE)
cat("Saved genes expressed in ≥10 cells in cluster15 to:", out_10, "\n")







suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

BASEDIR   <- "results/cluster15_reclustering_exploration"
OBJ15_RDS <- file.path(BASEDIR, "cluster15_integrated_reclustered.rds")

obj15 <- readRDS(OBJ15_RDS)
DefaultAssay(obj15) <- "RNA"

rna_assay <- obj15[["RNA"]]

## ---- レイヤー確認 ----
layer_names <- Layers(rna_assay)
cat("Layers in RNA assay:\n")
print(layer_names)

## ---- data.Silique / data.Flower を個別に取得 ----
mat_s <- NULL
mat_f <- NULL

if ("data.Silique" %in% layer_names) {
  mat_s <- LayerData(rna_assay, "data.Silique")
  cat("data.Silique dim:", paste(dim(mat_s), collapse=" x "), "\n")
}
if ("data.Flower" %in% layer_names) {
  mat_f <- LayerData(rna_assay, "data.Flower")
  cat("data.Flower dim:", paste(dim(mat_f), collapse=" x "), "\n")
}

if (is.null(mat_s) && is.null(mat_f)) {
  stop("data.Silique / data.Flower が見つかりません。")
}

## ---- 各レイヤーで「発現細胞数」を数える ----
ncell_s <- 0
ncell_f <- 0

if (!is.null(mat_s)) {
  ncell_s <- Matrix::rowSums(mat_s > 0)
}
if (!is.null(mat_f)) {
  ncell_f <- Matrix::rowSums(mat_f > 0)
}

# ベクトル長を揃える（行名ベースで合わせる）
all_genes <- union(rownames(mat_s), rownames(mat_f))
ncell_s_full <- setNames(integer(length(all_genes)), all_genes)
ncell_f_full <- setNames(integer(length(all_genes)), all_genes)

if (!is.null(mat_s)) ncell_s_full[rownames(mat_s)] <- ncell_s
if (!is.null(mat_f)) ncell_f_full[rownames(mat_f)] <- ncell_f

ncell_tot <- ncell_s_full + ncell_f_full

## ---- クラスタ15全体の細胞数 ----
n_cells_total <- (if (!is.null(mat_s)) ncol(mat_s) else 0) +
  (if (!is.null(mat_f)) ncol(mat_f) else 0)
cat("Total cells in cluster 15:", n_cells_total, "\n")

## ---- 10% 以上で発現している遺伝子だけ抽出 ----
min_prop  <- 0.10           # ここが 10%
min_cells <- ceiling(min_prop * n_cells_total)
cat("Threshold (cells):", min_cells, "cells (", min_prop*100, "% )\n")

genes_10pct <- names(ncell_tot)[ncell_tot >= min_cells]
cat("Genes expressed in ≥10% of cluster15 cells:", length(genes_10pct), "\n")

head(genes_10pct, 20)

## ---- CSV に保存 ----
out_10pct <- file.path(BASEDIR, "cluster15_genes_expressed_10pct.csv")
write.csv(data.frame(gene = genes_10pct),
          file = out_10pct, row.names = FALSE)
cat("Saved genes (≥10% of cells) to:", out_10pct, "\n")






#ヒートマップ

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(pheatmap)   # 無ければ install.packages("pheatmap")
})

BASEDIR   <- "results/cluster15_reclustering_exploration"
OBJ15_RDS <- file.path(BASEDIR, "cluster15_integrated_reclustered.rds")

obj15 <- readRDS(OBJ15_RDS)
DefaultAssay(obj15) <- "RNA"

rna_assay   <- obj15[["RNA"]]
layer_names <- Layers(rna_assay)
cat("Layers in RNA assay:\n"); print(layer_names)

## 1) Silique / Flower の data レイヤーを取り出す
mat_s <- if ("data.Silique" %in% layer_names) LayerData(rna_assay, "data.Silique") else NULL
mat_f <- if ("data.Flower"  %in% layer_names) LayerData(rna_assay, "data.Flower")  else NULL

if (is.null(mat_s) && is.null(mat_f)) {
  stop("data.Silique / data.Flower が見つかりません。")
}

cat("data.Silique dim:", if (!is.null(mat_s)) paste(dim(mat_s), collapse=" x ") else "NULL", "\n")
cat("data.Flower  dim:", if (!is.null(mat_f)) paste(dim(mat_f), collapse=" x ") else "NULL", "\n")

## 2) genes_use に行を絞る（見つかったものだけ）
rn_s <- if (!is.null(mat_s)) rownames(mat_s) else character(0)
rn_f <- if (!is.null(mat_f)) rownames(mat_f) else character(0)

g_ok_s <- if (!is.null(mat_s)) genes_use[genes_use %in% rn_s] else character(0)
g_ok_f <- if (!is.null(mat_f)) genes_use[genes_use %in% rn_f] else character(0)
g_ok   <- sort(unique(c(g_ok_s, g_ok_f)))

cat("Genes present in Silique/Flower layers:", length(g_ok), "\n")
if (length(g_ok) == 0) stop("genes_use が data.* のどちらにも見つかりません。")

mat_s_sub <- if (!is.null(mat_s)) mat_s[g_ok, , drop = FALSE] else NULL
mat_f_sub <- if (!is.null(mat_f)) mat_f[g_ok, , drop = FALSE] else NULL

## 3) 列方向に結合して「遺伝子×細胞」行列を作る
mats_list <- list(mat_f_sub, mat_s_sub)
mats_list <- mats_list[!sapply(mats_list, is.null)]   # NULL は除外
expr_mat  <- do.call(cbind, mats_list)

# 念のため、obj15 のセルと一致しているか確認
cat("ncol(expr_mat):", ncol(expr_mat), " / ncol(obj15):", ncol(obj15), "\n")

## 4) seurat_clusters で列を並べ替え（サブクラスタ順）
clusters <- obj15$seurat_clusters
tissues  <- obj15$._tissue

cells_order <- colnames(obj15)[order(clusters, tissues)]
# expr_mat の列を、この順番に並べ替え
cells_order <- intersect(cells_order, colnames(expr_mat))
expr_mat    <- expr_mat[, cells_order, drop = FALSE]

## 5) 列アノテーション（サブクラスタ & tissue）を付ける
annot_col <- data.frame(
  cluster = clusters[cells_order],
  tissue  = tissues[cells_order]
)
rownames(annot_col) <- cells_order

## 6) pheatmap で描画
pheatmap(
  expr_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annot_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  main = "Cluster 15 — Top genes (>=10% cells)"
)

# 保存したければ
OUTPNG <- file.path(BASEDIR, "cluster15_heatmap_top10pct_topN_pheatmap.png")
png(OUTPNG, width = 1200, height = 1600, res = 150)
pheatmap(
  expr_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annot_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  main = "Cluster 15 — Top genes (>=10% cells)"
)
dev.off()
cat("Saved heatmap:", OUTPNG, "\n")



#ヒートマップ
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(pheatmap)
  library(ggplot2)
  library(scales)  # hue_pal() 用
})

library(viridis)

## ===== パスなど基本設定 =====
BASEDIR    <- "results/cluster15_reclustering_exploration"
OBJ15_RDS  <- file.path(BASEDIR, "cluster15_integrated_reclustered.rds")
GENE10CSV  <- file.path(BASEDIR, "cluster15_genes_expressed_10pct.csv")

stopifnot(file.exists(OBJ15_RDS), file.exists(GENE10CSV))

## ===== オブジェクト & 遺伝子リスト読み込み =====
obj15 <- readRDS(OBJ15_RDS)
DefaultAssay(obj15) <- "RNA"

genes_10pct <- read.csv(GENE10CSV, stringsAsFactors = FALSE)$gene
cat("10% genes:", length(genes_10pct), "\n")

## ===== RNA assay / layer 確認 =====
rna_assay   <- obj15[["RNA"]]
layer_names <- Layers(rna_assay)
cat("Layers in RNA assay:\n"); print(layer_names)

mat_s <- if ("data.Silique" %in% layer_names) LayerData(rna_assay, "data.Silique") else NULL
mat_f <- if ("data.Flower"  %in% layer_names) LayerData(rna_assay, "data.Flower")  else NULL

if (is.null(mat_s) && is.null(mat_f)) {
  stop("data.Silique / data.Flower が見つかりません。")
}

cat("data.Silique dim:", if (!is.null(mat_s)) paste(dim(mat_s), collapse=" x ") else "NULL", "\n")
cat("data.Flower  dim:", if (!is.null(mat_f)) paste(dim(mat_f), collapse=" x ") else "NULL", "\n")

## ===== 各遺伝子の「発現細胞数」を計算 =====
ncell_s <- if (!is.null(mat_s)) Matrix::rowSums(mat_s > 0) else integer(0)
ncell_f <- if (!is.null(mat_f)) Matrix::rowSums(mat_f > 0) else integer(0)

all_genes <- union(names(ncell_s), names(ncell_f))
ncell_s_full <- setNames(integer(length(all_genes)), all_genes)
ncell_f_full <- setNames(integer(length(all_genes)), all_genes)
if (!is.null(mat_s)) ncell_s_full[names(ncell_s)] <- ncell_s
if (!is.null(mat_f)) ncell_f_full[names(ncell_f)] <- ncell_f
ncell_tot <- ncell_s_full + ncell_f_full

## ===== 10%遺伝子の中から「発現細胞数トップ100」を選ぶ =====
ncell_sub <- ncell_tot[genes_10pct]
ncell_sub <- ncell_sub[!is.na(ncell_sub)]

topN <- 100  # ★ここを 100 にしました
genes_top <- names(sort(ncell_sub, decreasing = TRUE))[1:min(topN, length(ncell_sub))]
cat("Top genes used in heatmap:", length(genes_top), "\n")

## ===== rownames（バージョン付き）に対応した genes_use を作る =====
rn <- rownames(obj15[["RNA"]])
genes_use <- character(0)
for (g in genes_top) {
  hit <- grep(paste0("^", g, "(\\.\\d+)?$"), rn, value = TRUE)
  if (length(hit) > 0) genes_use <- c(genes_use, hit[1])
}
genes_use <- unique(genes_use)
cat("genes_use (with version):", length(genes_use), "\n")
print(genes_use)

## ===== Silique / Flower 行列を genes_use でサブセット & 結合 =====
rn_s <- if (!is.null(mat_s)) rownames(mat_s) else character(0)
rn_f <- if (!is.null(mat_f)) rownames(mat_f) else character(0)

g_ok_s <- if (!is.null(mat_s)) genes_use[genes_use %in% rn_s] else character(0)
g_ok_f <- if (!is.null(mat_f)) genes_use[genes_use %in% rn_f] else character(0)
g_ok   <- sort(unique(c(g_ok_s, g_ok_f)))
cat("Genes present in data.Silique / data.Flower:", length(g_ok), "\n")

mat_s_sub <- if (!is.null(mat_s)) mat_s[g_ok, , drop = FALSE] else NULL
mat_f_sub <- if (!is.null(mat_f)) mat_f[g_ok, , drop = FALSE] else NULL

mats_list <- list(mat_f_sub, mat_s_sub)
mats_list <- mats_list[!sapply(mats_list, is.null)]
expr_mat  <- do.call(cbind, mats_list)

cat("expr_mat dim:", paste(dim(expr_mat), collapse=" x "), "\n")
cat("ncol(expr_mat) vs ncol(obj15):", ncol(expr_mat), ncol(obj15), "\n")

## ===== 列順とアノテーション（cluster / tissue） =====
Idents(obj15) <- obj15$seurat_clusters
clusters <- obj15$seurat_clusters
tissues  <- obj15$._tissue

cells_order <- colnames(obj15)[order(clusters, tissues)]
cells_order <- intersect(cells_order, colnames(expr_mat))
expr_mat    <- expr_mat[, cells_order, drop = FALSE]

annot_col <- data.frame(
  cluster = clusters[cells_order],
  tissue  = tissues[cells_order]
)
rownames(annot_col) <- cells_order

## ===== アノテーションカラー =====
# クラスター：Seurat (=ggplot) デフォルト順
cluster_levels  <- levels(Idents(obj15))
cluster_palette <- hue_pal()(length(cluster_levels))
names(cluster_palette) <- cluster_levels

# tissue：指定色（ggplot デフォルト1/2）
tissue_cols <- c(
  flower  = "#F8766D",  # サーモンピンク
  silique = "#00BFC4"   # ターコイズブルー
)

ann_colors <- list(
  cluster = cluster_palette,
  tissue  = tissue_cols
)

## ===== 発現値のカラースケール（0〜4に固定 & 紫→黒→黄） =====
hm_colors <- viridis::viridis(100)

# 発現値を 0〜6 の範囲にクリップ
expr_plot <- as.matrix(expr_mat)
expr_plot[expr_plot < 0] <- 0
expr_plot[expr_plot > 6] <- 6

# breaks を 0〜6 に合わせる（色100段 × 区切り101）
hm_breaks <- seq(0, 6, length.out = length(hm_colors) + 1)

## ===== pheatmap で描画 & 保存 =====
pheatmap(
  expr_plot,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  annotation_col    = annot_col,
  annotation_colors = ann_colors,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  color             = hm_colors,
  breaks            = hm_breaks,  # ★0〜4に固定
  # scale           = "none",    # 明示しておいてもOK
  main              = "Cluster 15 — Top 100 genes (>=10% cells, expr 0–4)"
)

OUTPNG <- file.path(BASEDIR, "cluster15_heatmap_top10pct_top100_0to4.png")
png(OUTPNG, width = 1200, height = 1600, res = 150)
pheatmap(
  expr_plot,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  annotation_col    = annot_col,
  annotation_colors = ann_colors,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  color             = hm_colors,
  breaks            = hm_breaks,
  main              = "Cluster 15 — Top 100 genes (>=10% cells, expr 0–4)"
)
dev.off()
cat("Saved heatmap with custom colors:", OUTPNG, "\n")



