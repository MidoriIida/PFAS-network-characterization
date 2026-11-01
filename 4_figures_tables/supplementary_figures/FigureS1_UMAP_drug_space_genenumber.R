##########################################
#   Plot according to 900 drug space     #
#              45 PFAS                   #
#   バブルサイズ = PFAS関連遺伝子数      #
##########################################
set.seed(123)

library(ggplot2)   # lemonより先にロード（+.gg上書き防止）
library(dplyr)
library(umap)
library(reshape2)
library(cowplot)
library(openxlsx)
library(tibble)
# library(lemon)  ← 除外

# ============================================================
# 1. データ読み込み
# ============================================================

## UMAP・クラスタリング結果
load("/Users/iida/Documents/Lab_Iida/2023_Iida/Interactome/Human_PPI_PFAS/Analysis_240819/list_umap_PFAS_drug.dat")
load("/Users/iida/Documents/Lab_Iida/2023_Iida/Interactome/Human_PPI_PFAS/Analysis_240819/umap_data_non.dat")

## PFAS略称
in_f <- "/Users/iida/Documents/論文/Published/2023_Maeda/参考文献/PFOA_abbrivation_list_pfasmaster-2023-01-31.xlsx"
abbri   <- read.xlsx(in_f)
Acronym <- setNames(abbri$Acronym, abbri$CASRN)

## PFAS関連遺伝子データ
load("/Users/iida/Documents/Lab_Iida/2023_Iida/Interactome/Human_PPI_PFAS/Analysis_230701/0.make_data/human_PFAS_genes/human_PFAS_genes_230712_iida.dat")

# ============================================================
# 2. 階層的クラスタリング
# ============================================================

hc       <- hclust(dist(list_umap$distnon$layout), method = "complete")
clusters <- cutree(hc, k = 2)

list_cluster <- list(hc = hc, clusters = clusters)

# ============================================================
# 3. PFAS関連遺伝子数の計算
# ============================================================

cas_35 <- rownames(umap_data_non)   # 35 PFAS の CAS番号

n_genes <- sapply(cas_35, function(cas) {
  length(unique(list_PFAS_genes$Merged[list_PFAS_genes$Merged$CASRN == cas, "GeneID"]))
})

cat("遺伝子数サマリー:\n")
print(summary(n_genes))

# ============================================================
# 4. プロット用データフレームの作成
# ============================================================

umap_df          <- as.data.frame(list_umap$distnon$layout)
umap_df$cluster  <- as.factor(clusters)
umap_df$n_genes  <- n_genes[rownames(umap_df)]
umap_df$log_genes <- log10(umap_df$n_genes + 1)   # log10変換（1-1956の幅を圧縮）
umap_df$acronym  <- Acronym[rownames(umap_df)]
umap_df$accolors <- ifelse(umap_df$cluster == 1, "#9e4a44", "#006669")

cat("\nlog10遺伝子数サマリー:\n")
print(summary(umap_df$log_genes))

# ============================================================
# 5. 可視化
# ============================================================

p <-
  ggplot(umap_df, aes(x = V1, y = V2, color = cluster, size = log_genes)) +
  geom_point() +
  scale_size(range = c(2, 10),
             name  = "Associated\ngenes (log10)") +
  scale_color_manual(values = c("#F8766D", "#00BFC4"),
                     name   = "cluster") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_bw() +
  geom_text(aes(label = acronym),
            color = umap_df$accolors,
            size  = 3,
            vjust = 0.5, hjust = 0.5) +
  theme(
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text  = element_text(size = 10),
    legend.title = element_text(size = 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))
p

# ============================================================
# 6. 保存
# ============================================================

out_dir  <- "/Users/iida/Documents/論文/2024_Iida/第二稿/Code_for_Manuscript/Supplymentary_Figures/"
ggsave(file.path(out_dir, "FigureS1_UMAP_drug_space_genenumber.pdf"),
       plot = p, width = 8, height = 6, units = "in", device = "pdf")
cat("Saved: FigureS1_UMAP_drug_space_genenumber.pdf\n")
