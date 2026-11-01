##########################################
#   Plot 45 PFAS into 148 disease space  #
#   バブルサイズ = PFAS関連遺伝子数      #
#   ※ lemon非依存版                     #
##########################################
set.seed(224)   # 元スクリプトと同じseedで再現性を確保

library(ggplot2)   # lemonより先にロード（+.gg上書き防止）
library(dplyr)
library(umap)
library(openxlsx)
library(tibble)
# library(lemon)  ← 除外

# ============================================================
# 1. データ読み込み
# ============================================================

## disease space UMAP結果
load("/path/to/your/local/directory/list_umap_45PFAS_148disease.dat")

## PFAS略称
in_f  <- "/path/to/your/local/directory/PFOA_abbrivation_list_pfasmaster-2023-01-31.xlsx"
abbri   <- read.xlsx(in_f)
Acronym <- setNames(abbri$Acronym, abbri$CASRN)

## PFAS関連遺伝子データ
load("/path/to/your/local/directory/human_PFAS_genes_230712_iida.dat")

# ============================================================
# 2. 階層的クラスタリング（k = 2）
# ============================================================

hc       <- hclust(dist(list_umap$distnon$layout), method = "complete")
clusters <- cutree(hc, k = 2)

cat("クラスター割り当て数:\n")
print(table(clusters))

# ============================================================
# 3. PFAS関連遺伝子数の計算
# ============================================================

cas_disease <- rownames(list_umap$distnon$layout)   # 44 PFAS の CAS番号

n_genes <- sapply(cas_disease, function(cas) {
  length(unique(list_PFAS_genes$Merged[list_PFAS_genes$Merged$CASRN == cas, "GeneID"]))
})

cat("\n遺伝子数サマリー:\n")
print(summary(n_genes))

# ============================================================
# 4. プロット用データフレームの作成
# ============================================================

umap_df           <- as.data.frame(list_umap$distnon$layout)
umap_df$cluster   <- as.factor(clusters)
umap_df$n_genes   <- n_genes[rownames(umap_df)]
umap_df$log_genes <- log10(umap_df$n_genes + 1)   # log10変換
umap_df$acronym   <- Acronym[rownames(umap_df)]
umap_df$accolors  <- ifelse(umap_df$cluster == 1, "#577a00", "#9159ba")  # 元スクリプトのテキスト色

# ============================================================
# 5. 可視化
# ============================================================

p <-
  ggplot(umap_df, aes(x = V1, y = V2, color = cluster, size = log_genes)) +
  geom_point() +
  scale_size(range = c(2, 10),
             name  = "Associated\ngenes (log₁₀)") +
  scale_color_manual(values = c("#7CAE00", "#C77CFF"),   # 元スクリプトの色を維持
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

out_dir <- "/path/to/your/local/directory/Code_for_Manuscript"
ggsave(file.path(out_dir, "FigureS5_UMAP_disease_space.pdf"),
       plot = p, width = 8, height = 6, units = "in", device = "pdf")
cat("Saved: FigureS5_UMAP_disease_space.pdf\n")
