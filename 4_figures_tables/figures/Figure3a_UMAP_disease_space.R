##########################################
#   Plot 45 PFAS into 148 disease space  #
#   バブルサイズ = nF（フッ素原子数）    #
#   ※ lemon非依存版                     #
##########################################
set.seed(224)   # 元スクリプトと同じseedで再現性を確保

library(ggplot2)   # lemonより先にロード（+.gg上書き防止）
library(dplyr)
library(umap)
library(openxlsx)
# library(lemon)  ← 除外

# ============================================================
# 1. データ読み込み
# ============================================================

## disease space UMAP結果
load("/Users/iida/Documents/Lab_Iida/2023_Iida/Interactome/Human_PPI_PFAS/Analysis_250707/list_umap_45PFAS_148disease.dat")

## PFAS略称
in_f  <- "/Users/iida/Documents/論文/Published/2023_Maeda/参考文献/PFOA_abbrivation_list_pfasmaster-2023-01-31.xlsx"
abbri   <- read.xlsx(in_f)
Acronym <- setNames(abbri$Acronym, abbri$CASRN)

## フッ素原子数（nF）
in_d  <- "/Users/iida/Documents/Lab_Ishibashi/original_data/"
in_f3 <- "hsPPARalpha_PFASs_wash_descriptors_Sscore_IC50_data.xlsx"
discriptors <- read.xlsx(paste0(in_d, in_f3))
nF <- setNames(discriptors$a_nF, discriptors$CASRN)

# ============================================================
# 2. 階層的クラスタリング（k = 2）
# ============================================================

hc       <- hclust(dist(list_umap$distnon$layout), method = "complete")
clusters <- cutree(hc, k = 2)

cat("クラスター割り当て数:\n")
print(table(clusters))

# ============================================================
# 3. プロット用データフレームの作成
# ============================================================

umap_df          <- as.data.frame(list_umap$distnon$layout)
umap_df$cluster  <- as.factor(clusters)
umap_df$nF       <- nF[rownames(umap_df)]
# nFが取得できなかった場合は元スクリプトと同じ値で補完
umap_df[is.na(umap_df$nF), "nF"] <- c(5, 6, 9, 5, 4, 6, 5, 6, 17, 4, 6, 6)
umap_df$acronym  <- Acronym[rownames(umap_df)]
umap_df$accolors <- ifelse(umap_df$cluster == 1, "#577a00", "#9159ba")

# ============================================================
# 4. 可視化
# ============================================================

p <-
  ggplot(umap_df, aes(x = V1, y = V2, color = cluster, size = nF)) +
  geom_point() +
  scale_size(range = c(2, 10), name = "nF") +
  scale_color_manual(values = c("#7CAE00", "#C77CFF"),
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
# 5. 保存
# ============================================================

out_dir <- "/Users/iida/Documents/論文/2024_Iida/第二稿/Code_for_Manuscript"
ggsave(file.path(out_dir, "Figure3a_UMAP_disease_space.pdf"),
       plot = p, width = 8, height = 6, units = "in", device = "pdf")
cat("Saved: Figure3a_UMAP_disease_space.pdf\n")
