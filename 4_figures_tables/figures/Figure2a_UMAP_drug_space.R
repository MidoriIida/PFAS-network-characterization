##########################################
#   Plot according to 900 drug space     #
#              45 PFAS                   #
# 　  　    クラスタリング               #
##########################################
set.seed(123)  # 結果の再現性のため
library(ggplot2)
library(dplyr)
library(umap)
library(dplyr)
library(lemon)
library(reshape2)
library(cowplot)
library(openxlsx)

### 1. load ###################################################################
load("/path/to/your/local/directory/list_umap_PFAS_drug.dat")
load("/path/to/your/local/directory/umap_data_non.dat")

in_f <- "/path/to/your/local/directory/PFOA_abbrivation_list_pfasmaster-2023-01-31.xlsx"
abbri <- read.xlsx(in_f)
Acronym <- abbri$Acronym
names(Acronym) <- abbri$CASRN

### 3. 階層的クラスタリング #################################################################
#次元削減したデータに対してk-meansクラスタリングを適用
hc <- hclust(dist(list_umap$distnon$layout), method = "complete")

# デンドログラムのプロット
plot(hc)

# デンドログラムを3つのクラスタにカット
clusters <- cutree(hc, k = 2)

#save
list_cluster <- list(hc = hc, clusters = clusters)

### 4. クラスタリングの結果を色分けしてプロット ###############################################
umap_df         <- as.data.frame(list_umap$distnon$layout)

merged_df <- left_join(
  umap_df    %>% tibble::rownames_to_column("id"),
  umap_data_non %>% tibble::rownames_to_column("id"),
  by = "id"
)

umap_df         <- as.data.frame(list_umap$distnon$layout)
umap_df$cluster <- as.factor(clusters)
umap_df$nF      <- merged_df$nF
umap_df$acronym <- Acronym[rownames(umap_df)]
umap_df$accolors <- ifelse(umap_df$cluster == 1, "#9e4a44", "#006669")

p <- 
  ggplot(umap_df, aes(x = V1, y = V2, color = cluster, size = nF)) +
  geom_point() +　scale_size(range = c(3, 12)) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) +
  labs(x = "UMAP 1", y = "UMAP 2") + #title = "UMAP projection with hierarchical clusters", 
  theme_bw() + geom_text(aes(label = acronym), color = umap_df$accolors, size = 4, vjust = 0.5, hjust = 0.5) +  # ラベルをプロット
  theme(    axis.text.x  = element_text(size = 10),  # X軸の文字サイズ
            axis.text.y  = element_text(size = 10),  # Y軸の文字サイズ
            axis.title.x = element_text(size = 12),  # X軸のラベルサイズと太字設定
            axis.title.y = element_text(size = 12),   # Y軸のラベルサイズと太字設定
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10) ) +
  guides(color = guide_legend(override.aes = list(size = 3))) # 凡例の点のサイズを固定で5に設定
p

#ggsave("./path/to/intermediate/data/plot_UMAP.eps", plot = p, width = 24, height = 23, units = "in", device = "eps")
ggsave("Figure2a_UMAP_drug_space.pdf", plot = p, width = 8, height = 6, units = "in", device = "pdf")




