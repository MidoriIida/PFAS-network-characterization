##########################################
#   Plot according to 900 drug space     #
#              45 PFAS                   #
#   　クラスターの意味 (ATC)             #
##########################################
set.seed(123)  # 結果の再現性のため
library(ggplot2)
library(dplyr)
library(umap)
library(dplyr)
library(lemon)
library(reshape2)
library(cowplot)
library(stringr)

### 1. load ###################################################################
load("./Analysis_240819/list_cluster.dat")

in_d2  <- "/path/to/your/local/directory/orignail_data"
in_f2 <- "/20230830/result_enrichment_ATC_code_KEGGdrugBackground_20230825.csv"

ATC        <- read.csv(paste0(in_d2, in_f2), row.names = 1)

#FDR補正
# データを行列に変換
p_matrix <- as.matrix(ATC)

# p値を1次元ベクトル化 → FDR補正 → 行列に戻す
p_fdr <- matrix(
  p.adjust(as.vector(p_matrix), method = "fdr"),
  nrow = nrow(p_matrix),
  ncol = ncol(p_matrix)
)

# 行名・列名を戻す
rownames(p_fdr) <- rownames(p_matrix)
colnames(p_fdr) <- colnames(p_matrix)

### 2. cluster ################################################################
clusters <- list_cluster$clusters

c1 <- names(clusters[clusters == 1])
c2 <- names(clusters[clusters == 2])

cas <- sapply(strsplit(gsub("^X", "", colnames(p_fdr)), "\\.\\."), function(x){gsub("\\.", "-",x[1])})
names(cas) <- cas

### 3. 階層的クラスタリングとATCの関係 #######################################################
c1_ATC <- p_fdr[cas %in% c1,]
c2_ATC <- p_fdr[cas %in% c2,]

c1_ATC_sig <- apply(c1_ATC, 1, function(x){sum(na.omit(x) < 0.05)})
c2_ATC_sig <- apply(c2_ATC, 1, function(x){sum(na.omit(x) < 0.05)})

ATC_sig    <- cbind(c1_ATC_sig, c2_ATC_sig)
y <-  as.data.frame(ATC_sig[which(rowSums(ATC_sig) != 0),]) %>%
  arrange(c2_ATC_sig, desc(c1_ATC_sig))
df_ATC_sig <- melt(as.matrix(y))
df_ATC_sig$Label <- str_wrap(df_ATC_sig[,1], width = 30)
df_ATC_sig$group <- gsub("c1_ATC_sig", "cluster 1", df_ATC_sig$Var2)
df_ATC_sig$group <- gsub("c2_ATC_sig", "cluster 2", df_ATC_sig$group)
levels(df_ATC_sig$group) <- c("cluster 1", "cluster 2")

colnames(df_ATC_sig)   <-  c("ATC", "cluster", "count", "Label", "group")

### 4. バープロットを作成 ######################################
g_ATC <- 
  ggplot(df_ATC_sig, aes(x = ATC, y = count, fill = group)) +  # x軸にギア、y軸に数を設定
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフの高さをデータの値に基づいて設定
  labs(#title = "Significant ATC class Counts by cluster", 
    x = "ATC class", y = "Count",  fill = "Cluster") +
  theme_bw() 

### 5. バープロットを調整（横） ######################################
g_ATC_h <- g_ATC + scale_x_discrete(labels = df_ATC_sig$Label) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = c(0.95, 0.95),  # 凡例をプロット内の右上に配置
        legend.justification = c("right", "top"))  # 凡例のアンカーポイントを右上に設定
g_ATC_h

### 5. バープロットを調整（縦） ######################################
g_ATC_v <- g_ATC +  coord_flip()  + scale_x_discrete(labels = df_ATC_sig$Label)
g_ATC_v

### 6. save ##########################################################
list_graph_enriched_ATC <- list(h = g_ATC_h, v = g_ATC_v, n = g_ATC)
save(list_graph_enriched_ATC, file = "./Analysis_240819/6.Enrichment_analysis/list_graph_enriched_ATC_fdr.dat")


