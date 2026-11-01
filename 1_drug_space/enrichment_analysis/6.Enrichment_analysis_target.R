##########################################
#   Plot according to 900 drug space     #
#              45 PFAS                   #
#   　クラスターの意味 (target)          #
##########################################
set.seed(123)  # 結果の再現性のため
library(ggplot2)
library(dplyr)
library(umap)
library(dplyr)
library(lemon)
library(reshape2)
library(cowplot)

### 1. load ###################################################################
load("./Analysis_240819/list_cluster.dat")

in_d2  <- "/path/to/your/local/directory/orignail_data"
in_f1 <- "/result_enrichment_Target_class_20230811.csv"

Target      <- read.csv(paste0(in_d2, in_f1), row.names = 1)

### 2. cluster ################################################################
clusters <- list_cluster$clusters

c1 <- names(clusters[clusters == 1])
c2 <- names(clusters[clusters == 2])

cas <- sapply(strsplit(gsub("^X", "", colnames(Target)), "\\.\\."), function(x){gsub("\\.", "-",x[1])})
names(cas) <- cas

### 3. 階層的クラスタリングとTargetの関係 ######################################
c1_Target <- Target[cas %in% c1,]
c2_Target <- Target[cas %in% c2,]

c1_Target_sig <- apply(c1_Target, 1, function(x){sum(na.omit(x) < 0.05)})
c2_Target_sig <- apply(c2_Target, 1, function(x){sum(na.omit(x) < 0.05)})

# データフレームを作成
Target_sig <- cbind(c1_Target_sig, c2_Target_sig)
x <- Target_sig[which(rowSums(Target_sig) != 0),]
df_Target_sig <- melt(x[order(x[,2]),])
df_Target_sig$group <- gsub("c1_Target_sig", "cluster 1", df_Target_sig$Var2)
df_Target_sig$group <- gsub("c2_Target_sig", "cluster 2", df_Target_sig$group)
levels(df_Target_sig$group) <- c("cluster 1", "cluster 2")

colnames(df_Target_sig) <- c("Target", "cluster", "count", "group")

### 4. バープロットを作成 ######################################
g_Target <-
  ggplot(df_Target_sig, aes(x = Target, y = count, fill = group)) +  # x軸にギア、y軸に数を設定
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフの高さをデータの値に基づいて設定
  labs(#title = "Significant Target Counts by cluster", 
       x = "Drug target", y = "Count",  fill = "Cluster") +
  theme_bw() 

### 5. バープロットを調整（横） ######################################
g_Target_h <- g_Target  + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                                legend.position = c(0.95, 0.95),  # 凡例をプロット内の右上に配置
                                legend.justification = c("right", "top"))  # 凡例のアンカーポイントを右上に設定
g_Target_h

### 5. バープロットを調整（縦） ######################################
g_Target_v <- g_Target + coord_flip()  + 
  theme( legend.position = c(1, 0.05),  # 凡例をプロット内の右上に配置
         legend.justification = c("right", "bottom"))
g_Target_v

### 6. save ##########################################################
list_graph_enriched_target <- list(h = g_Target_h, v = g_Target_v, n= g_Target)
save(list_graph_enriched_target, file = "./Analysis_240819/6.Enrichment_analysis/list_graph_enriched_target.dat")

