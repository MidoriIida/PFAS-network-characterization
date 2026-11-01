##########################################
#   Plot according to 900 drug space     #
#              45 PFAS                   #
#      クラスターの意味 (ATC)            #
#   ※ lemon非依存版 → .datを作り直す   #
##########################################
set.seed(123)

library(ggplot2)   # ggplot2を最初にロード（lemonによる+.gg上書きを防ぐ）
library(dplyr)
library(reshape2)
library(cowplot)
library(stringr)
# library(lemon)  ← 使用していないため除外
# library(umap)   ← 使用していないため除外

### 1. load ###################################################################
load("/path/to/your/local/directory/list_cluster.dat")

in_d2 <- "/path/to/your/local/directory/orignail_data"
in_f2 <- "/20230830/result_enrichment_ATC_code_KEGGdrugBackground_20230825.csv"
ATC   <- read.csv(paste0(in_d2, in_f2), row.names = 1)

### 2. cluster ################################################################
clusters <- list_cluster$clusters

c1 <- names(clusters[clusters == 1])
c2 <- names(clusters[clusters == 2])

cas <- sapply(strsplit(gsub("^X", "", colnames(ATC)), "\\.\\."),
              function(x) { gsub("\\.", "-", x[1]) })
names(cas) <- cas

### 3. 各クラスターで有意なATCの集計 ##########################################
# cas は colnames(ATC) 由来 → 列（PFAS）をサブセットする
c1_ATC <- ATC[, cas %in% c1]
c2_ATC <- ATC[, cas %in% c2]

cat("cluster 1 size:", length(c1), "/ cluster 2 size:", length(c2), "\n")

# apply(..., 1, ...) = 行ごと（ATCカテゴリごと）に有意なPFAS数をカウント
c1_ATC_sig_raw <- apply(c1_ATC, 1, function(x) { sum(na.omit(x) < 0.05) })
c2_ATC_sig_raw <- apply(c2_ATC, 1, function(x) { sum(na.omit(x) < 0.05) })

# クラスターサイズで正規化（割合: %）
c1_ATC_sig <- c1_ATC_sig_raw / length(c1) * 100
c2_ATC_sig <- c2_ATC_sig_raw / length(c2) * 100

ATC_sig <- cbind(c1_ATC_sig, c2_ATC_sig)
y <- as.data.frame(ATC_sig[which(rowSums(ATC_sig) != 0), ]) %>%
  arrange(c1_ATC_sig - c2_ATC_sig)  # 上: cluster1-cluster2の差が大 / 下: cluster2優位

df_ATC_sig       <- melt(as.matrix(y))
df_ATC_sig$Label <- str_wrap(df_ATC_sig[, 1], width = 30)
df_ATC_sig$group <- gsub("c1_ATC_sig", "cluster 1", df_ATC_sig$Var2)
df_ATC_sig$group <- gsub("c2_ATC_sig", "cluster 2", df_ATC_sig$group)
colnames(df_ATC_sig) <- c("ATC", "cluster", "count", "Label", "group")

cat("集計完了: ATC数 =", length(unique(df_ATC_sig$ATC)), "\n")

### 4. バープロット作成 #########################################################
g_ATC <-
  ggplot(df_ATC_sig, aes(x = ATC, y = count, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(breaks = function(x) seq(0, ceiling(max(x)), by = 10),
                     minor_breaks = NULL,
                     limits = c(0, NA)) +
  labs(x = "ATC class", y = "Percentage of PFAS with significant enrichment (%)", fill = "Cluster") +
  theme_bw()

# 横向き（x軸ラベルを回転）
g_ATC_h <- g_ATC +
  scale_x_discrete(labels = df_ATC_sig$Label) +
  theme(axis.text.x         = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position      = "none")

# 縦向き（coord_flip）
g_ATC_v <- g_ATC +
  coord_flip() + 
  scale_x_discrete(labels = df_ATC_sig$Label) +
  theme(legend.position      = c(0.95, 0.95),
        legend.justification = c("right", "top"))

### 5. save ##################################################################
ggsave("Figure2b_enrichment_ATC.pdf", plot = g_ATC_v, width = 8, height = 12, units = "in", device = "pdf")
