##########################################
#   Plot according to 900 drug space     #
#              45 PFAS                   #
#   Cluster interpretation (ATC)          #
##########################################
set.seed(123)  # for reproducibility
library(ggplot2)
library(dplyr)
library(umap)
library(dplyr)
library(lemon)
library(reshape2)
library(cowplot)
library(stringr)

### 1. load ###################################################################
load("./path/to/intermediate/data/list_cluster.dat")

# ATC enrichment p-values (computed by collaborators; Supplementary Data 7)
ATC <- read.csv("../../data/Supplementary_Data_7_Enrichment_ATC.txt", row.names = 1)

# FDR correction
p_matrix <- as.matrix(ATC)
p_fdr <- matrix(
  p.adjust(as.vector(p_matrix), method = "fdr"),
  nrow = nrow(p_matrix),
  ncol = ncol(p_matrix)
)
rownames(p_fdr) <- rownames(p_matrix)
colnames(p_fdr) <- colnames(p_matrix)

### 2. cluster ################################################################
clusters <- list_cluster$clusters

c1 <- names(clusters[clusters == 1])
c2 <- names(clusters[clusters == 2])

make_ATC_plot <- function(mat) {
  cas <- sapply(strsplit(gsub("^X", "", colnames(mat)), "\\.\\."), function(x){gsub("\\.", "-", x[1])})
  names(cas) <- cas

  c1_ATC <- mat[cas %in% c1,]
  c2_ATC <- mat[cas %in% c2,]

  c1_ATC_sig <- apply(c1_ATC, 1, function(x){sum(na.omit(x) < 0.05)})
  c2_ATC_sig <- apply(c2_ATC, 1, function(x){sum(na.omit(x) < 0.05)})

  ATC_sig <- cbind(c1_ATC_sig, c2_ATC_sig)
  y <- as.data.frame(ATC_sig[which(rowSums(ATC_sig) != 0),]) %>%
    arrange(c2_ATC_sig, desc(c1_ATC_sig))
  df_ATC_sig <- melt(as.matrix(y))
  df_ATC_sig$Label <- str_wrap(df_ATC_sig[,1], width = 30)
  df_ATC_sig$group <- gsub("c1_ATC_sig", "cluster 1", df_ATC_sig$Var2)
  df_ATC_sig$group <- gsub("c2_ATC_sig", "cluster 2", df_ATC_sig$group)
  levels(df_ATC_sig$group) <- c("cluster 1", "cluster 2")
  colnames(df_ATC_sig) <- c("ATC", "cluster", "count", "Label", "group")

  g_ATC <-
    ggplot(df_ATC_sig, aes(x = ATC, y = count, fill = group)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "ATC class", y = "Count", fill = "Cluster") +
    theme_bw()

  g_ATC_h <- g_ATC + scale_x_discrete(labels = df_ATC_sig$Label) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = c(0.95, 0.95),
          legend.justification = c("right", "top"))

  g_ATC_v <- g_ATC + coord_flip() + scale_x_discrete(labels = df_ATC_sig$Label)

  list(h = g_ATC_h, v = g_ATC_v, n = g_ATC)
}

### 3. Raw p-values ##########################################################
list_graph_enriched_ATC <- make_ATC_plot(ATC)
list_graph_enriched_ATC$h
save(list_graph_enriched_ATC, file = "./path/to/intermediate/data/6.Enrichment_analysis/list_graph_enriched_ATC.dat")

### 4. FDR-corrected p-values ################################################
list_graph_enriched_ATC_fdr <- make_ATC_plot(p_fdr)
list_graph_enriched_ATC_fdr$h
save(list_graph_enriched_ATC_fdr, file = "./path/to/intermediate/data/6.Enrichment_analysis/list_graph_enriched_ATC_fdr.dat")
