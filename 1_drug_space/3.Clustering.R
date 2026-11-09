##########################################
#   Plot according to 900 drug space     #
#              45 PFAS                   #
#   Hierarchical clustering on UMAP      #
##########################################
set.seed(123)  # for reproducibility
library(ggplot2)
library(dplyr)

### 1. load ###################################################################
# list_umap produced by 1_drug_space/2.Calculate_45PFAS_on_900drug_space.R
load("./path/to/intermediate/data/list_umap_PFAS_drug.dat")

### 2. Hierarchical clustering on UMAP layout #################################
hc <- hclust(dist(list_umap$distnon$layout), method = "complete")

# Dendrogram
plot(hc)

# Cut into 2 clusters
clusters <- cutree(hc, k = 2)

list_cluster <- list(hc = hc, clusters = clusters)
save(list_cluster, file = "./path/to/intermediate/data/list_cluster.dat")

### 3. Plot UMAP with cluster colors ##########################################
umap_df         <- as.data.frame(list_umap$distnon$layout)
umap_df$cluster <- as.factor(clusters)

p <- ggplot(umap_df, aes(x = V1, y = V2, color = cluster)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_bw()
p
