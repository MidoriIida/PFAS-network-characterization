##########################################
#   Plot according to 900 drug space     #
#              45 PFAS                   #
#          Clustering                   #
##########################################
set.seed(123)  # for reproducibility
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

### 3. Hierarchical clustering #################################################################
# Apply k-means clustering to dimensionality-reduced data
hc <- hclust(dist(list_umap$distnon$layout), method = "complete")

# Plot dendrogram
plot(hc)

# Cut dendrogram into 2 clusters
clusters <- cutree(hc, k = 2)

#save
list_cluster <- list(hc = hc, clusters = clusters)

### 4. Color-coded plot of clustering results ###############################################
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
  theme_bw() + geom_text(aes(label = acronym), color = umap_df$accolors, size = 4, vjust = 0.5, hjust = 0.5) +  # Add labels to plot
  theme(    axis.text.x  = element_text(size = 10),  # x-axis text size
            axis.text.y  = element_text(size = 10),  # y-axis text size
            axis.title.x = element_text(size = 12),  # x-axis label size
            axis.title.y = element_text(size = 12),   # y-axis label size
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10) ) +
  guides(color = guide_legend(override.aes = list(size = 3))) # fix legend point size to 5
p

#ggsave("./path/to/intermediate/data/plot_UMAP.eps", plot = p, width = 24, height = 23, units = "in", device = "eps")
ggsave("Figure2a_UMAP_drug_space.pdf", plot = p, width = 8, height = 6, units = "in", device = "pdf")




