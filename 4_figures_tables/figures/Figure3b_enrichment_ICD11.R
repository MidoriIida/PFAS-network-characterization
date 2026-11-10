##########################################
#   Figure 3b                            #
#   ICD-11 enrichment bar plot           #
#   cluster 1 vs cluster 2              #
#   lemon-free version                           #
##########################################
set.seed(224)

library(ggplot2)   # load before lemon
library(dplyr)
library(umap)
library(reshape2)
library(stringr)
# library(lemon)  # excluded

# ============================================================
# 1. Load data
# ============================================================

## ICD-11 enrichment results (rows = PFAS, cols = ICD-11 chapters)
load("./path/to/intermediate/data/df_results_enrichment_KEGG_ICD11.dat")

## UMAP + clustering (same settings as Figure 3a)
load("./path/to/intermediate/data/list_umap_45PFAS_148disease.dat")
hc       <- hclust(dist(list_umap$distnon$layout), method = "complete")
clusters <- cutree(hc, k = 2)
umap_df  <- as.data.frame(list_umap$distnon$layout)
umap_df$cluster <- as.factor(clusters)

# ============================================================
# 2. Count significant ICD-11 enrichments per cluster
# ============================================================

cl1 <- rownames(umap_df[umap_df$cluster == 1, ])
cl2 <- rownames(umap_df[umap_df$cluster == 2, ])

cat("cluster 1:", length(cl1), "PFAS\n")
cat("cluster 2:", length(cl2), "PFAS\n")

# Count PFAS significantly enriched (p < 0.05) per ICD-11 category
cluster_1_count <- colSums(df_results[cl1, ] < 0.05, na.rm = TRUE)
cluster_2_count <- colSums(df_results[cl2, ] < 0.05, na.rm = TRUE)

# Normalize by cluster size (proportion: %)
# Use proportions because cluster 1 and 2 have different sizes
cluster_1_pct <- cluster_1_count / length(cl1) * 100
cluster_2_pct <- cluster_2_count / length(cl2) * 100

cat("cluster 1 size:", length(cl1), "/ cluster 2 size:", length(cl2), "\n")

# Create data frame
df_plot <- data.frame(
  code    = rep(names(cluster_1_pct), 2),
  count   = c(cluster_1_pct, cluster_2_pct),
  group   = rep(c("cluster 1", "cluster 2"),
                c(length(cluster_1_pct), length(cluster_2_pct))),
  group_draw_order = factor(
    rep(c("cluster 1", "cluster 2"),
        c(length(cluster_1_pct), length(cluster_2_pct))),
    levels = c("cluster 2", "cluster 1")  # cluster 1 in foreground
  ),
  code_wrapped = str_wrap(rep(names(cluster_1_pct), 2), width = 40)
)

# Exclude categories with zero count in both clusters
df_plot <- df_plot[df_plot$code %in%
                     names(which(tapply(df_plot$count, df_plot$code, sum) > 0)), ]

cat("Number of ICD-11 categories to plot:", length(unique(df_plot$code)), "\n")

# ============================================================
# 3. Visualization
# ============================================================

g_b <-
  ggplot(df_plot, aes(x = code_wrapped, y = count, fill = group)) +
  geom_bar(stat = "identity", position = "dodge",
           aes(group = group_draw_order)) +
  scale_fill_manual(values = c("cluster 1" = "#7CAE00",
                               "cluster 2" = "#C77CFF"),
                    name = "Cluster") +
  scale_y_continuous(breaks = function(x) seq(0, ceiling(max(x)), by = 10),
                     minor_breaks = NULL,
                     limits = c(0, NA)) +
  labs(x = "ICD-11 chapter", y = "Percentage of PFAS with \n significant enrichment (%)") +
  theme_bw() +
  coord_flip() +
  theme(
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 9),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text  = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "none"   # no legend (consistent with Figure 3c; output separately)
  )
g_b

# ============================================================
# 4. Save
# ============================================================

out_dir <- "/path/to/your/local/directory/Code_for_Manuscript"
ggsave(file.path(out_dir, "Figure3b_enrichment_ICD11.pdf"),
       plot = g_b, width = 5, height = 8, units = "in", device = "pdf")
cat("Saved: Figure3b_enrichment_ICD11.pdf\n")
