##########################################
#   Plot 45 PFAS into 148 disease space  #
#   Bubble size = number of PFAS-associated genes  #
#   lemon-free version                           #
##########################################
set.seed(224)   # same seed as original script for reproducibility

library(ggplot2)   # load before lemon to prevent +.gg override
library(dplyr)
library(umap)
library(openxlsx)
library(tibble)
# library(lemon)  # excluded

# ============================================================
# 1. Load data
# ============================================================

## Disease-space UMAP results
load("./path/to/intermediate/data/list_umap_45PFAS_148disease.dat")

## PFAS abbreviations
in_f  <- "../../data/PFOA_abbrivation_list_pfasmaster-2023-01-31.xlsx"
abbri   <- read.xlsx(in_f)
Acronym <- setNames(abbri$Acronym, abbri$CASRN)

## PFAS-associated gene data
load("./path/to/intermediate/data/list_PFAS_genes.dat")

# ============================================================
# 2. Hierarchical clustering (k = 2)
# ============================================================

hc       <- hclust(dist(list_umap$distnon$layout), method = "complete")
clusters <- cutree(hc, k = 2)

cat("Cluster assignment counts:\n")
print(table(clusters))

# ============================================================
# 3. Compute number of PFAS-associated genes
# ============================================================

cas_disease <- rownames(list_umap$distnon$layout)   # CAS numbers for 44 PFAS

n_genes <- sapply(cas_disease, function(cas) {
  length(unique(list_PFAS_genes$Merged[list_PFAS_genes$Merged$CASRN == cas, "GeneID"]))
})

cat("\nGene count summary:\n")
print(summary(n_genes))

# ============================================================
# 4. Create data frame for plotting
# ============================================================

umap_df           <- as.data.frame(list_umap$distnon$layout)
umap_df$cluster   <- as.factor(clusters)
umap_df$n_genes   <- n_genes[rownames(umap_df)]
umap_df$log_genes <- log10(umap_df$n_genes + 1)   # log10 transform
umap_df$acronym   <- Acronym[rownames(umap_df)]
umap_df$accolors  <- ifelse(umap_df$cluster == 1, "#577a00", "#9159ba")  # text colors from original script

# ============================================================
# 5. Visualization
# ============================================================

p <-
  ggplot(umap_df, aes(x = V1, y = V2, color = cluster, size = log_genes)) +
  geom_point() +
  scale_size(range = c(2, 10),
             name  = "Associated\ngenes (log₁₀)") +
  scale_color_manual(values = c("#7CAE00", "#C77CFF"),   # colors from original script
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
# 6. Save
# ============================================================

out_dir <- "/path/to/your/local/directory/Code_for_Manuscript"
ggsave(file.path(out_dir, "FigureS5_UMAP_disease_space.pdf"),
       plot = p, width = 8, height = 6, units = "in", device = "pdf")
cat("Saved: FigureS5_UMAP_disease_space.pdf\n")
