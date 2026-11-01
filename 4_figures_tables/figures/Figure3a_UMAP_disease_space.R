##########################################
#   Plot 45 PFAS into 148 disease space  #
#   Bubble size = nF (number of fluorine atoms)  #
#   lemon-free version                           #
##########################################
set.seed(224)   # same seed as original script for reproducibility

library(ggplot2)   # load before lemon to prevent +.gg override
library(dplyr)
library(umap)
library(openxlsx)
# library(lemon)  # excluded

# ============================================================
# 1. Load data
# ============================================================

## Disease-space UMAP results
load("/path/to/your/local/directory/list_umap_45PFAS_148disease.dat")

## PFAS abbreviations
in_f  <- "/path/to/your/local/directory/PFOA_abbrivation_list_pfasmaster-2023-01-31.xlsx"
abbri   <- read.xlsx(in_f)
Acronym <- setNames(abbri$Acronym, abbri$CASRN)

## Number of fluorine atoms (nF)
in_d  <- "/path/to/your/local/directory/original_data"
in_f3 <- "hsPPARalpha_PFASs_wash_descriptors_Sscore_IC50_data.xlsx"
discriptors <- read.xlsx(paste0(in_d, in_f3))
nF <- setNames(discriptors$a_nF, discriptors$CASRN)

# ============================================================
# 2. Hierarchical clustering (k = 2)
# ============================================================

hc       <- hclust(dist(list_umap$distnon$layout), method = "complete")
clusters <- cutree(hc, k = 2)

cat("Cluster assignment counts:\n")
print(table(clusters))

# ============================================================
# 3. Create data frame for plotting
# ============================================================

umap_df          <- as.data.frame(list_umap$distnon$layout)
umap_df$cluster  <- as.factor(clusters)
umap_df$nF       <- nF[rownames(umap_df)]
# Impute nF if missing (same values as original script)
umap_df[is.na(umap_df$nF), "nF"] <- c(5, 6, 9, 5, 4, 6, 5, 6, 17, 4, 6, 6)
umap_df$acronym  <- Acronym[rownames(umap_df)]
umap_df$accolors <- ifelse(umap_df$cluster == 1, "#577a00", "#9159ba")

# ============================================================
# 4. Visualization
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
# 5. Save
# ============================================================

out_dir <- "/path/to/your/local/directory/Code_for_Manuscript"
ggsave(file.path(out_dir, "Figure3a_UMAP_disease_space.pdf"),
       plot = p, width = 8, height = 6, units = "in", device = "pdf")
cat("Saved: Figure3a_UMAP_disease_space.pdf\n")
