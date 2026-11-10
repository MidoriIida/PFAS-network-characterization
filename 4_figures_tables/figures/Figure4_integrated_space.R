##########################################
#   Figure 4                             #
#   Integrated Drug-space &              #
#   Disease-space cluster scores         #
#   y: Drug-space cluster 1 score        #
#   x: Disease-space cluster 1 score     #
#      (x reversed: left=cl1, right=cl2) #
#                                        #
#   Group 1 (drug cl1 + dis cl1): n=15  #
#   Group 2 (drug cl1 + dis cl2): n=1   #
#   Group 3 (drug cl2 + dis cl1): n=10  #
#   Group 4 (drug cl2 + dis cl2): n=9   #
##########################################
set.seed(123)

library(ggplot2)
library(dplyr)
library(ggrepel)
library(openxlsx)

# ============================================================
# 1. Load data
# ============================================================

## Drug-space UMAP
load("./path/to/intermediate/data/list_umap_PFAS_drug.dat")
list_umap_drug <- list_umap

## Disease-space UMAP
load("./path/to/intermediate/data/list_umap_45PFAS_148disease.dat")
list_umap_dis <- list_umap

## Pre-computed clusters (hierarchical clustering not re-run)
drug_cl_dat  <- get(load("./path/to/intermediate/data/list_cluster.dat"))
dis_cl_dat   <- get(load("./path/to/intermediate/data/list_cluster_disease.dat"))

clusters_drug    <- drug_cl_dat$clusters    # named vector, 35 PFAS
clusters_disease <- dis_cl_dat$clusters     # named vector, 44 PFAS

cat("Drug-space cluster sizes:\n");    print(table(clusters_drug))
cat("Disease-space cluster sizes:\n"); print(table(clusters_disease))

## PFAS abbreviations
in_f  <- "../../data/PFOA_abbrivation_list_pfasmaster-2023-01-31.xlsx"
abbri <- read.xlsx(in_f)
abbri <- abbri[!is.na(abbri$Acronym), ]
abbri <- abbri[!duplicated(abbri$CASRN), ]

# ============================================================
# 2. Create data frame from UMAP coordinates + cluster assignments
# ============================================================

drug_umap_df <- as.data.frame(list_umap_drug$distnon$layout) %>%
  tibble::rownames_to_column("CAS") %>%
  rename(UMAP1_drug = V1, UMAP2_drug = V2) %>%
  mutate(cluster_drug = clusters_drug[CAS])

disease_umap_df <- as.data.frame(list_umap_dis$distnon$layout) %>%
  tibble::rownames_to_column("CAS") %>%
  rename(UMAP1_dis = V1, UMAP2_dis = V2) %>%
  mutate(cluster_disease = clusters_disease[CAS])

# ============================================================
# 3. Compute cluster centroids
# ============================================================

centroids_drug <- drug_umap_df %>%
  group_by(cluster_drug) %>%
  summarise(cx = mean(UMAP1_drug), cy = mean(UMAP2_drug), .groups = "drop")

centroids_disease <- disease_umap_df %>%
  group_by(cluster_disease) %>%
  summarise(cx = mean(UMAP1_dis), cy = mean(UMAP2_dis), .groups = "drop")

# ============================================================
# 4. Compute Cluster 1 score
#    score = d2 / (d1 + d2): closer to cluster 1 centroid -> score closer to 1
# ============================================================

drug_scored <- drug_umap_df %>%
  rowwise() %>%
  mutate(
    d1_drug = sqrt((UMAP1_drug - centroids_drug$cx[centroids_drug$cluster_drug == 1])^2 +
                   (UMAP2_drug - centroids_drug$cy[centroids_drug$cluster_drug == 1])^2),
    d2_drug = sqrt((UMAP1_drug - centroids_drug$cx[centroids_drug$cluster_drug == 2])^2 +
                   (UMAP2_drug - centroids_drug$cy[centroids_drug$cluster_drug == 2])^2),
    drug_c1_score = d2_drug / (d1_drug + d2_drug)
  ) %>%
  ungroup()

disease_scored <- disease_umap_df %>%
  rowwise() %>%
  mutate(
    d1_dis = sqrt((UMAP1_dis - centroids_disease$cx[centroids_disease$cluster_disease == 1])^2 +
                  (UMAP2_dis - centroids_disease$cy[centroids_disease$cluster_disease == 1])^2),
    d2_dis = sqrt((UMAP1_dis - centroids_disease$cx[centroids_disease$cluster_disease == 2])^2 +
                  (UMAP2_dis - centroids_disease$cy[centroids_disease$cluster_disease == 2])^2),
    disease_c1_score = d2_dis / (d1_dis + d2_dis)
  ) %>%
  ungroup()

# ============================================================
# 5. Create integrated data frame
# ============================================================

df_scored <- inner_join(
  drug_scored    %>% select(CAS, cluster_drug,    drug_c1_score),
  disease_scored %>% select(CAS, cluster_disease, disease_c1_score),
  by = "CAS"
) %>%
  left_join(abbri[, c("CASRN", "Acronym")], by = c("CAS" = "CASRN")) %>%
  mutate(
    # Manually fill labels when Acronym is NA
    Acronym = case_when(
      CAS == "112281-77-3" ~ "Tetraconazole",
      TRUE ~ Acronym
    ),
    # Group definition: drug_cluster-disease_cluster
    #   1-1 → Group 1 (n=15)
    #   1-2 -> Group 2 (n=1, ambiguous)
    #   2-1 → Group 3 (n=10)
    #   2-2 → Group 4 (n=9)
    group = case_when(
      cluster_drug == 1 & cluster_disease == 1 ~ "Group 1",
      cluster_drug == 1 & cluster_disease == 2 ~ "Group 2",
      cluster_drug == 2 & cluster_disease == 1 ~ "Group 3",
      cluster_drug == 2 & cluster_disease == 2 ~ "Group 4"
    )
  )

cat("\nGroup distribution:\n")
print(table(df_scored$group))

# ============================================================
# 6. Visualization
#    y-axis: drug-space cl1 score (top = drug cl1)
#    x-axis: disease-space cl1 score (left = dis cl1, right = dis cl2: reversed)
#    -> upper-left = high cl1 similarity (Group 1)
#    -> lower-right = low cl1 similarity (Group 4)
# ============================================================

# Boundary positions (midpoints of score gaps between groups)
# drug:    Group1/2 min=0.7118, Group3/4 max=0.4990 → midpoint
# disease: Group1/3 min=0.3991, Group2/4 max=0.3525 → midpoint
vline <- (0.3991 + 0.3525) / 2   # disease-space boundary approx 0.376 (reversed axis)
hline <- (0.7118 + 0.4990) / 2   # drug-space boundary approx 0.605

shape_values <- c("Group 1" = 15,   # ■
                  "Group 2" = 3,    # +
                  "Group 3" = 16,   # ●
                  "Group 4" = 17)   # ▲

p <-
  ggplot(df_scored, aes(x = disease_c1_score, y = drug_c1_score,
                        shape = group)) +
  geom_hline(yintercept = hline, linetype = "dashed", color = "gray75") +
  geom_vline(xintercept = vline, linetype = "dashed", color = "gray75") +
  geom_point(size = 4) +
  geom_text_repel(
    aes(label = Acronym),
    size          = 4,
    max.overlaps  = 20,
    box.padding   = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size  = 0.2
  ) +
  scale_shape_manual(values = shape_values) +
  scale_x_reverse(limits = c(1, 0)) +   # left = high dis cl1, right = high dis cl2
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x     = "Disease-space cluster 1 score",
    y     = "Drug-space cluster 1 score",
    shape = NULL
  ) +
  theme_bw() +
  theme(
    legend.position  = "right",
    axis.text.x      = element_text(size = 10),
    axis.text.y      = element_text(size = 10),
    axis.title.x     = element_text(size = 12),
    axis.title.y     = element_text(size = 12)
  )
p

# ============================================================
# 7. Save
# ============================================================

out_dir <- "/path/to/your/local/directory/Code_for_Manuscript"
ggsave(file.path(out_dir, "Figure4_integrated_space.pdf"),
       plot = p, width = 7, height = 7, units = "in", device = "pdf")
cat("Saved: Figure4_integrated_space.pdf\n")

write.csv(df_scored, file.path(out_dir, "Figure4_integrated_scores.csv"),
          row.names = FALSE, quote = FALSE)
cat("Saved: Figure4_integrated_scores.csv\n")
