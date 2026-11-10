##########################################
#   Figure 5                             #
#   4-group classification               #
#   Master list PFAS (n ~ 8,831)         #
#   y: Drug-space  (top=cl1, bot=cl2)   #
#   x: Disease-space (left=cl1, right=cl2)#
#                                        #
#   Group 1 (drug cl1 × dis cl1)        #
#   Group 2 (drug cl1 × dis cl2)        #
#   Group 3 (drug cl2 × dis cl1)        #
#   Group 4 (drug cl2 × dis cl2)        #
##########################################
library(ggplot2)
library(dplyr)

# ============================================================
# 1. Load data
# ============================================================

df_M <- read.csv("./path/to/intermediate/data/Prediction_pfas_masterlist.csv",
                 header = FALSE)
colnames(df_M) <- c("label", "SMILES", "CASRN",
                    "pred_MOA", "prob_cl1_MOA", "prob_cl2_MOA", "c6", "c7", "c8")

df_O <- read.csv("./path/to/intermediate/data/final_predictions_XGBoost_Morgan_chems_PFASmastelist.csv")
df_O$pred_outcome <- df_O$predicted_label + 1   # 0/1 → 1/2

# ============================================================
# 2. Merge and filter for high-confidence predictions
# ============================================================

df <- merge(df_M[, c("CASRN", "pred_MOA", "prob_cl1_MOA")],
            df_O[, c("CASRN", "pred_outcome", "probability")],
            by = "CASRN")

df_h <- df[(df$prob_cl1_MOA > 0.7 | df$prob_cl1_MOA < 0.3) &
           (df$probability   > 0.7 | df$probability   < 0.3), ]

cat("High-confidence PFAS n =", nrow(df_h), "\n")

# ============================================================
# 3. Define groups (same naming as Figure 4)
#   pred_MOA:    1=drug cl1, 2=drug cl2
#   pred_outcome:1=dis  cl1, 2=dis  cl2
# ============================================================

df_h <- df_h %>%
  mutate(group = case_when(
    pred_MOA == 1 & pred_outcome == 1 ~ "Group 1",
    pred_MOA == 1 & pred_outcome == 2 ~ "Group 2",
    pred_MOA == 2 & pred_outcome == 1 ~ "Group 3",
    pred_MOA == 2 & pred_outcome == 2 ~ "Group 4"
  ))

counts <- df_h %>%
  count(pred_MOA, pred_outcome, group) %>%
  mutate(
    # x-axis: disease-space (cl1=left, cl2=right)
    dis_label  = factor(pred_outcome,
                        levels = c(1, 2),
                        labels = c("Disease-space\ncluster 1",
                                   "Disease-space\ncluster 2")),
    # y-axis: drug-space (cl1=top, cl2=bottom) -> reverse factor levels
    drug_label = factor(pred_MOA,
                        levels = c(2, 1),   # putting 2 first: bottom=cl2, top=cl1
                        labels = c("Drug-space\ncluster 2",
                                   "Drug-space\ncluster 1")),
    label_text = paste0(group, "\n(n = ", formatC(n, big.mark=","), ")")
  )

print(counts[, c("group", "n")])

# ============================================================
# 4. Compute mosaic plot coordinates
#    x-axis: disease-space (left=cl1, right=cl2), width proportional to column total n
#    y-axis: drug-space (top=cl1, bottom=cl2), height proportional to within-column proportion
# ============================================================

total <- sum(counts$n)

dis1_total <- sum(counts$n[counts$pred_outcome == 1])
dis2_total <- sum(counts$n[counts$pred_outcome == 2])

# Column boundaries (x): dis cl1 = 0 ~ w1, dis cl2 = w1 ~ 1
w1 <- dis1_total / total

mosaic_df <- counts %>%
  rowwise() %>%
  mutate(
    xmin = ifelse(pred_outcome == 1, 0,  w1),
    xmax = ifelse(pred_outcome == 1, w1, 1),
    # drug cl1 = upper portion (y: h_lo ~ 1), drug cl2 = lower (y: 0 ~ h_lo)
    h_cl1 = ifelse(pred_outcome == 1,
                   counts$n[counts$pred_MOA==1 & counts$pred_outcome==1] / dis1_total,
                   counts$n[counts$pred_MOA==1 & counts$pred_outcome==2] / dis2_total),
    ymin  = ifelse(pred_MOA == 1, 1 - h_cl1, 0),
    ymax  = ifelse(pred_MOA == 1, 1,          1 - h_cl1),
    xmid  = (xmin + xmax) / 2,
    ymid  = (ymin + ymax) / 2,
    pct   = round(n / total * 100, 1),
    label_text = paste0(group, "\n(n = ", formatC(n, big.mark=",", format="d"), ")")
  ) %>%
  ungroup()

# ============================================================
# 5. Visualization
# ============================================================

# Axis label positions
x_breaks <- c(w1 / 2, w1 + (1 - w1) / 2)
y_breaks <- c(
  (counts$n[counts$pred_MOA==2 & counts$pred_outcome==1] / dis1_total) / 2,   # Group 3 center
  1 - (counts$n[counts$pred_MOA==1 & counts$pred_outcome==1] / dis1_total) / 2 # Group 1 center
)

p <-
  ggplot(mosaic_df) +
  geom_rect(aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax,
                fill = n),
            color = "white", linewidth = 2) +
  geom_text(aes(x = xmid, y = ymid, label = label_text),
            size = 5, fontface = "bold",
            color = ifelse(mosaic_df$n > 2000, "white", "black")) +
  scale_fill_gradient(low = "grey90", high = "grey30") +
  scale_x_continuous(
    breaks = x_breaks,
    labels = c("Disease-space\ncluster 1", "Disease-space\ncluster 2"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = c("Drug-space\ncluster 2", "Drug-space\ncluster 1"),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 12)
  )
p

# ============================================================
# 5. Save
# ============================================================

out_dir <- "/path/to/your/local/directory/Code_for_Manuscript"
ggsave(file.path(out_dir, "Figure5_four_groups_masterlist.pdf"),
       plot = p, width = 6, height = 6, units = "in", device = "pdf")
cat("Saved: Figure5_four_groups_masterlist.pdf\n")
