##########################################
#   Figure S9a                           #
#   Drug-space classifier                #
#   Prediction probability histogram     #
#   for PFAS masterlist (n = 10,974)     #
##########################################

library(ggplot2)
library(dplyr)

# ============================================================
# 1. Load data
# ============================================================

df_pred <- read.csv(
  "/path/to/your/local/directory/Prediction_pfas_masterlist.csv",
  header = FALSE
)
df_pred <- df_pred[, 1:6]
colnames(df_pred) <- c("url", "SMILES", "CASRN",
                       "prediction_label", "prob_cls1", "prob_cls2")

cat("Total PFAS:", nrow(df_pred), "\n")
cat("Cluster 1:", sum(df_pred$prediction_label == 1), "\n")
cat("Cluster 2:", sum(df_pred$prediction_label == 2), "\n")
cat("Boundary (0.3-0.7):",
    sum(df_pred$prob_cls1 > 0.3 & df_pred$prob_cls1 < 0.7), "\n")

# ============================================================
# 2. Visualization
# ============================================================

p <-
  ggplot(df_pred, aes(x = prob_cls1)) +
  geom_histogram(binwidth = 0.02, fill = "grey60", color = "white") +
  geom_vline(xintercept = 0.3, linetype = "dashed",
             color = "firebrick", linewidth = 0.7) +
  geom_vline(xintercept = 0.7, linetype = "dashed",
             color = "firebrick", linewidth = 0.7) +
  annotate("rect", xmin = 0.3, xmax = 0.7, ymin = -Inf, ymax = Inf,
           fill = "firebrick", alpha = 0.07) +
  annotate("text", x = 0.5, y = Inf,
           label = paste0("Uncertainty region\n(n = ",
                          sum(df_pred$prob_cls1 > 0.3 & df_pred$prob_cls1 < 0.7),
                          ")"),
           vjust = 1.3, size = 3.5, color = "firebrick") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  labs(
    x = "Prediction probability for cluster 1",
    y = "Number of PFAS"
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
p

# ============================================================
# 3. Save
# ============================================================

out_dir <- "/path/to/your/local/directory/Supplymentary_Figures"
ggsave(file.path(out_dir, "FigureS9a_drug_space_probability_histogram.pdf"),
       plot = p, width = 6, height = 4, units = "in", device = "pdf")
cat("Saved: FigureS9a_drug_space_probability_histogram.pdf\n")
