##############################################################
#  Supplementary Data 9                                      #
#  4-group classification results for all EPA PFAS           #
#  Drug-space (RandomForest) × Disease-space (XGBoost)       #
#  n = 10,923 (all PFAS with valid SMILES in both models)    #
#  high_confidence = TRUE: prob > 0.7 or < 0.3 in both      #
#                          spaces (n = 8,831; used Fig. 5)   #
##############################################################
library(dplyr)

base <- "/path/to/your/local/directory/Human_PPI_PFAS"

# ============================================================
# 1. データ読み込み
# ============================================================

# Drug-space: RandomForest prediction on PFAS master list
# Columns (no header): DTXSID, SMILES, CASRN, pred_MOA, prob_cl1_MOA, prob_cl2_MOA, ...
df_M <- read.csv(file.path(base, "Analysis_250412/Prediction_pfas_masterlist.csv"),
                 header = FALSE)
colnames(df_M) <- c("DTXSID", "SMILES", "CASRN",
                    "pred_MOA", "prob_cl1_MOA", "prob_cl2_MOA",
                    "c7", "c8", "c9")

# Disease-space: XGBoost prediction on PFAS master list
# probability = P(cluster 2); pred_disease_cl1 = 1 - probability
df_O <- read.csv(file.path(base,
                 "Analysis_251010/final_predictions_XGBoost_Morgan_chems_PFASmastelist.csv"))
df_O$pred_disease_cluster <- df_O$predicted_label + 1   # 0/1 → 1/2
df_O$prob_disease_cl1     <- 1 - df_O$probability       # P(cluster 1)

# ============================================================
# 2. CASRN で結合 (merge → 8,831 PFAS)
# ============================================================
df <- merge(
  df_M[, c("CASRN", "SMILES", "pred_MOA", "prob_cl1_MOA")],
  df_O[, c("CASRN", "pred_disease_cluster", "prob_disease_cl1")],
  by = "CASRN"
)

cat("Total PFAS after merge:", nrow(df), "\n")  # 10,923

# ============================================================
# 3. グループ定義
#   Group 1: drug cl1 × disease cl1
#   Group 2: drug cl1 × disease cl2
#   Group 3: drug cl2 × disease cl1
#   Group 4: drug cl2 × disease cl2
# ============================================================
df <- df %>%
  mutate(
    group = case_when(
      pred_MOA == 1 & pred_disease_cluster == 1 ~ "Group 1",
      pred_MOA == 1 & pred_disease_cluster == 2 ~ "Group 2",
      pred_MOA == 2 & pred_disease_cluster == 1 ~ "Group 3",
      pred_MOA == 2 & pred_disease_cluster == 2 ~ "Group 4"
    ),
    # high_confidence: prob > 0.7 or < 0.3 in both spaces (same filter as Figure 5)
    high_confidence = (prob_cl1_MOA > 0.7 | prob_cl1_MOA < 0.3) &
                      (prob_disease_cl1 > 0.7 | prob_disease_cl1 < 0.3)
  )

cat("Group counts (all):\n")
print(table(df$group))
cat("\nGroup counts (high-confidence, n =", sum(df$high_confidence), "):\n")
print(table(df$group[df$high_confidence]))

# ============================================================
# 4. 出力列を整理・保存
# ============================================================
out_df <- df %>%
  select(
    CASRN,
    SMILES,
    pred_drug_cluster    = pred_MOA,
    prob_drug_cl1        = prob_cl1_MOA,
    pred_disease_cluster,
    prob_disease_cl1,
    group,
    high_confidence
  ) %>%
  arrange(group, CASRN)

out_path <- "/path/to/your/local/directory/Supplementary_Data_9_Classification_results_for_all_EPA_PFAS.txt"
write.table(out_df, file = out_path,
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

cat("Saved:", out_path, "\n")
cat("Total rows:", nrow(out_df), "\n")
cat("High-confidence rows:", sum(out_df$high_confidence), "\n")
