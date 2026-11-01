############################
#  Supplementary Data 3    #
#        drug genes        #
#  900 drugs × 21,441 genes#
#  271,419 interactions    #
############################
library(stringr)

### 1. load ######################

# 900薬物リスト（Takahashiラボ提供の近接性解析データから）
in_drugs <- "~/Documents/Lab_others/Lab_Takemoto/2023_Takahashi/orignail_data/20230703/summary_PFAS_appdrug_new_20230703.txt"
PFAS_drugs <- read.table(in_drugs, header=FALSE, sep="\t",
                         comment.char="", quote='"')
drug_names <- unique(PFAS_drugs$V2)
cat("ユニーク薬物数:", length(drug_names), "\n")  # 900

# CTD (accessed Jun 29, 2023)
in_ctd <- "/Users/iida/Documents/Lab_Iida/Orignal_data/CTD/CTD_chem_gene_ixns.csv"
ctd <- read.csv(in_ctd, skip=29, header=FALSE)
colnames(ctd) <- c("ChemicalName", "ChemicalID", "CasRN", "GeneSymbol",
                   "GeneID", "GeneForms", "Organism", "OrganismID",
                   "Interaction", "InteractionActions", "PubMedIDs")

### 2. フィルタリング ######################
# Homo sapiens のみ
ctd_human <- ctd[ctd$OrganismID == "9606", ]

# 900薬物のみ抽出
df_ctd_drugs <- ctd_human[ctd_human$ChemicalName %in% drug_names, ]

cat("interactions:", nrow(df_ctd_drugs), "\n")         # 271,419
cat("薬物数:", length(unique(df_ctd_drugs$ChemicalID)), "\n")  # 900
cat("遺伝子数:", length(unique(df_ctd_drugs$GeneID)), "\n")    # 21,441

### 3. 保存 ######################
out_df <- df_ctd_drugs[, c("ChemicalName", "ChemicalID", "CasRN",
                            "GeneSymbol", "GeneID")]

write.table(out_df,
            file = "./Supplementary_Data_3_ApprovedDruggenes_tentative.txt",
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

cat("Saved: Supplementary_Data_3_ApprovedDruggenes_tentative.txt\n")
