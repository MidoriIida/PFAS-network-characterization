############################
#  Supplementary Data 3    #
#        drug genes        #
#  900 drugs x 21,441 genes#
#  271,419 interactions    #
############################

### 1. Load ######################
# df_ctd_drugs produced by 0_data_preparation/3.drug_data_load.R
in_f <- "./path/to/intermediate/data/df_ctd_drugs.dat"
load(in_f)

cat("Interactions:", nrow(df_ctd_drugs), "\n")                          # 271,419
cat("Unique drugs (ChemicalID):", length(unique(df_ctd_drugs$ChemicalID)), "\n")  # 900
cat("Unique genes (GeneID):    ", length(unique(df_ctd_drugs$GeneID)), "\n")      # 21,441

### 2. Save ######################
write.table(df_ctd_drugs[, c("ChemicalName", "ChemicalID", "CasRN",
                              "GeneSymbol", "GeneID")],
            file = "./Supplementary_Data_3_ApprovedDruggenes.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

cat("Saved: Supplementary_Data_3_ApprovedDruggenes.txt\n")
