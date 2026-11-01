############################
#  Supplementary Data 3    #
#        drug genes        #
#  900 drugs × 21,441 genes#
#  271,419 interactions    #
############################
library(stringr)

### 1. load ######################

# List of 900 approved drugs (from proximity analysis dataset)
in_drugs <- "/path/to/your/local/directory/summary_PFAS_appdrug_new_20230703.txt"
PFAS_drugs <- read.table(in_drugs, header=FALSE, sep="\t",
                         comment.char="", quote='"')
drug_names <- unique(PFAS_drugs$V2)
cat("Unique drug count:", length(drug_names), "\n")  # 900

# CTD (accessed Jun 29, 2023)
in_ctd <- "/path/to/your/local/directory/CTD_chem_gene_ixns.csv"
ctd <- read.csv(in_ctd, skip=29, header=FALSE)
colnames(ctd) <- c("ChemicalName", "ChemicalID", "CasRN", "GeneSymbol",
                   "GeneID", "GeneForms", "Organism", "OrganismID",
                   "Interaction", "InteractionActions", "PubMedIDs")

### 2. Filtering ######################
# Homo sapiens only
ctd_human <- ctd[ctd$OrganismID == "9606", ]

# Extract only the 900 drugs
df_ctd_drugs <- ctd_human[ctd_human$ChemicalName %in% drug_names, ]

cat("interactions:", nrow(df_ctd_drugs), "\n")         # 271,419
cat("Drug count:", length(unique(df_ctd_drugs$ChemicalID)), "\n")  # 900
cat("Gene count:", length(unique(df_ctd_drugs$GeneID)), "\n")    # 21,441

### 3. Save ######################
out_df <- df_ctd_drugs[, c("ChemicalName", "ChemicalID", "CasRN",
                            "GeneSymbol", "GeneID")]

write.table(out_df,
            file = "./Supplementary_Data_3_ApprovedDruggenes_tentative.txt",
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

cat("Saved: Supplementary_Data_3_ApprovedDruggenes_tentative.txt\n")
