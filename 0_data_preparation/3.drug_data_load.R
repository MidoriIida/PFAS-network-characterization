##################################################
#                                                #
#   Approved drug-gene interaction data          #
#                                                #
#   Data source:                                 #
#     CTD (Comparative Toxicogenomics Database)  #
#     http://ctdbase.org/                        #
#     File: CTD_chem_gene_ixns.csv               #
#     Filtered to:                               #
#       - Homo sapiens (OrganismID = 9606)       #
#       - Gene form = protein                    #
#       - 900 approved drugs used in the         #
#         network proximity analysis             #
#         (drug list from Cheng et al. 2019      #
#          Nat Commun; CAS/name-matched to CTD)  #
#                                                #
#   Input:  ../../data/Supplementary_Data_3_     #
#            ApprovedDruggenes.txt               #
#   Output: df_ctd_drugs.dat                     #
#     data frame: ChemicalName / ChemicalID /    #
#     CasRN / GeneSymbol / GeneID               #
#                                                #
##################################################
library(data.table)

### 1. Load Supplementary Data 3 ######################################
in_f  <- "../../data/Supplementary_Data_3_ApprovedDruggenes.txt"
out_f <- "./path/to/intermediate/data/"

df_ctd_drugs <- as.data.frame(
  fread(in_f, sep = "\t", header = TRUE, colClasses = "character")
)
df_ctd_drugs$GeneID <- as.integer(df_ctd_drugs$GeneID)

cat("Drug-gene interactions loaded:", nrow(df_ctd_drugs), "\n")
cat("Unique drugs (ChemicalID):    ", length(unique(df_ctd_drugs$ChemicalID)), "\n")
cat("Unique genes (GeneID):        ", length(unique(df_ctd_drugs$GeneID)), "\n")

### 2. Save ###########################################################
save(df_ctd_drugs, file = paste0(out_f, "df_ctd_drugs.dat"))
