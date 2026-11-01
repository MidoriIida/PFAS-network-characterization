##################################################
#                                                #
#   PFAS-associated gene data preparation        #
#                                                #
#   Data sources:                                #
#     - CTD (Comparative Toxicogenomics DB)      #
#       http://ctdbase.org/                      #
#       File: CTD_chem_gene_ixns.csv             #
#       Human protein interactions filtered      #
#       (Organism = "Homo sapiens",              #
#        GeneForms = "protein")                  #
#       PFAS identified by cross-referencing     #
#       EPA PFAS master list (pfasmaster-        #
#       2023-01-16) via chemical name and CAS    #
#                                                #
#     - LINCS L1000 (level 5 signatures)         #
#       https://lincsproject.org/                #
#       File: L1000_level5_signatures_978genes   #
#       Top/bottom 5% of differentially          #
#       expressed landmark genes per PFAS        #
#       (matched by InChIKey to EPA list)        #
#                                                #
#   Output: list_PFAS_genes (list)               #
#     $Merged : data frame, one row per          #
#               PFAS-gene pair, with columns     #
#               CASRN / GeneID / PREFERRED.NAME  #
#               / ChemicalName / ChemicalID      #
#                                                #
#   Reads: ../../data/Supplementary_Data_2_      #
#           PFASgenes.txt                        #
#   Saves: list_PFAS_genes.dat                   #
#                                                #
##################################################

### 1. Load Supplementary Data 2 ######################################
# Supplementary_Data_2_PFASgenes.txt is space-delimited with no header.
# Columns: CASRN  GeneID  PREFERRED.NAME  ChemicalName  ChemicalID
# Note: PREFERRED.NAME may contain spaces; ChemicalName and ChemicalID
#       never contain spaces (parsed from the end of each line).

in_f  <- "../../data/Supplementary_Data_2_PFASgenes.txt"
out_f <- "./path/to/intermediate/data/"

lines <- readLines(in_f)

df_PFAS_genes <- do.call(rbind, lapply(lines, function(l) {
  tok <- strsplit(l, " ")[[1]]
  data.frame(
    CASRN          = tok[1],
    GeneID         = as.integer(tok[2]),
    PREFERRED.NAME = paste(tok[3:(length(tok) - 2)], collapse = " "),
    ChemicalName   = tok[length(tok) - 1],
    ChemicalID     = tok[length(tok)],
    stringsAsFactors = FALSE
  )
}))

cat("PFAS-gene pairs loaded:", nrow(df_PFAS_genes), "\n")
cat("Unique PFAS (CASRN):  ", length(unique(df_PFAS_genes$CASRN)), "\n")
cat("Unique genes (GeneID):", length(unique(df_PFAS_genes$GeneID)), "\n")

### 2. Build list_PFAS_genes for downstream analyses ##################
# Downstream scripts (1_drug_space, 2_disease_space) access
# list_PFAS_genes$Merged to obtain per-PFAS gene sets.

list_PFAS_genes <- list(Merged = df_PFAS_genes)

### 3. Save ###########################################################
save(list_PFAS_genes,
     file = paste0(out_f, "list_PFAS_genes.dat"))
