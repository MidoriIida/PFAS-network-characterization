##################################################
#                                                #
#   Disease-gene association data                #
#                                                #
#   Data sources:                                #
#     OMIM (Online Mendelian Inheritance in Man) #
#       https://www.omim.org/                    #
#     GWAS Catalog                               #
#       https://www.ebi.ac.uk/gwas/              #
#     ClinVar                                    #
#       https://www.ncbi.nlm.nih.gov/clinvar/    #
#     PheWAS Catalog                             #
#       https://phewascatalog.org/               #
#     Gene IDs converted via Ensembl BioMart     #
#       (biomaRt R package)                      #
#                                                #
#   152 diseases x EntrezGeneIDs                 #
#                                                #
#   Input:  ../../data/Supplementary_Data_4_     #
#            DiseaseGenes.txt                    #
#   Output: list_disease_genes.dat               #
#     named list: disease name -> integer vector #
#     of EntrezGeneIDs                           #
#                                                #
##################################################

### 1. Load Supplementary Data 4 ######################################
in_f  <- "../../data/Supplementary_Data_4_DiseaseGenes.txt"
out_f <- "./path/to/intermediate/data/"

lines <- readLines(in_f)
lines <- lines[-1]  # remove header ("DiseaseNames EntrezGeneID")

df_disease_genes <- do.call(rbind, lapply(lines, function(l) {
  tok <- strsplit(l, " ")[[1]]
  data.frame(
    DiseaseNames = paste(tok[seq_len(length(tok) - 1)], collapse = " "),
    EntrezGeneID = as.integer(tok[length(tok)]),
    stringsAsFactors = FALSE
  )
}))

### 2. Reconstruct named list ##########################################
list_disease_genes <- tapply(df_disease_genes$EntrezGeneID,
                             df_disease_genes$DiseaseNames,
                             function(x) as.integer(unique(x)))
list_disease_genes <- as.list(list_disease_genes)

cat("Unique diseases: ", length(list_disease_genes), "\n")
cat("Total gene-disease pairs:", nrow(df_disease_genes), "\n")

### 3. Save ###########################################################
save(list_disease_genes, file = paste0(out_f, "list_disease_genes.dat"))
