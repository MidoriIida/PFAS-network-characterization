###############################
#                             #
#  Disease genes data make    #
#                             #
###############################
library(dplyr)


### 1. load #############################
load(file = "./path/to/intermediate/data/0.make_data/causal_genes/list_disease_EntrezIDs.dat")
load(file = "./path/to/intermediate/data/0.make_data/causal_genes/list_disease_genes_PheWAS_GWAS.dat")

### 2. Align data types #############################
list_disease_genes_PheWAS_GWAS <- sapply(list_disease_genes_PheWAS_GWAS, as.integer)

### 3. Merge ##################################
list_disease_genes <- c(list_disease_EntrezIDs, list_disease_genes_PheWAS_GWAS)

save(list_disease_genes, file = "./path/to/intermediate/data/0.data_make/human_disease_genes/human_disease_genes_230714_iida.dat")

### 4. Output ##################################
df_genes <- as.data.frame(unlist(list_disease_genes))
df_genes$Disease_name <- rep(names(list_disease_genes), sapply(list_disease_genes, length))
colnames(df_genes)[1] <- "GeneID"
rownames(df_genes)    <- 1:length(df_genes$GeneID)
df_genes <- df_genes[,c(2,1)]

write.xlsx(df_genes, file = "./path/to/intermediate/data/0.data_make/human_disease_genes/human_disease_genes_230714_iida.xlsx")
