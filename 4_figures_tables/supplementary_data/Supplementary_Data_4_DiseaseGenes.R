############################
#  Supplementary Data 4    #
#      disease genes       #
############################
library(igraph)

### 1. load ####################
# list_disease_genes produced by 0_data_preparation/4.disease_data_load.R
in_f <- "./path/to/intermediate/data/list_disease_genes.dat"
load(in_f)

### 5. save ##########################
df_disease <- data.frame(unlist(list_disease_genes))
cleaned <- sub("[0-9]+$", "", rownames(df_disease))
df_disease <- as.data.frame(cbind("DiseaseNames"=cleaned, "EntrezGeneID" = df_disease$unlist.list_disease_genes.))

write.table(df_disease, file = "./Supplementary_Data_4_DiseaseGenes.txt", 
            quote = F, row.names = F, col.names = T)

