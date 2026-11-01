############################
#  Supplementary Data 4    #
#      disease genes       #
############################
library(igraph)

### 1. load ####################
in_f <- "/Users/iida/Documents/Lab_Iida/2023_Iida/Interactome/Human_PPI/Analysis_230701/0.make_data/human_disease_genes/human_disease_genes_230714_iida.dat"
load(in_f)

### 5. save ##########################
df_disease <- data.frame(unlist(list_disease_genes))
cleaned <- sub("[0-9]+$", "", rownames(df_disease))
df_disease <- as.data.frame(cbind("DiseaseNames"=cleaned, "EntrezGeneID" = df_disease$unlist.list_disease_genes.))

write.table(df_disease, file = "./Supplementary_Data_4_DiseaseGenes.txt", 
            quote = F, row.names = F, col.names = T)

