############################
#  Supplementary Data 2    #
#       PFAS genes         #
############################
library(igraph)

### 1. load #################################
in_f <- "/Users/iida/Documents/Lab_Iida/2023_Iida/Interactome/Human_PPI/Analysis_230701/0.make_data/human_PFAS_genes/human_PFAS_genes_230712_iida.dat"
load(in_f)

### 2. output ##################
write.table(list_PFAS_genes$Merged, file = "./Supplementary_Data_2_PFASgenes.txt", 
            quote = F, row.names = F, col.names = F)

