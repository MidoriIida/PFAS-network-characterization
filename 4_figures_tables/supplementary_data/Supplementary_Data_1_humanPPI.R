############################
#  Supplementary Data 1    #
#        Human PPI         #
############################
library(igraph)

### 1. load ####################
in_f <- "/Users/iida/Documents/Lab_Iida/2023_Iida/Interactome/Human_PPI/Analysis_230701/0.make_data/human_PPI/human_PPI_230710_iida.dat"
load(in_f)

### 2. output ##################
df_edge <- as_edgelist(g_PPI_ext, names = TRUE)
write.table(df_edge, file = "./Supplementary_Data_1_humanPPI.txt", 
            quote = F, row.names = F, col.names = F)
