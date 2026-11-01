############################
#  Supplementary Data 1    #
#        Human PPI         #
############################
library(igraph)

### 1. load ####################
in_f <- "/path/to/your/local/directory/human_PPI_230710_iida.dat"
load(in_f)

### 2. output ##################
df_edge <- as_edgelist(g_PPI_ext, names = TRUE)
write.table(df_edge, file = "./Supplementary_Data_1_humanPPI.txt", 
            quote = F, row.names = F, col.names = F)
