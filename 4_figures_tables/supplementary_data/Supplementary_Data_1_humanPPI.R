############################
#  Supplementary Data 1    #
#        Human PPI         #
############################
library(igraph)

### 1. load ####################
# g_comp produced by 0_data_preparation/1.PPI_data_load.R
load("./path/to/intermediate/data/g_comp.dat")

### 2. output ##################
df_edge <- as_edgelist(g_comp, names = TRUE)
write.table(df_edge, file = "../../data/Supplementary_Data_1_humanPPI.txt",
            quote = F, row.names = F, col.names = F)
