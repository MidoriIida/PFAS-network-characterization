############################
#                          #
#     PPI data check       #
#                          #
############################
library(igraph)

### 1. load #################
in_f <- "/Users/iida/Documents/Lab_Yamanishi/2022_Adverse/Analysis_220517/g_PPI_ext.dat"
load(in_f)

g_PPI_ext <- upgrade_graph(g_PPI_ext)

### 2. output #################
df <- as_edgelist(g_PPI_ext)

save(g_PPI_ext, file = "./Analysis_230701/0.data_check/human_PPI_230710_iida.dat")
write.table(df, file = "./Analysis_230701/0.data_check/human_PPI_230710_iida.txt", row.names = F, col.names = F, quote = F)
