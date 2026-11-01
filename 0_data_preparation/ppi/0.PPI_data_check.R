############################
#                          #
#     PPI data check       #
#                          #
############################
library(igraph)

### 1. load #################
in_f <- "/path/to/your/local/directory/g_PPI_ext.dat"
load(in_f)

g_PPI_ext <- upgrade_graph(g_PPI_ext)

### 2. output #################
df <- as_edgelist(g_PPI_ext)

save(g_PPI_ext, file = "./path/to/intermediate/data/0.data_check/human_PPI_230710_iida.dat")
write.table(df, file = "./path/to/intermediate/data/0.data_check/human_PPI_230710_iida.txt", row.names = F, col.names = F, quote = F)
