############################
#                          #
#      PPI data make       #
#                          #
############################
library(igraph)
source("function.R")

### 1. load #################
out_f <- "./path/to/intermediate/data/"
in_f <- "../../data/g_PPI_ext.dat"   # pre-built PPI network (see DATA_SOURCES.md)
load(in_f)

g_PPI_ext <- upgrade_graph(g_PPI_ext)

### 2. Extract largest connected component ##########################
cls <- clusters(g_PPI_ext, "weak")
g_comp <- delete.vertices(g_PPI_ext, subset(V(g_PPI_ext),
                          cls$membership != which(cls$csize == max(cls$csize))[[1]]))

cat("PPI network (full):", vcount(g_PPI_ext), "nodes,", ecount(g_PPI_ext), "edges\n")
cat("PPI network (LCC): ", vcount(g_comp),    "nodes,", ecount(g_comp),    "edges\n")

### 3. Get node names in the network #######
vg_cop <- V(g_comp)
names_vg_comp <- names(V(g_comp))

### 4. Assign degree bins ########################
deg      <- degree(g_comp)
cut_degs <- Binning(deg, nums = 100)
names(cut_degs) <- names(deg)

### 5. Save ##########################################
#save(vg_cop,   file = paste0(out_f, "g_vg_cop.dat")) # saving vg_cop changes the data structure
save(cut_degs, file = paste0(out_f, "g_cut_degs.dat"))
save(g_comp,   file = paste0(out_f, "g_comp.dat"))

