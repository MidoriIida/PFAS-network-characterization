############################
#                          #
#      PPI data make       #
#                          #
############################
library(igraph)
source("function.R")

### 1. load #################
out_f <- "/path/to/your/local/directory/0.make_data"
in_f <- "/path/to/your/local/directory/g_PPI_ext.dat"
load(in_f)

g_PPI_ext <- upgrade_graph(g_PPI_ext)

### 2.最大成分を抽出 ########################## 
cls <- clusters(g_PPI_ext,"weak")
g_comp <- delete.vertices(g_PPI_ext,subset(V(g_PPI_ext),
                                   + cls$membership!=which(cls$csize==max(cls$csize))[[1]]))

### 3. ネットワークのノードの名前を指定 #######
vg_cop <- V(g_comp)
names_vg_comp <- names(V(g_comp))

### 4. degree binを指定 ########################
deg      <- degree(g_comp)
cut_degs <- Binning(deg, nums = 100)
names(cut_degs) <- names(deg)

### 5. 保存 ##########################################
#save(vg_cop,   file = paste0(out_f, "g_vg_cop.dat")) #保存するとデータが変わってしまう
save(cut_degs, file = paste0(out_f, "g_cut_degs.dat"))
save(g_comp,   file = paste0(out_f, "g_comp.dat"))

