##################################################
#                                                #
#   Human PPI network construction               #
#                                                #
#   Data source:                                 #
#     Integrated human PPI network built from    #
#     multiple databases (see DATA_SOURCES.md):  #
#       - APID (http://apid.dep.usal.es/)        #
#       - BioGRID (https://thebiogrid.org/)      #
#       - CORUM (https://mips.helmholtz-         #
#           muenchen.de/corum/)                  #
#       - HuRI (http://www.interactome-atlas.org/)
#       - InnateDB (https://www.innatedb.com/)   #
#       - IntAct / MINT / MatrixDB /             #
#         DIP / HPRD / PhosphoSitePlus           #
#     Filtered to human proteins with at least   #
#     2 experimental evidence entries.           #
#     Vertex IDs = NCBI Gene IDs (Entrez).       #
#                                                #
#   Input:  ../../data/Supplementary_Data_1_     #
#            humanPPI.txt  (edge list, no header)#
#   Outputs (intermediate):                      #
#     g_comp.dat     igraph LCC object           #
#     g_cut_degs.dat degree-bin vector           #
#     node_names.dat character vector of LCC IDs #
#                                                #
##################################################
library(igraph)
source("../../function.R")

### 1. Load edge list ########################################
in_f  <- "../../data/Supplementary_Data_1_humanPPI.txt"
out_f <- "./path/to/intermediate/data/"

df_edges <- read.table(in_f, header = FALSE,
                       col.names = c("GeneID1", "GeneID2"),
                       colClasses = "character")

### 2. Build undirected graph ################################
g_PPI <- graph_from_data_frame(df_edges, directed = FALSE)
g_PPI <- simplify(g_PPI)   # remove self-loops and duplicate edges

cat("PPI network (full):", vcount(g_PPI), "nodes,", ecount(g_PPI), "edges\n")

### 3. Extract largest connected component (LCC) #############
cls    <- clusters(g_PPI, "weak")
g_comp <- induced_subgraph(g_PPI,
           V(g_PPI)[cls$membership == which.max(cls$csize)])

cat("PPI network (LCC): ", vcount(g_comp), "nodes,", ecount(g_comp), "edges\n")

### 4. Degree bins (used by proximity z-score computation) ###
deg      <- degree(g_comp)
cut_degs <- Binning(deg, nums = 100)
names(cut_degs) <- names(deg)

### 5. Save ##################################################
node_names <- names(V(g_comp))   # character vector of LCC gene IDs

save(g_comp,     file = paste0(out_f, "g_comp.dat"))
save(cut_degs,   file = paste0(out_f, "g_cut_degs.dat"))
save(node_names, file = paste0(out_f, "node_names.dat"))
