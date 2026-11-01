##########################################
#                                        #
#                                        #
#         Make disease module            #
#                                        #
#                                        #
##########################################
source("function.R")
library(igraph)

### 1. load ##############################
out_f <- "/path/to/your/local/directory/0.make_data"
load("./Analysis_230701/0.make_data/g_vg_cop.dat")    #human interactome
load("./Analysis_230701/0.make_data/human_disease_genes/human_disease_genes_230714_iida.dat")

### 2. PPIにある遺伝子だけに絞る ##############################
list_disease_module <- sapply(list_disease_genes, function(x){x[x %in% names(vg_cop)]})

save(list_disease_module, file = paste0(out_f,"list_disease_module.dat"))
