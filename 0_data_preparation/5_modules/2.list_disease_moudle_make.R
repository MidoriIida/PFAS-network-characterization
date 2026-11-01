##########################################
#                                        #
#                                        #
#         Make disease module            #
#                                        #
#                                        #
##########################################
source("function.R")
library(igraph)

### 1. Load ##############################
out_f <- "./path/to/intermediate/data/"
load("./path/to/intermediate/data/node_names.dat")            # LCC gene IDs from 1_ppi
load("./path/to/intermediate/data/list_disease_genes.dat")    # from 3_disease_genes

### 2. Filter to genes present in the PPI network ##############################
list_disease_module <- sapply(list_disease_genes, function(x){x[x %in% node_names]})

save(list_disease_module, file = paste0(out_f, "list_disease_module.dat"))
