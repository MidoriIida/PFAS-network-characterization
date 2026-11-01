##########################################
#                                        #
#                                        #
#         Make chemical module            #
#                                        #
#                                        #
##########################################
source("function.R")
library(igraph)

### 1. Load ##############################
out_f <- "./path/to/intermediate/data/"
load("./path/to/intermediate/data/node_names.dat")       # LCC gene IDs from 1_ppi
load("./path/to/intermediate/data/list_PFAS_genes.dat")  # from 2_pfas_genes

### 2. Convert to per-PFAS gene list ##############################
list_chems <- tapply(list_PFAS_genes$Merged$GeneID, list_PFAS_genes$Merged$CASRN, list)

### 3. Filter to genes present in the PPI network ##############################
list_chem_module <- sapply(list_chems, function(x){x[x %in% node_names]})

save(list_chem_module, file = paste0(out_f, "list_chemical_module.dat"))


