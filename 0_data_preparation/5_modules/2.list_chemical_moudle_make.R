##########################################
#                                        #
#                                        #
#         Make chemical module            #
#                                        #
#                                        #
##########################################
source("function.R")
library(igraph)

### 1. load ##############################
out_f <- "/path/to/your/local/directory/0.make_data"
load("./path/to/intermediate/data/0.make_data/g_vg_cop.dat")    #human interactome
load("./path/to/intermediate/data/0.make_data/human_PFAS_genes/human_PFAS_genes_230712_iida.dat")

### 2. Convert back to list_chem format ##############################
list_chems <- tapply(list_PFAS_genes$Merged$GeneID, list_PFAS_genes$Merged$CASRN, list)


### 2. Filter to genes present in the PPI network ##############################
list_chem_module <- sapply(list_chems, function(x){x[x %in% names(vg_cop)]})

save(list_chem_module, file = paste0(out_f,"list_chemical_module.dat"))


