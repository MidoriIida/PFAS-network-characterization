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
out_f <- "/Users/iida/Documents/Lab_Iida/2023_Iida/Interactome/Human_PPI/Analysis_230701/0.make_data/"
load("./Analysis_230701/0.make_data/g_vg_cop.dat")    #human interactome
load("./Analysis_230701/0.make_data/human_PFAS_genes/human_PFAS_genes_230712_iida.dat")

### 2. list_chemに戻す ##############################
list_chems <- tapply(list_PFAS_genes$Merged$GeneID, list_PFAS_genes$Merged$CASRN, list)


### 2. PPIにある遺伝子だけに絞る ##############################
list_chem_module <- sapply(list_chems, function(x){x[x %in% names(vg_cop)]})

save(list_chem_module, file = paste0(out_f,"list_chemical_module.dat"))


