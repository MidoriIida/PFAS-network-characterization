############################
#  Supplementary Data 2    #
#       PFAS genes         #
############################

### 1. Load ################################
# list_PFAS_genes produced by 0_data_preparation/2.PFAS_genes_load.R
in_f <- "./path/to/intermediate/data/list_PFAS_genes.dat"
load(in_f)

### 2. Output ##################
write.table(list_PFAS_genes$Merged, file = "../../data/Supplementary_Data_2_PFASgenes.txt",
            quote = F, row.names = F, col.names = F)

