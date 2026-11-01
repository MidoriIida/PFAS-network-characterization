############################
#  Supplementary Data 5    #
#   PFAS-drug distances    #
############################

#### 1. load ###############################
# list_dist_mat_drug_PFAS produced by 1_drug_space/1.compute_network_proximity_between_PFAS_drug.R
load("./path/to/intermediate/data/list_dist_mat_drug_PFAS.dat")

# flatten list to data frame
data <- do.call(rbind, list_dist_mat)
colnames(data) <- c("ID", "Name", "dc_real", "sd(dc_null)", "mean(dc_null)", "Zdc(dc_null)")

write.table(data, file = "../../data/Supplementary_Data_5_Zdc_PFAS_Drug.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")

#### 2. for manuscript ######################
data$`Zdc(dc_null)` <- as.numeric(data$`Zdc(dc_null)`)

PFOA <- data[grep("335-67-1", data$ID), ]
PFOA <- PFOA[order(PFOA$`Zdc(dc_null)`), ]
PFOA[PFOA$`Zdc(dc_null)` > 2, ]

PFOS <- data[grep("1763-23-1", data$ID), ]
PFOS <- PFOS[order(PFOS$`Zdc(dc_null)`), ]
PFOS[PFOS$`Zdc(dc_null)` > 2, ]
