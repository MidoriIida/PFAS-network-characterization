############################
#  Supplementary Data 6    #
# PFAS-disease distances   #
############################

#### 1. load ###############################
# list_dist_mat produced by 2_disease_space/1.compute_network_proximity_between_PFAS_disease.R
load("./path/to/intermediate/data/list_dist_mat_disease_PFAS.dat")

# flatten list to data frame
data <- do.call(rbind, list_dist_mat)
colnames(data) <- c("ID", "Name", "dc_real", "mean(dc_null)", "sd(dc_null)", "Zdc(dc_null)")
data$dc_real        <- as.numeric(data$dc_real)
data$`Zdc(dc_null)` <- as.numeric(data$`Zdc(dc_null)`)

write.table(data, file = "./path/to/intermediate/data/Supplementary_Data_6_Zdc_PFAS_disease.txt",
            quote = F, row.names = F, col.names = T)

#### 2. for manuscript ######################
PFOA <- data[grep("335-67-1", data$ID),]
PFOA <- PFOA[order(PFOA$`Zdc(dc_null)`),]
PFOA[PFOA$`Zdc(dc_null)` > 2, ]

PFOS <- data[grep("1763-23-1", data$ID),]
PFOS <- PFOS[order(PFOS$`Zdc(dc_null)`),]
PFOS[PFOS$`Zdc(dc_null)` > 1, ]
