############################
#  Supplementary Data 6    #
# PFAS-disease distances   #
############################
library(readxl)

#### 1. load ###############################
in_f <- "/path/to/your/local/directory/summary_PFAS_diseaes_iidaData_20230716.xlsx"
data <- read_excel(in_f)
colnames(data)[c(5,6)] <- c("sd(dc_null)", "Zdc(dc_null)")

write.table(data, file = "./path/to/intermediate/data/Supplementary_Data_6_Zdc_PFAS_disease.txt",
            quote = F, row.names = F, col.names = T)

#### 2. for manuscript ######################
PFOA <- data[grep("335-67-1", data$ID),]
PFOA <- PFOA[order(PFOA$`Zdc(dc_null)`),]
PFOA[PFOA$`Zdc(dc_null)` > 2, ]

PFOS <- data[grep("1763-23-1", data$ID),]
PFOS <- PFOS[order(PFOS$`Zdc(dc_null)`),]
PFOS[PFOS$`Zdc(dc_null)` > 1, ]
