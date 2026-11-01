#################################
#          Table 3              #
#   Drug ATC for legacy PFAS    #
#################################
library(readxl)
library(openxlsx)

### 1. load data  ####################
in_d2  <- "/path/to/your/local/directory/orignail_data"
in_f2 <- "/20230830/result_enrichment_ATC_code_KEGGdrugBackground_20230825.csv"

ATC        <- read.csv(paste0(in_d2, in_f2), row.names = 1)

### 2. extract ####################
legacyPFAS <- c("375-22-4", "375-85-9", "335-67-1", "375-95-1", "335-76-2", "2058-94-8", "307-55-1",
                "375-73-5", "355-46-4", "375-92-8", "1763-23-1", "754-91-6")
AbbPFAS <- c(
  "PFBA",  #Perfluorobutyric acid C4HF7O2 4
  "PFHpA", #Perfluoroheptanoic acid C6F13COOH    6
  "PFOA",  #Perfluorooctanoic acid C8F15O2X 8
  "PFNA" , #Perfluorononanoic acid C9HF17O2 9
  "PFDA",  #perfluoro-n-decanoic acid C9F19COOH 10
  "PFUdA" ,#Perfluoroundecanoic acid C11HF21O2 11
  "PFDoA", #Perfluorododecanoic acid C12HF23O2 12
  "PFBS" ,#Perfluorobutanesulfonic acid C4F9SO3H 4
  "PFHxS" ,#Perfluorohexanesulfonic Acid C6HF13O3S 16
  "PFHpS" ,# Perfluoroheptanesulfonic acid C7HF15O3S 7
  "PFOS" ,#perfluorooctanesulfonic acid C8HF17O3S 8
  "PFOSA" #Perfluorooctanesulfonamide C8H2F17NO2S 8
)

data <- ATC[,sapply(gsub("-", ".", legacyPFAS), function(x){grep(x, colnames(ATC))})]
colnames(data) <- AbbPFAS

res <- round(data[order(apply(data, 1, function(x){mean(na.omit(x))})),] , digits = 3)

### 3. output ####################
write.xlsx(res, file = "./TableS1_Drug_ATC_Enrichment_Analysis.xlsx", rowNames = T)
