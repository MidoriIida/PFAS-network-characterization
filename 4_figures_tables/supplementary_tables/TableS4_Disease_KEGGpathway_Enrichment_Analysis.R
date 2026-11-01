######################################
#          Table 6                   #
# Enrichment analysis of KEGGpathway #
######################################
library(readxl)
library(openxlsx)

### 1. load #########################
#"/path/to/your/local/directory/df_fdr_enrichment_KEGG_pathway.dat"
load("/path/to/your/local/directory/df_results_enrichment_KEGG_pathway.dat")
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

### 2. data make  ###########################
data <- t(df_results[legacyPFAS, ])
colnames(data) <- AbbPFAS

res <- round(data[order(apply(data, 1, function(x){mean(na.omit(x))})),], digits = 3)

### 3. output ####################
write.xlsx(res, file = "./TableS4_Disease_KEGGpathway_Enrichment_Analysis.xlsx", rowNames = T)
