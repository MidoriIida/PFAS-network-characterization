#########################
#                       #
#    Check the PFAS     #
#         CTD           #
#########################
#PFASいくつ入ってる？
library(tidyverse)
library(tictoc)
library(fs)
library(vroom)
library(data.table)
library(openxlsx)

### 1. load ##############
#EPAのPFAS list
in_f <- "/path/to/your/local/directory/Chemical List pfasmaster-2023-01-16/Chemical List pfasmaster-2023-01-16.csv"
df_EPA_PFAS <- read.csv(in_f)

#drug-protein
in_f <- "/path/to/your/local/directory/CTD_chem_gene_ixns.csv"
ctd <- fread(paste0(in_f), skip=29, header = F) #2502394 interactions
colnames(ctd) <- c("ChemicalName",	"ChemicalID",	"CasRN",	"GeneSymbol",
                   "GeneID", "GeneForms",	"Organism",	"OrganismID",
                   "Interaction",	"InteractionActions",	"PubMedIDs")

#CTD_chems
in_f2 <- "/path/to/your/local/directory/CTD_chemicals.csv"
CTD_chems <- fread(in_f2, skip=29, header = F)
colnames(CTD_chems) <- c("ChemicalName","ChemicalID","CasRN","Definition","ParentIDs","TreeNumbers","ParentTreeNumbers","Synonyms")

#Takemoto lab PFAS
load("./Analysis_230701/test.dat")


### 2. Chem-genes data selection ######################
#Chem-genes
ctd_human <- ctd[ctd$Organism == "Homo sapiens",]        #1120037 interactions
ctd_hPPI <- ctd_human[ctd_human$GeneForms == "protein",] #364122 interactions for 9914 chemicals

save(ctd_hPPI, file = "./Analysis_230701/0.data_make/ctd_hPPI.dat")

### 3. PFAS data selection ############################
#EPAに名前のあったPFASのうちCTDにも名前があったPFAS
targetID <- sort(unique(c(which(CTD_chems$ChemicalName %in% df_EPA_PFAS$PREFERRED.NAME), 
              which(CTD_chems$CasRN %in% df_EPA_PFAS$CASRN)))) #140 PFAS

CTD_PFAS_EPA    <- CTD_chems[targetID,] #140 PFAS

### 4. PFAS-genes data selection #######################
CTD_PFAS_genes <- ctd_hPPI[(ctd_hPPI$ChemicalName %in% CTD_PFAS_EPA$ChemicalName),] #18 PFAS

### 5. 漏れがないか確認 ################################
#CTDデータベースでフッ素化合物とされているもの
CTD_PFAS_D005466 <- CTD_chems[grep("D005466",CTD_chems$ParentIDs),]
CTD_PFAS_D006845 <- CTD_chems[grep("D006845",CTD_chems$ParentIDs),]
CTD_PFAS_D005464 <- CTD_chems[grep("D005464",CTD_chems$ParentIDs),]
CTD_PFAS_D004737 <- CTD_chems[grep("D004737",CTD_chems$ChemicalID),]

#CTD databaseでChem-geneがわかっているPFAS
z <- unique(c(CTD_PFAS_EPA$ChemicalName, 
              #CTD_PFAS_EPAcas$ChemicalName, 
              CTD_PFAS_D005466$ChemicalName, 
              CTD_PFAS_D006845$ChemicalName,
              CTD_PFAS_D005464$ChemicalName,
              CTD_PFAS_D004737$ChemicalName))
CTD_PFAS_genes_conf <- ctd_hPPI[ctd_hPPI$ChemicalName %in% z,]
length(unique(CTD_PFAS_genes_conf$ChemicalName)) #50 PFAS

### 6. check CTD ##############################################
#上記50個のPFASを一旦外に出して、CAS numberを入れる
df <- as.data.frame(table(CTD_PFAS_genes_conf$ChemicalName))
write.xlsx(df, file = "./Analysis_230701/0.data_make/df_CTD_EPA_PFAS_20230712.xlsx")

#上記PFASのうちCASが確認できたPFASは49個だった
df_conf <- read.xlsx("./Analysis_230701/0.data_make/df_CTD_EPA_PFAS_20230712.xlsx")

#上記PFASのうちCASがEPAリストに確認できたPFASは41個だった
df_conf_EPA <- df_conf[sapply(strsplit(df_conf$CAS, "\\|"), function(x){sum(x %in% df_EPA_PFAS$CASRN)}) > 0,]
df_conf_EPA$CAS_1 <- sapply(strsplit(df_conf_EPA$CAS, "\\|"), function(x){x[1]})

  #漏れがないか確認
sum(unique(CTD_PFAS_genes$ChemicalName) %in% df_conf_EPA$Var1) #CASでとれてきたやつは全てある

### 7. df_CTDを作成 ############################################
df_conf_CTD <- ctd_hPPI[ctd_hPPI$ChemicalName  %in% df_conf_EPA$Var1,]
#unique(df_conf_CTD$ChemicalName) #41 PFAS
#sum(sapply(strsplit(df_conf$CAS, "\\|"), function(x){sum(x %in% df_EPA_PFAS$CASRN)}) ) #46 PFAS based on CAS

### 8. merge ############################################
#EPA dataを取得
df_EPA_PFAS_conf      <- merge(df_EPA_PFAS, df_conf_EPA[, c(1,4)], by.x = "CASRN", by.y = "CAS_1")
df_EPA_PFAS_target_ctd <- merge(df_conf_CTD, df_EPA_PFAS_conf, by.x = "ChemicalName", by.y = "Var1")

save(df_EPA_PFAS_target_ctd, file="./Analysis_230701/0.data_make//df_EPA_PFAS_target_ctd.dat")
