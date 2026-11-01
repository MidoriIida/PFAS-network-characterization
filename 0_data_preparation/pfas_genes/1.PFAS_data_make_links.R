#########################
#                       #
#    Check the PFAS     #
#        LINCS          #
#########################
#PFASいくつ入ってる？
library(tidyverse)
library(tictoc)
library(fs)
library(vroom)
library(data.table)

### 1. load ##############
#EPAのPFAS list
in_f <- "/path/to/your/local/directory/Chemical List pfasmaster-2023-01-16/Chemical List pfasmaster-2023-01-16.csv"
df_EPA_PFAS <- read.csv(in_f)

#Links drugs
in_f2 <- "/path/to/your/local/directory/Original data/LINCS2/L1000_level5_signatures_978genes.tabs"
df_LINKS <- fread(in_f2)

#drug-protein
in_f <- "/path/to/your/local/directory/CTD_chem_gene_ixns.csv"
ctd <- read.csv(paste0(in_f), skip=29, header = F) #2502394 interactions
colnames(ctd) <- c("ChemicalName",	"ChemicalID",	"CasRN",	"GeneSymbol",
                   "GeneID", "GeneForms",	"Organism",	"OrganismID",
                   "Interaction",	"InteractionActions",	"PubMedIDs")



### 2. check the LINCS ##############
LINCS_inchi <- unique(df_LINKS$inchi_key)[nchar(unique(df_LINKS$inchi_key)) == 27]
target <- df_EPA_PFAS$INCHIKEY[df_EPA_PFAS$INCHIKEY %in% LINCS_inchi]

length(target) #8

### 3. dataを取得 ##############################################################
#EPA data
df_EPA_PFAS_target <- df_EPA_PFAS[df_EPA_PFAS$INCHIKEY %in% LINCS_inchi,]

#LINCS data
df_LINCS_target <- df_LINKS[df_LINKS$inchi_key %in% target,]
df_LINCS_target_exp <- df_LINCS_target[,10:dim(df_LINCS_target)[2]]

df_LINCS_target_exp_mean <- apply(df_LINCS_target_exp, 2, function(x){tapply(x, df_LINCS_target$inchi_key, mean)})


ceiling(dim(df_LINCS_target_exp_mean)[2]*0.05)

PFAS_LINCS_Proteins <- NULL
for(i in 1:dim(df_LINCS_target_exp_mean)[1]){
  x <- sort(df_LINCS_target_exp_mean[i,],decreasing = T)[ceiling(dim(df_LINCS_target_exp_mean)[2]*0.05)]
  y <- sort(df_LINCS_target_exp_mean[i,],decreasing = F)[ceiling(dim(df_LINCS_target_exp_mean)[2]*0.05)]
  
  genes <- colnames(df_LINCS_target_exp_mean)
  genes_x <- genes[df_LINCS_target_exp_mean[i,] >= as.numeric(x)]
  genes_y <- genes[df_LINCS_target_exp_mean[i,] <= as.numeric(y)]
  
  genes_xy <- gsub("hsa.", "",c(genes_x, genes_y))
  
  GeneID <- as.data.frame(cbind(rownames(df_LINCS_target_exp_mean)[i], genes_xy))
  PFAS_LINCS_Proteins <- rbind(PFAS_LINCS_Proteins, GeneID)
}

### 4. dataをmerge ##############################################################
df_EPA_PFAS_target_Lincs <- merge(df_EPA_PFAS_target, PFAS_LINCS_Proteins, by.x = "INCHIKEY", by.y = "V1")

save(df_EPA_PFAS_target_Lincs, file="./Analysis_230701/0.data_check/df_EPA_PFAS_target_Lincs.dat")

