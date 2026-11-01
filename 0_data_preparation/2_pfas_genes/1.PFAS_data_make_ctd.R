#########################
#                       #
#    Check the PFAS     #
#         CTD           #
#########################
# Check how many PFAS are included
library(tidyverse)
library(tictoc)
library(fs)
library(vroom)
library(data.table)
library(openxlsx)

### 0. Download source data ##############
# Run this section once to download the three required input files.
# Set `data_dir` to the folder where you want to save the files.
data_dir <- "/path/to/your/local/directory"

# --- (a) EPA PFAS master list ---
# Source: EPA CompTox Chemicals Dashboard
# https://www.epa.gov/chemical-research/per-and-polyfluoroalkyl-substances-pfas
pfas_url  <- "https://comptox.epa.gov/dashboard/chemical-lists/pfasmaster/xlsx"
pfas_file <- file.path(data_dir, "Chemical List pfasmaster-2023-01-16.csv")
# download.file(pfas_url, destfile = pfas_file)   # uncomment to download

# --- (b) CTD chemical-gene interactions ---
# Source: Comparative Toxicogenomics Database (http://ctdbase.org/)
# Register at ctdbase.org, then download CTD_chem_gene_ixns.csv from the Downloads page.
ctd_ixns_file <- file.path(data_dir, "CTD_chem_gene_ixns.csv")
# download.file("http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz", destfile = paste0(ctd_ixns_file, ".gz"))
# system(paste("gunzip", paste0(ctd_ixns_file, ".gz")))   # uncomment to download and decompress

# --- (c) CTD chemicals vocabulary ---
# Source: Comparative Toxicogenomics Database (http://ctdbase.org/)
# Download CTD_chemicals.csv from the Downloads page.
ctd_chems_file <- file.path(data_dir, "CTD_chemicals.csv")
# download.file("http://ctdbase.org/reports/CTD_chemicals.csv.gz", destfile = paste0(ctd_chems_file, ".gz"))
# system(paste("gunzip", paste0(ctd_chems_file, ".gz")))   # uncomment to download and decompress

### 1. Load ##############
# EPA PFAS master list
df_EPA_PFAS <- read.csv(pfas_file)

# CTD chemical-gene interactions (human proteins only selected in step 2)
ctd <- fread(ctd_ixns_file, skip=29, header = F) #2502394 interactions
colnames(ctd) <- c("ChemicalName", "ChemicalID", "CasRN", "GeneSymbol",
                   "GeneID", "GeneForms", "Organism", "OrganismID",
                   "Interaction", "InteractionActions", "PubMedIDs")

# CTD chemicals vocabulary
CTD_chems <- fread(ctd_chems_file, skip=29, header = F)
colnames(CTD_chems) <- c("ChemicalName","ChemicalID","CasRN","Definition","ParentIDs","TreeNumbers","ParentTreeNumbers","Synonyms")


### 2. Chem-genes data selection ######################
#Chem-genes
ctd_human <- ctd[ctd$Organism == "Homo sapiens",]        #1120037 interactions
ctd_hPPI <- ctd_human[ctd_human$GeneForms == "protein",] #364122 interactions for 9914 chemicals

save(ctd_hPPI, file = "./path/to/intermediate/data/0.data_make/ctd_hPPI.dat")

### 3. PFAS data selection ############################
# PFAS found in both the EPA list and CTD
targetID <- sort(unique(c(which(CTD_chems$ChemicalName %in% df_EPA_PFAS$PREFERRED.NAME), 
              which(CTD_chems$CasRN %in% df_EPA_PFAS$CASRN)))) #140 PFAS

CTD_PFAS_EPA    <- CTD_chems[targetID,] #140 PFAS

### 4. PFAS-genes data selection #######################
CTD_PFAS_genes <- ctd_hPPI[(ctd_hPPI$ChemicalName %in% CTD_PFAS_EPA$ChemicalName),] #18 PFAS

### 5. Check for missing entries ################################
# Fluorinated compounds registered in CTD database
CTD_PFAS_D005466 <- CTD_chems[grep("D005466",CTD_chems$ParentIDs),]
CTD_PFAS_D006845 <- CTD_chems[grep("D006845",CTD_chems$ParentIDs),]
CTD_PFAS_D005464 <- CTD_chems[grep("D005464",CTD_chems$ParentIDs),]
CTD_PFAS_D004737 <- CTD_chems[grep("D004737",CTD_chems$ChemicalID),]

# PFAS with known chem-gene interactions in CTD database
z <- unique(c(CTD_PFAS_EPA$ChemicalName, 
              #CTD_PFAS_EPAcas$ChemicalName, 
              CTD_PFAS_D005466$ChemicalName, 
              CTD_PFAS_D006845$ChemicalName,
              CTD_PFAS_D005464$ChemicalName,
              CTD_PFAS_D004737$ChemicalName))
CTD_PFAS_genes_conf <- ctd_hPPI[ctd_hPPI$ChemicalName %in% z,]
length(unique(CTD_PFAS_genes_conf$ChemicalName)) #50 PFAS

### 6. check CTD ##############################################
# Export the 50 PFAS to add CAS numbers manually
df <- as.data.frame(table(CTD_PFAS_genes_conf$ChemicalName))
write.xlsx(df, file = "./path/to/intermediate/data/0.data_make/df_CTD_EPA_PFAS_20230712.xlsx")

# 49 of the above PFAS had confirmed CAS numbers
df_conf <- read.xlsx("./path/to/intermediate/data/0.data_make/df_CTD_EPA_PFAS_20230712.xlsx")

# 41 of the above PFAS had CAS numbers confirmed in the EPA list
df_conf_EPA <- df_conf[sapply(strsplit(df_conf$CAS, "\\|"), function(x){sum(x %in% df_EPA_PFAS$CASRN)}) > 0,]
df_conf_EPA$CAS_1 <- sapply(strsplit(df_conf_EPA$CAS, "\\|"), function(x){x[1]})

  # Verify no entries are missing
sum(unique(CTD_PFAS_genes$ChemicalName) %in% df_conf_EPA$Var1) # all CAS-matched PFAS are present

### 7. Create df_CTD ############################################
df_conf_CTD <- ctd_hPPI[ctd_hPPI$ChemicalName  %in% df_conf_EPA$Var1,]
#unique(df_conf_CTD$ChemicalName) #41 PFAS
#sum(sapply(strsplit(df_conf$CAS, "\\|"), function(x){sum(x %in% df_EPA_PFAS$CASRN)}) ) #46 PFAS based on CAS

### 8. merge ############################################
# Merge with EPA PFAS data
df_EPA_PFAS_conf      <- merge(df_EPA_PFAS, df_conf_EPA[, c(1,4)], by.x = "CASRN", by.y = "CAS_1")
df_EPA_PFAS_target_ctd <- merge(df_conf_CTD, df_EPA_PFAS_conf, by.x = "ChemicalName", by.y = "Var1")

save(df_EPA_PFAS_target_ctd, file="./path/to/intermediate/data/0.data_make//df_EPA_PFAS_target_ctd.dat")
