############################################
#                                          #      
#             make module                  #
#            causal genes                  #
#           Hyperlipidemia                #
############################################
#cheng et al 2019 used the 8 databases below
#to get genes which associated to disease
#1 OMIM
#2 (CTD)
#3 ClinVar
#4 GWAS Catalog
#5 (GWASdb)
#6 PheWAS Catalog
#7 (HuGE Navigator)
#8 (DisGeNET (PW diseaseNET))

#within the 8 database, 
#CTD wasn't the database for causal genes
#GWASdb & DisGeNET were out of service 
#HuGe Navigate is a Co-occurrence data
#Therefor I used 1, 3, 4, 6 databases below.

library(openxlsx)
library(readr)
library("readxl")
library(tidyverse)
library(tictoc)
library(fs)
library(vroom)
library(data.table)
library(stringi)
library(dplyr)

### 1. download OMIM genes ################
#1 access to https://www.omim.org/about
#2 go to https://www.omim.org/search/advanced/entry
# "liver cancer"
#3 choose MIM Number Prefix:* and +
#4 search
#5 Download as excel file
in_fa <- "/path/to/your/local/directory/OMIM-Entry-Retrieval_H01635_Hyperlipidemia.xlsx"

OMIM_H01635  <- read_excel(in_fa, skip = 4)
OMIM_H01635   <- na.omit(unique(c(OMIM_H01635$`Entrez Gene ID`, OMIM_H01635$`Entrez Gene ID`)))

### 2. download CinVar genes ################
#1 access to https://www.ncbi.nlm.nih.gov/clinvar/
#2 search the disease name through search window
#3 choose Download "
#4 Format: tabular  sort by: location
#5 create file

in_fa <- "/path/to/your/local/directory/clinvar_result_H01635_Hyperlipidemia.xlsx"

CinVar_H01635 <- read_excel(in_fa)

### 3. download GWAS catalog genes ################
#1 access to https://www.ebi.ac.uk/gwas/
#2 search the disease name through search window
# only "Hyperlipidemia" is registered
#3 Associations
#4 Click the arrow button
in_fa  <- "/path/to/your/local/directory/gwas-association-downloaded_2023-07-14-EFO_0004911-withChildTraits_H01635_Hyperlipidemia.tsv"

GWAS_H01635    <- fread(in_fa)
GWAS_H01635 <- GWAS_H01635[as.numeric(sapply(strsplit(as.character(GWAS_H01635$`P-VALUE`), "-"), function(x){x[2]})) >= 8,]

### 4. download PheWAS catalog genes ################
#1 access to https://phewascatalog.org
#2 choose "PheWAS of GWAS Catalog of SNPs"
#3 Download this full catalog 
#4 this operation could retrieve all the file

in_f   <- "/path/to/your/local/directory/phewas-catalog_2023-03-2.csv"
PheWAS <- read_csv(in_f)
PheWAS_a   　 <- PheWAS[grep("Hyperlipidemia",PheWAS$`phewas phenotype`),]     # only "Hyperlipidemia" is registered
PheWAS_H01635 <- PheWAS_a[PheWAS_a$`p-value` < 0.05,]


### 7. Collect gene names #######################################
list_VG_H01635 <- list(OMIM     = OMIM_H01635,
                       ClinVar  = unique(unlist(strsplit(as.character(as.matrix(na.omit(CinVar_H01635[grep("^NM_", CinVar_H01635$Name), "Gene(s)"]))), "\\|"))),
                       GWAS     = unique(gsub(" ", "",unlist(strsplit(unlist(strsplit(GWAS_H01635$`MAPPED_GENE`, ";" )), ","))))[unique(gsub(" ", "",unlist(strsplit(unlist(strsplit(GWAS_H01635$`MAPPED_GENE`, ";" )), ",")))) != "'-"], # gsub might have been a better approach
                       PheWAS   = unique(PheWAS_H01635$gene_name[PheWAS_H01635$gene_name != "NULL"]))


save(list_VG_H01635, file="./path/to/intermediate/data/Causal_genes_H01635.dat")















