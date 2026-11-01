############################################
#                                          #      
#             make module                  #
#            causal genes                  #
#           goiter                #
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
in_fa <- "/path/to/your/local/directory/OMIM-Entry-Retrieval_H00032_goiter.xlsx"

OMIM_H00032  <- read_excel(in_fa, skip = 4)
OMIM_H00032   <- na.omit(unique(c(OMIM_H00032$`Entrez Gene ID`, OMIM_H00032$`Entrez Gene ID`)))

### 2. download CinVar genes ################
#1 access to https://www.ncbi.nlm.nih.gov/clinvar/
#2 search the disease name through search window
#3 choose Download "
#4 Format: tabular  sort by: location
#5 create file

in_fa <- "/path/to/your/local/directory/clinvar_result_H00032_goiter.xlsx"

CinVar_H00032 <- read_excel(in_fa)

### 3. download GWAS catalog genes ################
#1 access to https://www.ebi.ac.uk/gwas/
#2 search the disease name through search window
#  "hepatocellular carcinoma"のみ登録があった
#3 Associations
#4 Click the arrow button
in_fa  <- "/path/to/your/local/directory/gwas-association-downloaded_2023-07-14-EFO_0004237-withChildTraits_H00032_goiter.tsv"

GWAS_H00032    <- fread(in_fa)
GWAS_H00032 <- GWAS_H00032[as.numeric(sapply(strsplit(as.character(GWAS_H00032$`P-VALUE`), "-"), function(x){x[2]})) >= 8,]

### 4. download PheWAS catalog genes ################
#1 access to https://phewascatalog.org
#2 choose "PheWAS of GWAS Catalog of SNPs"
#3 Download this full catalog 
#4 this operation could retrieve all the file

in_f   <- "/path/to/your/local/directory/phewas-catalog_2023-03-2.csv"
PheWAS <- read_csv(in_f)
PheWAS_a   　 <- PheWAS[grep("Simple goiter",PheWAS$`phewas phenotype`),]     #"goiter"のみ登録あり
PheWAS_H00032 <- PheWAS_a[PheWAS_a$`p-value` < 0.05,]


### 7. 遺伝子名を集める#######################################
list_VG_H00032 <- list(OMIM     = OMIM_H00032,
                       ClinVar  = unique(unlist(strsplit(as.character(as.matrix(na.omit(CinVar_H00032[grep("^NM_", CinVar_H00032$Name), "Gene(s)"]))), "\\|"))),
                       GWAS     = unique(gsub(" ", "",unlist(strsplit(unlist(strsplit(GWAS_H00032$`MAPPED_GENE`, ";" )), ","))))[unique(gsub(" ", "",unlist(strsplit(unlist(strsplit(GWAS_H00032$`MAPPED_GENE`, ";" )), ",")))) != "'-"], #gsubの方がよかったかも
                       PheWAS   = unique(PheWAS_H00032$gene_name[PheWAS_H00032$gene_name != "NULL"]))


save(list_VG_H00032, file="./Analysis_230701/Causal_genes_H00032.dat")















