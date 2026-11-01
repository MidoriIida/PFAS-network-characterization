####################################
#                                  #
#        Disease genes を          #
#        まとめてとりたい          #
#                                  #
####################################
library(tidyverse)
library(tictoc)
library(fs)
library(vroom)
library(data.table)
library(stringi)
library(dplyr)

### 1. OMIM database genes ################
#2018年のOMIM...
in_fomim  <- "/path/to/your/local/directory/morbidmap.xlsx"
in_fomim1 <- "/path/to/your/local/directory/mim2gene.txt"

OMIM      <- read_excel(in_fomim)
OMIM_gene <- fread(in_fomim1)

#ID 変換
OMIM$OMIM_GeneID  <- sapply(OMIM$`Gene Symbols`, function(genesyms){
  tem <- unique(unlist(sapply(strsplit(genesyms, ",")[[1]],  
                              function(x){OMIM_gene[(OMIM_gene$`Approved Gene Symbol (HGNC)`== x), "Entrez Gene ID (NCBI)"]})))
  if(length(tem) == 0){
    tem <- ""
  }
  return(tem)
})

#ID 変換2
result  <- stri_reverse(OMIM$Phenotype)
result2 <- sapply(strsplit(result, ","), function(x){x[1]})
result3 <- stri_reverse(sapply(strsplit(result2, " "), function(x){x[2]}))

r2 <- grep("[0-9]{6}",result2)
#r3 <- (1:length(result2))[!(1:length(result2) %in% r2)]

OMIM$MIMID <- ""
OMIM$MIMID[r2] <- result3[r2]

#必要なものだけ取ってくる
df_OMIM <- as.data.frame(t(sapply(split(OMIM[,c("MIMID", "OMIM_GeneID")], OMIM$`MIMID`), 
                                  function(x){tmp <- unique(x)
                                  tmp2 <- cbind(unique(tmp$MIMID),paste0(x$OMIM_GeneID, collapse = ","))
                                  return(tmp2)}
))[-1,])
colnames(df_OMIM) <- c("DiseaseMIM", "OMIM_GeneID")

### 2. PheWAS catalog genes ################
#1 access to https://phewascatalog.org
#2 choose "PheWAS of GWAS Catalog of SNPs"

in_f_phewas   <- "/path/to/your/local/directory/phewas-catalog_2023-03-2.csv"
PheWAS <- read_csv(in_f_phewas)

#ID 変換
PheWAS$PheWAS_GeneID  <- sapply(PheWAS$`gene_name`, function(genesyms){
  tem <- unique(unlist(sapply(strsplit(genesyms, ",")[[1]],  
                              function(x){OMIM_gene[(OMIM_gene$`Approved Gene Symbol (HGNC)`== x), "Entrez Gene ID (NCBI)"]})))
  if(length(tem) == 0){
    tem <- ""
  }
  return(tem)
})

#必要なものだけ取ってくる
list_PheWAS <- split(PheWAS[, c("phewas phenotype", "gene_name", "PheWAS_GeneID")], PheWAS$`phewas phenotype`) #1354 diseases

df_PheWAS <- as.data.frame(t(sapply(list_PheWAS , function(x){
  cbind(unique(x$"phewas phenotype"),
        paste0(x$gene_name, collapse = ","),
        paste0(x$PheWAS_GeneID, collapse = ",")
  )
})))

colnames(df_PheWAS) <- c("phewas phenotype", "gene_name", "PheWAS_GeneID")

### 3. GWAS catalog genes ################
#The NHGRI-EBI Catalog of human genome-wide association studies
#1 access to https://www.ebi.ac.uk/gwas/
#2 get "association all"

in_f_gwas <- "/path/to/your/local/directory/gwas_catalog_v1.0-associations_e109_r2023-07-05.tsv"
in_f_102  <- "/path/to/your/local/directory/gwas_catalog_v1.0.2-associations_e109_r2023-07-05.tsv"
GWAS     <- fread(in_f_gwas, header = T, sep = "\t", quote = "")
GWAS_102 <- fread(in_f_102, header = T, sep = "\t", quote = "")

#ID 変換
GWAS_102$GWAS_GeneID  <- sapply(GWAS$`MAPPED_GENE`, function(genesyms){
  tem <- unique(unlist(sapply(strsplit(genesyms, ",")[[1]],  
                              function(x){OMIM_gene[(OMIM_gene$`Approved Gene Symbol (HGNC)`== x), "Entrez Gene ID (NCBI)"]})))
  if(length(tem) == 0){
    tem <- ""
  }
  return(tem)
})

#必要なものだけ取ってくる
list_GWAS_102 <- split(GWAS_102[, c("MAPPED_TRAIT", "MAPPED_TRAIT_URI", "GWAS_GeneID")], GWAS_102$`MAPPED_TRAIT`) #8803 diseases

df_GWAS_102 <- as.data.frame(t(sapply(list_GWAS_102 , function(x){
  cbind(unique(x$MAPPED_TRAIT),
        paste0(unique(x$MAPPED_TRAIT_URI), collapse = ","),
        paste0(x$GWAS_GeneID, collapse = ",")
  )
}))[-1,])

colnames(df_GWAS_102) <- c("MAPPED_TRAIT", "MAPPED_TRAIT_URI", "GWAS_GeneID")

### 3. Clinver database genes ################
#ClinVar aggregates information about genomic variation and its relationship to human health.
#1 access tohttps://www.ncbi.nlm.nih.gov/clinvar/
#2 get "disease_names" "gene_condition_source_id"

in_f_clinver <- "/path/to/your/local/directory/gene_condition_source_id_20230714.txt"
Clinver <- fread(in_f_clinver, header = T, sep = "\t", quote = "")

#必要なところだけ抽出
list_Clinver <- split(Clinver[,c(1, 2, 5, 8)], Clinver$DiseaseName)
df_Clinver <- as.data.frame(t(sapply(list_Clinver , function(x){
  cbind(unique(x$DiseaseName),
        unique(x$DiseaseMIM),
        paste0(x$AssociatedGenes, collapse = ","),
        paste0(x$"#GeneID", collapse = ",")
  )
})))

colnames(df_Clinver) <- c("DiseaseName", "DiseaseMIM", "AssociatedGenes", "Clinver_GeneID")



### 5. dataを合体 ##############################
#ClinverとOMIMデータはOMIMで合体
Clinver_OMIM <- merge(df_Clinver, df_OMIM, by.x = "DiseaseMIM", by.y = "DiseaseMIM")

#PheWASとGWASデータは名前で合体
df_PheWAS$UP_DiseaseName <- toupper(df_PheWAS$`phewas phenotype`)
df_GWAS_102$UP_DiseaseName <- toupper(df_GWAS_102$MAPPED_TRAIT)

PheWAS_GWAS <- merge(df_GWAS_102, df_PheWAS, by.x = "UP_DiseaseName", by.y = "UP_DiseaseName")

#PheWAS&GWASとClinver&OMIMデータは名前で合体
Clinver_OMIM$UP_DiseaseName <- toupper(Clinver_OMIM$DiseaseName)

Clinver_OMIM_PheWAS_GWAS <- merge(Clinver_OMIM, PheWAS_GWAS, by.x = "UP_DiseaseName", by.y = "UP_DiseaseName")

#必要なところだけ抽出
result <- lapply(c("OMIM_GeneID", "Clinver_GeneID", "GWAS_GeneID", "PheWAS_GeneID"), 
            function(x){tmp <- strsplit(Clinver_OMIM_PheWAS_GWAS[[x]], ",")
            tmp2 <- sapply(tmp, function(x){x[x != ""]})
            return(tmp2)})
list_disease_genes <- sapply(1:dim(Clinver_OMIM_PheWAS_GWAS)[1], function(y){unique(c(result[[1]][[y]],result[[2]][[y]] ,result[[3]][[y]], result[[4]][[y]]))})
names(list_disease_genes) <- Clinver_OMIM_PheWAS_GWAS$DiseaseName


#必要なところだけ抽出
result <- lapply(c("GWAS_GeneID", "PheWAS_GeneID"), 
                 function(x){tmp <- strsplit(PheWAS_GWAS[[x]], ",")
                 tmp2 <- sapply(tmp, function(x){x[x != ""]})
                 return(tmp2)})
list_disease_genes_PheWAS_GWAS <- sapply(1:dim(PheWAS_GWAS)[1], function(y){unique(c(result[[1]][[y]],result[[2]][[y]]))})
names(list_disease_genes_PheWAS_GWAS) <- PheWAS_GWAS$`phewas phenotype`

save(list_disease_genes_PheWAS_GWAS, file = "./Analysis_230701/0.data_make/causal_genes/list_disease_genes_PheWAS_GWAS.dat")
