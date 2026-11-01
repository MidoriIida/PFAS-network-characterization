###############################
#                             #
#       ID conversion         #
#                             #
###############################
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("biomaRt")
library(biomaRt)
library(openxlsx)
library(enrichR)

### 1. load ####################
in_d <- "./path/to/intermediate/data/0.data_make/causal_genes/"
in_f <- list.files(in_d)

list_causal_genes <- lapply(paste0(in_d, in_f), function(x){get(load(x))})
names(list_causal_genes) <- gsub( "Causal_genes_","",gsub(".dat", "", in_f))

### 2. data set ####################
db <- useMart("ensembl")
hd <- useDataset("hsapiens_gene_ensembl", mart = db)

filters <- listFilters(hd)
filters[grep("symbol", filters[,1]),]

attributes <- listAttributes(hd)
attributes[grep("entrez", attributes[,1]),]

### 3. conversion ####################
My_getBM <- 
function(x){getBM(attributes = c("hgnc_symbol", "entrezgene_id", "description"), 
                  filters = "hgnc_symbol", values = x, mart = hd, useCache = FALSE)}

My_getBM_OMIM <- 
  function(x){getBM(attributes = c("hgnc_symbol", "entrezgene_id", "description"), 
                    filters = "entrezgene_id", values = x, mart = hd, useCache = FALSE)}



list_GeneIDs <- lapply(list_causal_genes, function(y){lapply(y[2:4], My_getBM)})
list_GeneIDs_OMIM  <-  lapply(list_causal_genes, function(y){lapply(y[1], My_getBM_OMIM)})
  
### 4. 合体 ####################
for(i in 1:length(list_GeneIDs_OMIM)){
  list_GeneIDs[[i]]$OMIM<- list_GeneIDs_OMIM[[i]]$OMIM
}
save(list_GeneIDs, file = "./path/to/intermediate/data/0.data_make/causal_genes/list_GeneIDs.dat")


### 5. unlist and unique ####################
list_Disease_GeneIDs <- sapply(list_GeneIDs, function(x){unique(unlist(sapply(x, function(x){x$entrezgene_id })))})
save(list_Disease_GeneIDs, file = "./path/to/intermediate/data/0.data_make/causal_genes/list_Disease_GeneIDs.dat")


### 6. 出力 ####################
for (i in names(list_GeneIDs)) {
  #write.xlsx(list_GeneIDs[[i]],          file = paste0("./path/to/intermediate/data/0.data_make/causal_genes/list_GeneIDs_", i, ".xlsx"))
  write.table(list_Disease_GeneIDs[[i]], file = paste0("./path/to/intermediate/data/0.data_make/causal_genes/list_Disease_GeneIDs_", i, ".txt"), quote = F, row.names = F, col.names = F)
}

### 4. check #####################
dbs   <- listEnrichrDbs()
my_genes <- function(y){unique(unlist(sapply(y, function(x){x$hgnc_symbol})))}
list_enr <- sapply(sapply(list_GeneIDs, my_genes), function(x){enrichr(x,"KEGG_2021_Human")})


