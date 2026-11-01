##########################
# Make proximity matrix  #
##########################
library(openxlsx)

### 1. load #########################
load(file = "./path/to/intermediate/data/list_dist_mat_disease_PFAS.dat")

### 2. get proximity ################
df_dist <- as.data.frame(sapply(list_dist_mat, function(x){as.numeric(x[,6])}))
rownames(df_dist) <- list_dist_mat[["103055-07-8"]]$disease

save(df_dist, file = "./path/to/intermediate/data/df_dist_45PFAS_disease.dat")

### 3. check ########################
PFOA <- df_dist[,"335-67-1"]
PFOS <- df_dist[,"1763-23-1"]

names(PFOA) <- names(PFOS) <- rownames(df_dist)

head(sort(PFOA[PFOA < -0.50]))
head(sort(PFOS[PFOS < -0.50]))
