##########################################
#   Plot according to 900 drug space     #
#              45 PFAS                   #
##########################################
set.seed(123)  # for reproducibility
library(tidyverse)
library(ggplot2)
library(umap)
library(dplyr)
library(lemon)
#define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

### 1. load ##############################################################
# list_dist_mat_drug_PFAS produced by 1_drug_space/1.compute_network_proximity_between_PFAS_drug.R
load("./path/to/intermediate/data/list_dist_mat_drug_PFAS.dat")

### 2. Make z matrix ####################################################
# rows = drugs, cols = PFAS (CASRN); then transpose to rows = PFAS, cols = drugs
df_z <- as.data.frame(sapply(list_dist_mat, function(x) as.numeric(x[, "z_dc"])))
rownames(df_z) <- list_dist_mat[[1]]$drug
df_z <- t(df_z)
rownames(df_z) <- names(list_dist_mat)
df_z <- as.data.frame(df_z)

# Check for missing values in the dataset
sum(is.na(df_z))
    
# Impute missing values in the dataset
df_z_0　 <- df_z %>% mutate(across(everything(), ~replace(., is.nan(.), 0))) %>% 
  mutate(across(everything(), ~replace(., is.na(.), 0)))       

# Remove columns with missing values if present
df_z_non <- df_z[,which(apply(df_z, 2, function(x){sum(is.na(x))}) == 0)]

### 3. Make similarity matrix #################################################
df_sim_norm_0     <- apply(df_z_0*(-1), 2, min_max_norm)
df_sim_norm_non   <- apply(df_z_non*(-1), 2, min_max_norm)

### 4. Convert similarity to distance (1 minus similarity) #############################
distance_matrix_0   <- 1 - df_sim_norm_0
distance_matrix_non <- 1 - df_sim_norm_non

### 5. Run UMAP ##############################################################
umap_sim_result_0   <- umap(df_sim_norm_0, n_components = 2)
umap_sim_result_non <- umap(df_sim_norm_non, n_components = 2)

umap_dist_result_0   <- umap(distance_matrix_0, n_components = 2)
umap_dist_result_non <- umap(distance_matrix_non, n_components = 2)

list_umap <- list(sim0    = umap_sim_result_0, 
                  simnon  = umap_sim_result_non,
                  dist0   = umap_dist_result_0,
                  distnon = umap_dist_result_non)

save(list_umap, file = "./path/to/intermediate/data/list_umap_PFAS_drug.dat")













