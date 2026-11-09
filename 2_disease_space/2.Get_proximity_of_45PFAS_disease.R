##############################################
#   Make proximity matrix and UMAP/t-SNE    #
#       45 PFAS × 148 disease space         #
##############################################
set.seed(123)  # for reproducibility
library(tidyverse)
library(ggplot2)
library(umap)
library(dplyr)
library(lemon)
library(Rtsne)

# define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

### 1. load #################################################################
# list_dist_mat produced by 1.compute_network_proximity_between_PFAS_disease.R
load(file = "./path/to/intermediate/data/list_dist_mat_disease_PFAS.dat")

### 2. Make z matrix ########################################################
df_dist <- as.data.frame(sapply(list_dist_mat, function(x){as.numeric(x[,"z_dc"])}))
rownames(df_dist) <- list_dist_mat[[1]]$disease

save(df_dist, file = "./path/to/intermediate/data/df_dist_45PFAS_disease.dat")

### 3. Check ################################################################
PFOA <- df_dist[,"335-67-1"]
PFOS <- df_dist[,"1763-23-1"]
names(PFOA) <- names(PFOS) <- rownames(df_dist)
head(sort(PFOA[PFOA < -0.50]))
head(sort(PFOS[PFOS < -0.50]))

### 4. Remove duplicate diseases ############################################
df_z <- df_dist[which(!duplicated(sapply(strsplit(rownames(df_dist), "_"), function(x){x[length(x)]}))),]

# Check for missing values
sum(is.na(df_z))

# Impute missing values
df_z_0   <- df_z %>% mutate(across(everything(), ~replace(., is.nan(.), 0))) %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))

# Remove columns with any missing values
df_z_non <- df_z[, which(apply(df_z, 2, function(x){sum(is.na(x))}) == 0)]

### 5. Make similarity matrix ###############################################
df_sim_norm_0   <- apply(df_z_0 * (-1), 2, min_max_norm)
df_sim_norm_non <- apply(df_z_non * (-1), 2, min_max_norm)

### 6. Convert similarity to distance (1 minus similarity) ##################
distance_matrix_0   <- 1 - df_sim_norm_0
distance_matrix_non <- 1 - df_sim_norm_non

### 7. Run UMAP #############################################################
list_mat <- list(
  sim0    = df_sim_norm_0,
  simnon  = df_sim_norm_non,
  dist0   = distance_matrix_0,
  distnon = distance_matrix_non
)

list_umap <- lapply(list_mat, function(x){umap(t(x), n_components = 2)})

save(list_umap, file = "./path/to/intermediate/data/list_umap_45PFAS_148disease.dat")
