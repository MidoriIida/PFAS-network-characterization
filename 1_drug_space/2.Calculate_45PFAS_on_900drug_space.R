##########################################
#   Plot according to 900 drug space     #
#              45 PFAS                   #
##########################################
set.seed(123)  # 結果の再現性のため
library(readxl)
library(openxlsx)
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
in_f <- "/path/to/your/local/directory/summary_PFAS_appdrug_iidaData_20230716.xlsx"
PFAS_drug <- read_excel(in_f, col_names = T)

### 2. Make z matrix ####################################################
df_z           <- as.data.frame(matrix(PFAS_drug$`sd(dc_null)...6`, ncol = 900))
rownames(df_z) <- sapply(strsplit(unique(PFAS_drug$ID), "\\|"), function(x){x[1]})
colnames(df_z) <- PFAS_drug$Name[1:900]

# データセットの欠損値の確認
sum(is.na(df_z))
    
# データセットの欠損値の補完
df_z_0　 <- df_z %>% mutate(across(everything(), ~replace(., is.nan(.), 0))) %>% 
  mutate(across(everything(), ~replace(., is.na(.), 0)))       

# データセットの欠損値がある場合、その列は削除
df_z_non <- df_z[,which(apply(df_z, 2, function(x){sum(is.na(x))}) == 0)]

### 3. Make similarity matrix #################################################
df_sim_norm_0     <- apply(df_z_0*(-1), 2, min_max_norm)
df_sim_norm_non   <- apply(df_z_non*(-1), 2, min_max_norm)

### 4. 類似度から距離への変換（1から類似度を引く）#############################
distance_matrix_0   <- 1 - df_sim_norm_0
distance_matrix_non <- 1 - df_sim_norm_non

### 5. UMAPの実行##############################################################
umap_sim_result_0   <- umap(df_sim_norm_0, n_components = 2)
umap_sim_result_non <- umap(df_sim_norm_non, n_components = 2)

umap_dist_result_0   <- umap(distance_matrix_0, n_components = 2)
umap_dist_result_non <- umap(distance_matrix_non, n_components = 2)

list_umap <- list(sim0    = umap_sim_result_0, 
                  simnon  = umap_sim_result_non,
                  dist0   = umap_dist_result_0,
                  distnon = umap_dist_result_non)

save(list_umap, file = "./Analysis_240819/list_umap_PFAS_drug.dat")













