#########################
#   Plot 45 PFAS into   #
#   148 disease space   #
#diseaseに重複がない状態#
#########################
set.seed(123)  # 結果の再現性のため
library(readxl)
library(openxlsx)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(umap)
library(dplyr)
library(lemon)
library(Rtsne)
#define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

### 1. load ##############################################################
load("./Analysis_250707/df_dist_45PFAS_disease.dat")

in_d  <- "/path/to/your/local/directory/original_data"
in_f3 <- "hsPPARalpha_PFASs_wash_descriptors_Sscore_IC50_data.xlsx"

discriptors <- read.xlsx(paste0(in_d, in_f3))

in_f <- "/path/to/your/local/directory/PFOA_abbrivation_list_pfasmaster-2023-01-31.xlsx"
abbri <- read.xlsx(in_f)

### 2. Make z matrix ####################################################
df_z <- df_dist[which(!duplicated(sapply(strsplit(rownames(df_dist), "_"), function(x){x[length(x)]}))),]

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
list_mat <- list(
  sim0    = df_sim_norm_0,
  simnon  = df_sim_norm_non,
  dist0   = distance_matrix_0,
  distnon = distance_matrix_non
)

list_umap <- lapply(list_mat, function(x){umap(t(x), n_components = 2)})

save(list_umap, file = "./Analysis_250707/list_umap_45PFAS_148disease.dat")

### 6. PCAの実行##############################################################
list_pca <- lapply(list_mat, function(x){prcomp(t(x), scale. = T)})

biplot(list_pca[[1]])
biplot(list_pca[[3]])

### 7. t-sneの実行##############################################################
list_tsne<- lapply(list_mat, function(x){Rtsne(t(x),pca=FALSE,perplexity=8,theta=0)})
list_pcatsne<- lapply(list_mat, function(x){Rtsne(t(x),pca=TRUE,perplexity=8,theta=0)})

tsne_value<-as.data.frame(list_tsne[[1]]$Y)
p<-ggplot(tsne_value,aes(V1,V2))+geom_point()+
  stat_ellipse(type = "norm",linetype=1)
p

tsne_value<-as.data.frame(list_tsne[[2]]$Y)
p<-ggplot(tsne_value,aes(V1,V2))+geom_point()+
  stat_ellipse(type = "norm",linetype=1)
p

tsne_value<-as.data.frame(list_tsne[[3]]$Y)
p<-ggplot(tsne_value,aes(V1,V2))+geom_point()+
  stat_ellipse(type = "norm",linetype=1)
p

tsne_value<-as.data.frame(list_pcatsne[[4]]$Y)
p<-ggplot(tsne_value,aes(V1,V2))+geom_point()+
  stat_ellipse(type = "norm",linetype=1)
p

