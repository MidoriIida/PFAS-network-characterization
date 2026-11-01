########################
# enrichment analysis  #
#  for signaling kegg  #
########################

MyFisher <- function(my_h, h_in_category, all_h){
  # 集計
  a <- sum(my_h %in% h_in_category)                     # ヒットしている数（重なり）
  b <- sum(!(my_h %in% h_in_category))                  # 自分のセットのうち非ヒット
  c <- sum(h_in_category %in% setdiff(all_h, my_h))     # 他のHにあり、my_hに無い
  d <- sum(!(all_h %in% union(my_h, h_in_category)))    # 全体のうち無関係な部分
  
  # 分割表
  contingency <- matrix(c(a, b, c, d), nrow = 2,
                        dimnames = list(MySet = c("InCat", "OutCat"),
                                        Background = c("InMySet", "OutMySet")))
  
  res <- fisher.test(contingency, alternative = "greater")  # 片側検定（enrichment）
  return(res$p.value) 
}

### 1. load data ##########################
#疾患のリストを読み込む
load(file = "/Users/iida/Documents/Lab_Iida/2023_Iida/Interactome/Human_PPI_PFAS/Analysis_250707/5.ID_conversionList_KEGG_disease/df_H_numbers.dat") #df_H_numbers

#Proximityのリストを読み込む
load(file = "./Analysis_250707/df_dist_45PFAS_disease.dat") #df_dist

# テキストファイルを読み込む（必要に応じて file.choose() などを使用）
lines <- readLines("/Users/iida/Documents/Lab_Iida/Orignal_data/KEGG/KEGG_2025/br08403_KEGG BRITE_ Human Diseases in ICD-11 Classification.txt")

### 2. make H list  ########################
# 結果を格納するリスト
list_brite <- list()
current_b <- NULL

for (line in lines) {
  if (grepl("^[0-9]{2} ", line)) {
    # B行から名前を抽出（例: "nt06031  Citrate cycle and pyruvate metabolism"）
    current_b <- line
    list_brite[[current_b]] <- character(0)
  } else if (grepl("^H[0-9]{5}", line) && !is.null(current_b)) {
    # C行から要素を抽出（例: "H01225  D-2-hydroxyglutaric aciduria"）
    #value <- trimws(sub("^C\\s+", "", line))
    h_number <- sub("(H\\d+).*", "\\1", line)
    list_brite[[current_b]] <- c(list_brite[[current_b]], h_number)
  }
}

save(list_brite, file = "./Analysis_250707/6.Enrichment_anaysis/list_brite.dat")
       
### 3. Make dist mat ##############################
df_z <- df_dist[which(!duplicated(sapply(strsplit(rownames(df_dist), "_"), function(x){x[length(x)]}))),]

### 4. Fishers exact testの準備  ########################
all_h <- unique(unlist(list_brite)) #背景集合（全体に含まれるユニークなH番号）

list_my_h <- sapply(1:dim(df_z)[2], function(i){na.omit(df_H_numbers$H_name[(df_z[,i] < -0.50)])})
names(list_my_h) <- colnames(df_z)

### 5. 丸っと全部のPFAS で　Fishers exact test###############
my_h <- unique(unlist(list_my_h))

df_res <- sapply(list_brite, function(h_in_category){
MyFisher(my_h = my_h, 
         h_in_category = h_in_category, 
         all_h = all_h)}) 

df_res_fdr <- p.adjust(df_res, method = "fdr")

### 5. 各物質ごとにFishers exact test########################
df_results <- sapply(list_brite, function(h_in_category){
  sapply(list_my_h, function(my_h){MyFisher(my_h = my_h, 
                                            h_in_category = h_in_category, 
                                            all_h = all_h)}) # フィッシャー検定の適用
})

#FDR補正
# データを行列に変換
p_matrix <- as.matrix(df_results)

# p値を1次元ベクトル化 → FDR補正 → 行列に戻す
p_fdr <- matrix(
  p.adjust(as.vector(p_matrix), method = "fdr"),
  nrow = nrow(p_matrix),
  ncol = ncol(p_matrix)
)

# 行名・列名を戻す
rownames(p_fdr) <- rownames(p_matrix)
colnames(p_fdr) <- colnames(p_matrix)

# データフレームに変換（元と同じ形）
df_fdr <- as.data.frame(p_fdr)


### 5. check the result ########################
list_res <- sapply(1:45, function(i){which((df_results[i,] < 0.05))})
colSums(df_results < 0.05) #あまり変わり映えしないようだ...
colSums(df_fdr < 0.05)     #いくつかは有意に異なる結果として出るけど、結構少なめ

save(df_results, file="/Users/iida/Documents/Lab_Iida/2023_Iida/Interactome/Human_PPI_PFAS/Analysis_250707/6.Enrichment_anaysis/df_results_enrichment_KEGG_ICD11.dat")
save(df_fdr,     file="./Analysis_250707/df_results_enrichment_KEGG_ICD11_fdr.dat")
