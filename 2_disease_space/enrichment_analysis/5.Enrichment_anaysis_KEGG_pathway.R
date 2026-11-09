########################
# enrichment analysis  #
#  for signaling kegg  #
########################

MyFisher <- function(my_h, h_in_category, all_h){
  # Contingency counts
  a <- sum(my_h %in% h_in_category)                     # hits (overlap)
  b <- sum(!(my_h %in% h_in_category))                  # non-hits in our set
  c <- sum(h_in_category %in% setdiff(all_h, my_h))     # in background but not in our set
  d <- sum(!(all_h %in% union(my_h, h_in_category)))    # neither in our set nor in category
  
  # Contingency table
  contingency <- matrix(c(a, b, c, d), nrow = 2,
                        dimnames = list(MySet = c("InCat", "OutCat"),
                                        Background = c("InMySet", "OutMySet")))
  
  res <- fisher.test(contingency, alternative = "greater")  # one-sided test (enrichment)
  return(res$p.value) 
}

### 1. load data ##########################
# Load proximity list (produced by 2_disease_space/2.Get_proximity_of_45PFAS.R)
load(file = "./path/to/intermediate/data/df_dist_45PFAS_disease.dat") #df_dist

# Build H-number lookup from disease rownames (format: "H00003_Acute myeloid leukemia")
df_H_numbers <- data.frame(
  targets = rownames(df_dist),
  H_name  = sub("_.*", "", rownames(df_dist)),
  stringsAsFactors = FALSE
)

# Load KEGG BRITE pathway classification (downloaded from https://www.genome.jp/brite/br08402)
lines <- readLines("./path/to/KEGG/data/br08402_KEGG_HumanDiseases_pathway.txt")

### 2. make H list  ########################
# List to store results
list_brite <- list()
current_b <- NULL

for (line in lines) {
  if (grepl("^B\\s+", line)) {
    # B rows
    current_b <- trimws(sub("^B\\s+", "", line))
    list_brite[[current_b]] <- character(0)
  } else if (grepl("^C\\s+H[0-9]{5}", line) && !is.null(current_b)) {
    # C rows
    h_number <- sub("^C\\s+(H[0-9]{5}).*", "\\1", line)
    list_brite[[current_b]] <- c(list_brite[[current_b]], h_number)
  }
}


### 3. Make dist mat ##############################
df_z <- df_dist[which(!duplicated(sapply(strsplit(rownames(df_dist), "_"), function(x){x[length(x)]}))),]

### 4. Prepare for Fisher's exact test ########################
all_h <- unique(unlist(list_brite)) # background set (all unique H numbers)

list_my_h <- sapply(1:dim(df_z)[2], function(i){na.omit(df_H_numbers$H_name[(df_z[,i] < -0.50)])})
names(list_my_h) <- colnames(df_z)

### 5. Fisher's exact test across all PFAS ###############
my_h <- unique(unlist(list_my_h))

df_res <- sapply(list_brite, function(h_in_category){
  MyFisher(my_h = my_h, 
           h_in_category = h_in_category, 
           all_h = all_h)}) 

df_res_fdr <- p.adjust(df_res, method = "fdr")

### 5. Fisher's exact test for each PFAS ########################
df_results <- sapply(list_brite, function(h_in_category){
  sapply(list_my_h, function(my_h){MyFisher(my_h = my_h, 
                                            h_in_category = h_in_category, 
                                            all_h = all_h)}) # apply Fisher's exact test
})

# FDR correction
# Convert data to matrix
p_matrix <- as.matrix(df_results)

# Flatten p-values to vector → FDR correction → reshape to matrix
p_fdr <- matrix(
  p.adjust(as.vector(p_matrix), method = "fdr"),
  nrow = nrow(p_matrix),
  ncol = ncol(p_matrix)
)

# Restore row and column names
rownames(p_fdr) <- rownames(p_matrix)
colnames(p_fdr) <- colnames(p_matrix)

# Convert back to data frame (same shape as original)
df_fdr <- as.data.frame(p_fdr)

### 5. check the result ########################
list_res <- sapply(1:45, function(i){which((df_results[i,] < 0.05))})

colSums(df_results < 0.05) # not many significant results
colSums(df_fdr < 0.05) # none

save(df_results, file="./path/to/intermediate/data/df_results_enrichment_KEGG_pathway.dat")
