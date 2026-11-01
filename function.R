####################################
#      shortest_path_length        #
#         within module A          #
####################################
within_shortest_path_length <- function(g_comp, dA, dB){
  mat_dst <- igraph::distances(g_comp, v  = as.character(dA), to = as.character(dB)) #d(s, t)
  mindst  <- apply(mat_dst, 1, function(x){min(x[x != 0])})
  return(sum(mindst)/(length(mindst)))
}

####################################
#      shortest_path_length        #
#      between module A & B        #
####################################
between_shortest_path_length <- function(g_comp, dA, dB){
  mat_dst <- igraph::distances(g_comp, v  = as.character(dA), to = as.character(dB)) #d(s, t)
  min_A <- apply(mat_dst, 1, function(x){min(x)})
  min_B <- apply(mat_dst, 2, function(x){min(x)})
  ans <- sum(c(min_A, min_B))/length(c(min_A, min_B))
  return(ans)
}

####################################
#       separation_measure         #
#      between module A & B        #
####################################
separation_measure <- function(g_comp, dA, dB){
  dAA <- within_shortest_path_length(g_comp, dA, dA)
  dBB <- within_shortest_path_length(g_comp, dB, dB)
  dAB <- between_shortest_path_length(g_comp, dA, dB)
  sAB <- dAB - (dAA + dBB)*0.5
  return(sAB)
}

####################################
#        closest_measure           #
#      between module A & B        #
####################################
closest_measure <- function(g_comp, dA, dB){
  mat_dst <- igraph::distances(g_comp, v  = as.character(dA), to = as.character(dB)) 
  min_B   <- apply(mat_dst, 2, function(x){min(x)})
  CsAB    <- sum(min_B)/length(min_B)
  return(CsAB)
}

####################################
#        shortest_measure          #
#      between module A & B        #
####################################
shortest_measure <- function(g_comp, dA, dB){
  mat_dst <- igraph::distances(g_comp, v  = as.character(dA), to = as.character(dB)) 
  ShsAB <- sum(mat_dst)/length(mat_dst)
  return(ShsAB)
}

####################################
#         kernel_measure          #
#      between module A & B        #
####################################
kernel_measure <- function(g_comp, dA, dB){
  mat_dst <- igraph::distances(g_comp, v  = as.character(dA), to = as.character(dB)) 
  A <- log10(apply(mat_dst, 1, function(x){mean(exp(-1*(x+1)))}))
  B <- log10(apply(mat_dst, 2, function(x){mean(exp(-1*(x+1)))}))
  kdAB <- -1/(length(A) + length(B))*(sum(A)+sum(B))
  return(kdAB)
}

####################################
#         centre_measure           #
#      between module A & B        #
####################################
centre_measure <- 
  function(g_comp, dA, dB){
    mat_dst_A  <- igraph::distances(g_comp, v  = as.character(dA), to = as.character(dA)) 
    dAu        <- apply(mat_dst_A, 1, sum)
    argmin_dA  <- which(dAu == min(dAu))
    
    mat_dst_B  <- igraph::distances(g_comp, v  = as.character(dB), to = as.character(dB)) 
    dBu        <- apply(mat_dst_B, 1, sum)
    argmin_dB  <- which(dBu == min(dBu))
    
    if(length(argmin_dA)==1 & length(argmin_dB)==1){
      dccAB <- igraph::distances(g_comp, v  = names(argmin_dA), to = names(argmin_dB)) 
    }else{
      dccAB <- mean(igraph::distances(g_comp, v  = names(argmin_dA), to = names(argmin_dB)))
    }
    
    return(dccAB)
  }

####################################
#   shortest_path_lengthを求める   #
####################################
shortest_path_length <- function(g_comp, dA, dB){
  dst     <- igraph::distances(g_comp, v  = names(V(g_comp)) %in% dA, to = names(V(g_comp)) %in% dB ) #d(s, t)
  y <- mean(apply(dst, 1, min))
  return(y)
}

##############################
#       Binを作る            #  
##############################
#cut degree vectors with optional number of nords
Binning <- function(deg, nums){
  y <- seq(0:max(deg))
  x <- cut(deg, breaks = y,  include.lowest = FALSE, right = FALSE)
  table_x <- (table(x))
  seqs <- c(1)
  while(TRUE){
    i <- seqs[length(seqs)]
    s <- 0
    while (s <= nums) {
      s <- s + table_x[y[i]]
      i <- i + 1
      if(y[i] > length(table_x)) break
      #print(i)
    }
    
    seqs <- c(seqs, i)
    if(y[i] > length(table_x)) break
  }
  
  test <- cut(deg, breaks = c(seqs), include.lowest = FALSE, right = FALSE)
  table_test <- table(test)
  if(table_test[length(table_test)] <= nums){
    seqs <- seqs[-(length(seqs)-1)]
  }
  
  cut_deg <- cut(deg, breaks = seqs, include.lowest = FALSE, right = FALSE)
  
  return(cut_deg)
}

##############################
#   ランダムにノードを抽出   #
##############################
ChooseNodes <- function(cut_degs, origene){
  list_cut_degs <- tapply(names(cut_degs), cut_degs, list)
  
  degs_dA <- na.omit(cut_degs[names(origene)]) #graphにあるもののみ
  freq_degs_dA <- table(degs_dA)
  
  ramd_dA <- c()
  for(k in 1:length(freq_degs_dA)){
    tmp <- freq_degs_dA[k]
    smp <- sample(list_cut_degs[[names(tmp)]], size=tmp, replace = FALSE)
    ramd_dA <- c(ramd_dA,smp)
  }
  
  ramd_dA <- as.character(ramd_dA)
  
  return(ramd_dA)
}
