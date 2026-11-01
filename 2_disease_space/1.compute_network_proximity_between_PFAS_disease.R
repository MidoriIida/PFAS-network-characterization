#######################
#   Randomization     #
#######################
library(igraph)
library(parallel)

# random seed
set.seed(123)

# get number of available cores
numCores <- detectCores()

# function for get candidates hash
get_candidates_hash <- function(node_degrees){
  sorted_node_degrees <- sort(node_degrees)
  candidates_hash <- new.env()
  for(v in names(sorted_node_degrees)){
    idx <- which(names(sorted_node_degrees) == v)
    end_idx <- min(idx + 50, length(sorted_node_degrees))
    start_idx <- max(idx - 50, 1)
    if(end_idx < idx + 50){
      start_idx <- start_idx - (idx + 50 - end_idx)
    }
    if(idx - 50 < 1){
      end_idx <- end_idx - (idx - 50 - 1)
    }
    
    sorted_node_degrees[start_idx:end_idx]
    k_min <- min(sorted_node_degrees[start_idx:end_idx])
    k_max <- max(sorted_node_degrees[start_idx:end_idx])
    
    candidates_hash[[v]] <- sample(names(sorted_node_degrees[sorted_node_degrees <= k_max & sorted_node_degrees >= k_min]))
  }
  return(candidates_hash)
}

# function for computing dc
compute_dc <- function(mtx, set_S, set_T){
  distance_set <- c()
  for(t in node_set_T){
    distance_set <- c(distance_set, min(distance_matrix[set_S, t]))
  }
  return(mean(distance_set))
}

# function for computing null dc
compute_dc_parallel <- function(n, distance_matrix, node_set_S, node_set_T, hash) {
  # get 2 randomly selected groups
  candicate_nodes <- mget(node_set_S, envir = hash)
  niter <- 0
  while(niter < 20){
    set_S_random <- unlist(lapply(candicate_nodes, function(x) sample(x, 1)))
    if(length(unique(set_S_random)) == length(node_set_S)){
      break
    }
    niter <- niter + 1
  }
  
  candicate_nodes <- mget(node_set_T, envir = hash)
  niter <- 0
  while(niter < 20){
    set_T_random <- unlist(lapply(candicate_nodes, function(x) sample(x, 1)))
    if(length(unique(set_T_random)) == length(node_set_T)){
      break
    }
    niter <- niter + 1
  }
  
  # compute dc
  dc_random <- compute_dc(distance_matrix, set_S_random, set_T_random)
  return(dc_random)
}

### 1.load data ############################
# PPI network (LCC) — produced by 0_data_preparation/1_ppi/1.PPI_data_make.R
load("./path/to/intermediate/data/g_comp.dat")
g <- g_comp
rm(g_comp)

# compute distance matrix of the network
distance_matrix <- distances(g)

# compute node degrees
node_degrees <- degree(g)

# load chemical-gene interactions (produced by 0_data_preparation/2_pfas_genes/1.PFAS_genes_load.R)
load("./path/to/intermediate/data/list_PFAS_genes.dat")
chemical_gene_pairs <- list_PFAS_genes$Merged
#rm(list_PFAS_genes)

# load disease-gene interactions
in3 =  "./path/to/intermediate/data/0.make_data/human_disease_genes/human_disease_genes_230714_iida.dat"
load(in3) #"human_disease_genes/human_disease_genes_230714_iida.dat"

# get disease and chemical names
disease_names <- names(list_disease_genes)
chemical_casrn <- unique(chemical_gene_pairs$CASRN)

### 2. compute distance between chemical-disease pair #########################
# get random candidate nodes for null hypothesis
random_candicates_hash <- get_candidates_hash(node_degrees)

# compute chemical-disease pair
list_dist_mat <-  vector(mode = "list", length = length(chemical_casrn))
names(list_dist_mat) <- chemical_casrn
dist_mat      <- matrix(data=NA, nrow = length(disease_names), ncol = 6)
rownames(dist_mat) <- disease_names
colnames(dist_mat) <- c("chemical_tag", "disease", "dc_real", "dc_mean", "dc_std", "z_dc")

for(chemical in chemical_casrn ){ #c("103055-07-8")
  chemical_tag <- paste(unique(chemical_gene_pairs[chemical_gene_pairs["CASRN"] == chemical, c("CASRN", "ChemicalID","ChemicalName")]), collapse = "||")
  node_set_T <- intersect(as.character(chemical_gene_pairs[chemical_gene_pairs["CASRN"] == chemical, "GeneID"]), V(g)$name)
  for(i in 1:length(disease_names)){ #
    disease <- disease_names[i]
    node_set_S <- intersect(as.character(list_disease_genes[disease][[1]]), V(g)$name)
    # compute dc between S and T
    dc_real <- compute_dc(distance_matrix, node_set_S, node_set_T)
    # compute dc between two randomly selected groups
    dc_random_set <- unlist(mclapply(1:1000, compute_dc_parallel, mc.cores = numCores, distance_matrix, node_set_S, node_set_T, random_candicates_hash))
    # compute z-score
    dc_mean <- mean(dc_random_set)
    dc_std <- sd(dc_random_set)
    z_dc <- (dc_real - dc_mean) / dc_std
    # print
    #cat(paste(chemical_tag, disease, dc_real, dc_mean, dc_std, z_dc, sep="\t"), "\n", sep="")
    dist_mat[i,] <- c(chemical_tag, disease, dc_real, dc_mean, dc_std, z_dc)
    print(i)
  }
  list_dist_mat[[chemical]] <- as.data.frame(dist_mat)
}
save(list_dist_mat, file = "./path/to/intermediate/data/list_dist_mat_disease_PFAS.dat")
#system("python3.9 notification.py")


