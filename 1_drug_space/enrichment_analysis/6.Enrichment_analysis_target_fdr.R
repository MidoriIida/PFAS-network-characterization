##########################################
#   Plot according to 900 drug space     #
#              45 PFAS                   #
#   Cluster interpretation (target)       #
##########################################
set.seed(123)  # for reproducibility
library(ggplot2)
library(dplyr)
library(umap)
library(dplyr)
library(lemon)
library(reshape2)
library(cowplot)

### 1. load ###################################################################
load("./path/to/intermediate/data/list_cluster.dat")

# Drug target class enrichment p-values (computed by collaborators; Supplementary Data 8)
Target      <- read.csv("../../data/Supplementary_Data_8_Enrichment_Target.txt", row.names = 1)

# FDR correction
# Convert data to matrix
p_matrix <- as.matrix(Target)

# Flatten p-values to vector → FDR correction → reshape to matrix
p_fdr <- matrix(
  p.adjust(as.vector(p_matrix), method = "fdr"),
  nrow = nrow(p_matrix),
  ncol = ncol(p_matrix)
)

# Restore row and column names
rownames(p_fdr) <- rownames(p_matrix)
colnames(p_fdr) <- colnames(p_matrix)


### 2. cluster ################################################################
clusters <- list_cluster$clusters

c1 <- names(clusters[clusters == 1])
c2 <- names(clusters[clusters == 2])

cas <- sapply(strsplit(gsub("^X", "", colnames(Target)), "\\.\\."), function(x){gsub("\\.", "-",x[1])})
names(cas) <- cas

### 3. Relationship between hierarchical clustering and drug targets ######################################
c1_Target <- Target[cas %in% c1,]
c2_Target <- Target[cas %in% c2,]

c1_Target_sig <- apply(c1_Target, 1, function(x){sum(na.omit(x) < 0.05)})
c2_Target_sig <- apply(c2_Target, 1, function(x){sum(na.omit(x) < 0.05)})

# Create data frame
Target_sig <- cbind(c1_Target_sig, c2_Target_sig)
x <- Target_sig[which(rowSums(Target_sig) != 0),]
df_Target_sig <- melt(x[order(x[,2]),])
df_Target_sig$group <- gsub("c1_Target_sig", "cluster 1", df_Target_sig$Var2)
df_Target_sig$group <- gsub("c2_Target_sig", "cluster 2", df_Target_sig$group)
levels(df_Target_sig$group) <- c("cluster 1", "cluster 2")

colnames(df_Target_sig) <- c("Target", "cluster", "count", "group")

### 4. Create bar plot ######################################
g_Target <-
  ggplot(df_Target_sig, aes(x = Target, y = count, fill = group)) +  # x-axis: category, y-axis: count
  geom_bar(stat = "identity", position = "dodge") +  # bar height based on data values
  labs(#title = "Significant Target Counts by cluster", 
    x = "Drug target", y = "Count",  fill = "Cluster") +
  theme_bw() 

### 5. Adjust bar plot (horizontal) ######################################
g_Target_h <- g_Target  + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                                legend.position = c(0.95, 0.95),  # place legend in the upper right inside the plot
                                legend.justification = c("right", "top"))  # anchor legend to upper-right corner
g_Target_h

### 5. Adjust bar plot (vertical) ######################################
g_Target_v <- g_Target + coord_flip()  + 
  theme( legend.position = c(1, 0.05),  # place legend in the upper right inside the plot
         legend.justification = c("right", "bottom"))
g_Target_v

### 6. save ##########################################################
list_graph_enriched_target <- list(h = g_Target_h, v = g_Target_v, n= g_Target)
save(list_graph_enriched_target, file = "./path/to/intermediate/data/6.Enrichment_analysis/list_graph_enriched_target_fdr.dat")

