##########################################
#   Plot according to 900 drug space     #
#              45 PFAS                   #
#      Cluster interpretation (target)   #
#   lemon-free version                   #
##########################################
set.seed(123)

library(ggplot2)   # load ggplot2 first to prevent lemon from overriding +.gg
library(dplyr)
library(reshape2)
library(cowplot)
# library(lemon)  # excluded (not used)

### 1. load ###################################################################
load("./path/to/intermediate/data/list_cluster.dat")

# Supplementary_Data_8 produced by 1_drug_space/enrichment_analysis/5.Enrichment_analysis_target.R
Target <- read.csv("../../data/Supplementary_Data_8_Enrichment_Target.txt", row.names = 1)

### 2. cluster ################################################################
clusters <- list_cluster$clusters

c1 <- names(clusters[clusters == 1])
c2 <- names(clusters[clusters == 2])

cas <- sapply(strsplit(gsub("^X", "", colnames(Target)), "\\.\\."),
              function(x) { gsub("\\.", "-", x[1]) })
names(cas) <- cas

### 3. Count significant drug targets per cluster ########################################
# cas from colnames(Target) -> subset columns (PFAS)
c1_Target <- Target[, cas %in% c1]
c2_Target <- Target[, cas %in% c2]

cat("cluster 1 size:", length(c1), "/ cluster 2 size:", length(c2), "\n")

# apply(..., 1, ...) = count significant PFAS per row (per target category)
c1_Target_sig_raw <- apply(c1_Target, 1, function(x) { sum(na.omit(x) < 0.05) })
c2_Target_sig_raw <- apply(c2_Target, 1, function(x) { sum(na.omit(x) < 0.05) })

# Normalize by cluster size (proportion: %)
c1_Target_sig <- c1_Target_sig_raw / length(c1) * 100
c2_Target_sig <- c2_Target_sig_raw / length(c2) * 100

Target_sig    <- cbind(c1_Target_sig, c2_Target_sig)
x             <- Target_sig[which(rowSums(Target_sig) != 0), ]
df_Target_sig <- melt(x[order(x[, 2]), ])
df_Target_sig$Label <- str_wrap(df_Target_sig[, 1], width = 30)
df_Target_sig$group <- gsub("c1_Target_sig", "cluster 1", df_Target_sig$Var2)
df_Target_sig$group <- gsub("c2_Target_sig", "cluster 2", df_Target_sig$group)

colnames(df_Target_sig) <- c("Target", "cluster", "count","Label", "group")

### 4. Create bar plot #########################################################
g_Target <-
  ggplot(df_Target_sig, aes(x = Target, y = count, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(breaks = function(x) seq(0, ceiling(max(x)), by = 10),
                     minor_breaks = NULL,
                     limits = c(0, NA)) +
  labs(x = "Drug target", y = "Percentage of PFAS with significant enrichment (%)", fill = "Cluster") +
  theme_bw()

# Horizontal (rotate x-axis labels)
g_Target_h <- g_Target +
  theme(axis.text.x         = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position      = c(0.95, 0.95),
        legend.justification = c("right", "top"))

# Vertical (coord_flip)
g_Target_v <- g_Target +
  coord_flip() +
  scale_x_discrete(labels = df_Target_sig$Label) +
  theme(legend.position      = c(0.95, 0.05),
        legend.justification = c("right", "bottom"))

### 5. save #######################################################################
ggsave("Figure2c_enrichment_target.pdf", plot = g_Target_v, width = 8, height = 6, units = "in", device = "pdf")
