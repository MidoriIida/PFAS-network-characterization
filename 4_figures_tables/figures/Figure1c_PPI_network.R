############################
#        Figure 1A         #
#   PPI Subnetwork of      #
#   PFOA & Drug modules    #
############################
#
# Purpose:
#   Visualize target genes of PFOA (or PFOS) and proximal/distal drugs
#   on the PPI network (cf. Cheng et al. 2019 Nat Commun Fig.1a)
#
# Required packages:
#   install.packages(c("igraph","ggraph","tidygraph","dplyr","readxl","stringr","ggrepel"))
#
# Output: Figure1A.pdf (12 x 10 inch)
############################

library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(stringr)
library(ggrepel)

# ============================================================
# 1. Load data
# ============================================================

## --- Path settings ---
in_f_ppi      <- "./path/to/intermediate/data/g_comp.dat"
in_f_pfas     <- "./path/to/intermediate/data/list_PFAS_genes.dat"
in_f_zscore   <- "../../data/Supplementary_Data_5_Zdc_PFAS_Drug.txt"
in_f_ctd      <- "./path/to/intermediate/data/df_ctd_drugs.dat"

## --- 1-1. PPI network (igraph object) ---
load(in_f_ppi)   # loads g_comp
g_PPI <- g_comp; rm(g_comp)
cat("PPI nodes:", vcount(g_PPI), "/ edges:", ecount(g_PPI), "\n")

## --- 1-2. PFAS-associated gene list ---
load(in_f_pfas)   # loads list_PFAS_genes
# list_PFAS_genes$Merged: list of gene IDs keyed by CAS number

## --- 1-3. PFAS-drug z-score data ---
# Supplementary_Data_5 columns: ID Name dc_real mean(dc_null) sd(dc_null) Zdc(dc_null)
zscore_data <- read.table(in_f_zscore, sep = "\t", header = TRUE,
                          quote = "", stringsAsFactors = FALSE)
colnames(zscore_data)[c(5, 6)] <- c("sd_dc_null", "Zdc")

## --- 1-4. CTD drug-gene interaction data (900 approved drugs) ---
# produced by 0_data_preparation/3.drug_data_load.R
load(in_f_ctd)   # loads df_ctd_drugs

# Keep only genes present in the PPI network
ctd_human <- df_ctd_drugs %>%
  filter(GeneID %in% as.integer(names(V(g_PPI))))

# ============================================================
# 2. Select target PFAS and drugs
# ============================================================

## --- 2-0. Build GeneID -> GeneSymbol lookup table from CTD ---
id2symbol <- ctd_human %>%
  select(GeneID, GeneSymbol) %>%
  distinct() %>%
  mutate(GeneID = as.character(GeneID))
# Function: convert ID vector to symbols (returns original ID when no match)
to_symbol <- function(ids) {
  syms <- id2symbol$GeneSymbol[match(as.character(ids), id2symbol$GeneID)]
  ifelse(is.na(syms), as.character(ids), syms)
}

## --- 2-1. PFAS-associated genes ---
# CAS numbers: PFOA = 335-67-1 / PFOS = 1763-23-1
pfas_cas   <- "335-67-1"   # PFOA
pfas_name  <- "PFOA"
pfoa_genes <- unique(list_PFAS_genes$Merged[list_PFAS_genes$Merged$CASRN %in% pfas_cas, "GeneID"])
pfoa_genes <- pfoa_genes[pfoa_genes %in% names(V(g_PPI))]
cat(pfas_name, "genes in PPI:", length(pfoa_genes), "\n")

## --- 2-2. Pre-compute per-drug CTD gene counts for filtering ---
#
# Issue: some drugs in CTD have many non-specific interactions (e.g., Ivermectin > 4000).
# Solution: apply both lower and upper bounds on gene count,
#   keeping drugs with 20-200 genes for visualization.

get_drug_genes <- function(drug_name, ctd_df) {
  ctd_df %>%
    filter(ChemicalName == drug_name) %>%
    pull(GeneID) %>%
    unique()
}

## Get PFOA z-score data
pfoa_zscore <- zscore_data %>%
  filter(str_detect(ID, fixed(pfas_cas))) %>%
  arrange(Zdc)

## Compute gene count for all drugs (may take a moment)
cat("Counting CTD gene counts for all drugs in PFOA z-score table...\n")
pfoa_zscore$n_genes <- sapply(pfoa_zscore$Name, function(nm) {
  length(get_drug_genes(nm, ctd_human))
})

## Filter to drugs with 20-200 associated genes
##   lower bound 20: minimum meaningful module size (same as original analysis)
##   upper bound 200: exclude non-specific multi-binding drugs; ensures comparable sizes
pfoa_zscore_filtered <- pfoa_zscore %>%
  filter(n_genes >= 20, n_genes <= 600)

cat("Drugs with 20-200 genes:", nrow(pfoa_zscore_filtered), "out of", nrow(pfoa_zscore), "\n")

# Most proximal drug (smallest Zdc)
drug_close_row  <- pfoa_zscore_filtered[1, ]
drug_close_name <- drug_close_row$Name

# Most distal drug (largest Zdc)
drug_far_row  <- pfoa_zscore_filtered[nrow(pfoa_zscore_filtered), ]
drug_far_name <- drug_far_row$Name

cat("Closest drug :", drug_close_name,
    "(Zdc =", round(drug_close_row$Zdc, 2),
    ", n_genes =", drug_close_row$n_genes, ")\n")
cat("Farthest drug:", drug_far_name,
    "(Zdc =", round(drug_far_row$Zdc, 2),
    ", n_genes =", drug_far_row$n_genes, ")\n")

## --- 2-3. Retrieve target genes for each drug from CTD ---
drug_close_genes <- get_drug_genes(drug_close_name, ctd_human)
drug_far_genes   <- get_drug_genes(drug_far_name,   ctd_human)

cat("Genes for", drug_close_name, ":", length(drug_close_genes), "\n")
cat("Genes for", drug_far_name,   ":", length(drug_far_genes),   "\n")

# ============================================================
# 3. Extract subnetwork
# ============================================================

## --- 3-0. Unify data types and verify presence in PPI ---
# CTD GeneIDs may be integers; convert to character to match igraph vertex names
ppi_node_names   <- names(V(g_PPI))
pfoa_genes       <- as.character(pfoa_genes)
drug_close_genes <- as.character(drug_close_genes)
drug_far_genes   <- as.character(drug_far_genes)
pfoa_genes       <- pfoa_genes[pfoa_genes %in% ppi_node_names]
drug_close_genes <- drug_close_genes[drug_close_genes %in% ppi_node_names]
drug_far_genes   <- drug_far_genes[drug_far_genes %in% ppi_node_names]

cat("PFOA genes in PPI:", length(pfoa_genes), "\n")
cat(drug_close_name, "genes in PPI:", length(drug_close_genes), "\n")
cat(drug_far_name,   "genes in PPI:", length(drug_far_genes),   "\n")

## --- 3-1. Restrict each module to top-N genes by PPI degree ---
# Too many module genes cause the 1-hop neighborhood to span the entire PPI,
# so we restrict to high-centrality proteins (top PPI degree) as representatives.
MAX_MODULE <- 30   # max genes per module (adjust as needed)

top_by_degree <- function(genes, graph, n) {
  degs  <- degree(graph, v = genes, mode = "all")
  genes[order(degs, decreasing = TRUE)][1:min(n, length(genes))]
}

pfoa_genes_top       <- top_by_degree(pfoa_genes,       g_PPI, MAX_MODULE)
drug_close_genes_top <- top_by_degree(drug_close_genes, g_PPI, MAX_MODULE)
drug_far_genes_top   <- top_by_degree(drug_far_genes,   g_PPI, MAX_MODULE)

all_module_top <- unique(c(pfoa_genes_top, drug_close_genes_top, drug_far_genes_top))
cat("Module genes after top-N filter:", length(all_module_top), "\n")

## --- 3-2. Background nodes: 1-hop neighbors of PFOA ---
#
# Design rationale:
#   Sampling background nodes from PFOA neighborhood
#   forms a "PFOA cloud" on the PPI layout.
#   - Proximal drug (z-score < 0): shares many neighbors with PFOA
#     -> naturally pulled into the PFOA cloud in the layout
#   - Distal drug (z-score > 0): few neighbors shared with PFOA
#     -> pushed to the periphery of the PFOA cloud in the layout
#   This makes proximal/distal relationships visually apparent in the layout.

MAX_BG_NODES <- 60   # max background nodes (adjustable)

# Get 1-hop neighbors of PFOA module genes
pfoa_vids     <- V(g_PPI)[ppi_node_names %in% pfoa_genes_top]
pfoa_1hop     <- names(unique(unlist(
  ego(g_PPI, order = 1, nodes = pfoa_vids, mode = "all")
)))

# Exclude drug module genes (shown in a different color later)
drug_all_top  <- unique(c(drug_close_genes_top, drug_far_genes_top))
bg_candidates <- setdiff(pfoa_1hop, c(pfoa_genes_top, drug_all_top))

# Select in order of high PPI degree (prioritize hub proteins)
bg_deg   <- degree(g_PPI, v = bg_candidates, mode = "all")
bg_nodes <- names(sort(bg_deg, decreasing = TRUE))[
  seq_len(min(MAX_BG_NODES, length(bg_candidates)))
]
cat("Background nodes (PFOA 1-hop neighborhood):", length(bg_nodes), "\n")

## --- 3-3. Create induced subgraph ---
keep_genes <- unique(c(all_module_top, bg_nodes))
g_sub <- induced_subgraph(g_PPI, vids = V(g_PPI)[ppi_node_names %in% keep_genes])
cat("Subnetwork nodes:", vcount(g_sub), "/ edges:", ecount(g_sub), "\n")

## --- 3-4. Extract largest connected component ---
comp        <- components(g_sub, mode = "weak")
largest_id  <- which.max(comp$csize)
g_sub       <- induced_subgraph(g_sub,
                 vids = V(g_sub)[comp$membership == largest_id])
cat("Largest component nodes:", vcount(g_sub),
    "/ edges:", ecount(g_sub), "\n")

# ============================================================
# 4. Assign node attributes
# ============================================================

node_df <- data.frame(
  gene = names(V(g_sub)),   # GeneID (character)
  stringsAsFactors = FALSE
) %>%
  mutate(
    symbol   = to_symbol(gene),   # GeneID → GeneSymbol
    is_pfas  = gene %in% pfoa_genes_top,
    is_close = gene %in% drug_close_genes_top,
    is_far   = gene %in% drug_far_genes_top,
    role = case_when(
      is_pfas & is_close        ~ paste0("Shared (", pfas_name, " & close drug)"),
      is_pfas & is_far          ~ paste0("Shared (", pfas_name, " & far drug)"),
      is_pfas                   ~ paste0(pfas_name, " targets"),
      is_close                  ~ drug_close_name,
      is_far                    ~ drug_far_name,
      TRUE                      ~ "Background"
    ),
    # Show gene symbol label for module nodes only; NA for background nodes
    label = if_else(role != "Background", symbol, NA_character_)
  )

## Set attributes on igraph object
V(g_sub)$role   <- node_df$role[match(names(V(g_sub)), node_df$gene)]
V(g_sub)$label  <- node_df$label[match(names(V(g_sub)), node_df$gene)]
V(g_sub)$symbol <- node_df$symbol[match(names(V(g_sub)), node_df$gene)]

# ============================================================
# 5. Visualization（ggraph）
# ============================================================

## Color palette (matching reference figure)
pfas_targets_label <- paste0(pfas_name, " targets")
shared_close_label <- paste0("Shared (", pfas_name, " & close drug)")
shared_far_label   <- paste0("Shared (", pfas_name, " & far drug)")

role_levels <- c(
  pfas_targets_label,
  drug_close_name,
  shared_close_label,
  drug_far_name,
  shared_far_label,
  "Background"
)

role_colors <- c(
  "#E8476A",   # PFAS targets: pink-red
  "#56B4E9",   # close drug:   sky blue
  "#F5C842",   # shared close: yellow
  "#2ECC71" ,  # far drug:     green
  "#C36CC9",   # shared far:   purple
  "#C8C8C8"    # Background:   gray
)
names(role_colors) <- role_levels

role_size <- c(
  "PFOA targets"               = 4.5,
  "Shared (PFOA & close drug)" = 4.5,
  "Shared (PFOA & far drug)"   = 4.5,
  "Background"                 = 2.5
)
role_size[drug_close_name] <- 4.5
role_size[drug_far_name]   <- 4.5

## Assign edge attributes on igraph before converting to tidygraph
## (activate(edges) + .N() is unstable across versions; handle directly in igraph)

# Get vertex roles as a named vector
node_roles_vec <- setNames(as.character(V(g_sub)$role), names(V(g_sub)))

# Look up the role at each edge endpoint
el         <- as_edgelist(g_sub, names = TRUE)
from_roles <- node_roles_vec[el[, 1]]
to_roles   <- node_roles_vec[el[, 2]]

# Classify edge types
E(g_sub)$edge_type <- case_when(
  # Edges within the same module
  from_roles == to_roles & from_roles == pfas_targets_label ~ pfas_targets_label,
  from_roles == to_roles & from_roles == drug_close_name    ~ drug_close_name,
  from_roles == to_roles & from_roles == drug_far_name      ~ drug_far_name,
  # # Cross-edges: PFOA target <-> Ketoprofen (blue)
  # (from_roles == pfas_targets_label & to_roles == drug_close_name) |
  # (from_roles == drug_close_name    & to_roles == pfas_targets_label) ~ drug_close_name,
  # # Cross-edges: PFOA target <-> Norfloxacin (green)
  # (from_roles == pfas_targets_label & to_roles == drug_far_name) |
  # (from_roles == drug_far_name      & to_roles == pfas_targets_label) ~ drug_far_name,
  TRUE                                                       ~ "other"
)

## Convert to tidygraph (edge_type carried over as igraph attribute)
tg <- as_tbl_graph(g_sub) %>%
  activate(nodes) %>%
  mutate(role  = factor(role, levels = role_levels),
         label = label)

## Edge color/width scale (unified with node colors)
edge_colors <- c(role_colors[pfas_targets_label],
                 role_colors[drug_close_name],
                 role_colors[drug_far_name])
names(edge_colors) <- c(pfas_targets_label, drug_close_name, drug_far_name)

edge_widths <- c(
  setNames(1.5, pfas_targets_label),
  setNames(1.5, drug_close_name),
  setNames(1.5, drug_far_name)
)

## Edge color/width scale (including "other"; defined in one block)
all_edge_colors <- c(edge_colors, "other" = "gray80")
all_edge_widths <- c(edge_widths, "other" = 0.3)
all_edge_alphas <- c(
  setNames(0.75, pfas_targets_label),
  setNames(0.75, drug_close_name),
  setNames(0.75, drug_far_name),
  "other" = 0.35
)

set.seed(42)
p <- ggraph(tg, layout = "kk") +   # Kamada-Kawai: nodes with more shared neighbors placed closer
  # Draw all edges in one layer (color, width, alpha controlled by scale)
  geom_edge_link(aes(colour = edge_type,
                     width  = edge_type,
                     alpha  = edge_type)) +
  scale_edge_colour_manual(values = all_edge_colors, guide = "none") +
  scale_edge_width_manual( values = all_edge_widths, guide = "none") +
  scale_edge_alpha_manual( values = all_edge_alphas, guide = "none") +
  # Nodes (color and size vary by role)
  geom_node_point(aes(color = role, size = role), show.legend = TRUE) +
  # Labels (module nodes only)
  geom_node_label(
    aes(label = label, color = role),
    size = 12, repel = TRUE,
    label.padding = unit(0.15, "lines"),
    label.size    = 0.2,
    show.legend   = FALSE,
    na.rm         = TRUE
  ) +
  # Scales
  scale_color_manual(values = role_colors, na.value = "#C8C8C8",
                     breaks = role_levels[role_levels != "Background"],
                     name = NULL) +
  scale_size_manual(values = role_size, guide = "none") +
  # Theme
  theme_graph(base_family = "Helvetica") +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    legend.key       = element_blank(),
    legend.box       = "horizontal"
  ) +
  guides(color = guide_legend(
    nrow = 2,
    override.aes = list(size = 5)
  ))

## Output
out_dir  <- "/path/to/your/local/directory/Code_for_Manuscript"
out_path <- file.path(out_dir, "Figure1A.pdf")
ggsave(out_path, plot = p, width = 21, height = 10, units = "in")
cat("Saved:", out_path, "\n")

# ============================================================
# 6. Cytoscape export
# ============================================================
# Two export options; choose whichever is more convenient.

## -------------------------------------------------------
## Option A: GraphML format (single file; recommended)
##   In Cytoscape: File > Import > Network from File > Figure1A.graphml
## -------------------------------------------------------

# Assign node attributes to igraph object (symbol, role, color)
V(g_sub)$color_hex <- role_colors[V(g_sub)$role]   # used for Cytoscape color settings

graphml_path <- file.path(out_dir, "Figure1A.graphml")
write_graph(g_sub, graphml_path, format = "graphml")
cat("GraphML saved:", graphml_path, "\n")

## -------------------------------------------------------
## Option B: CSV format (node table + edge table)
##   Cytoscape procedure:
##   1. File > Import > Network from File → Figure1A_edges.csv
##      (set source column = source node, target column = target node)
##   2. File > Import > Table from File  → Figure1A_nodes.csv
##      (merge into existing nodes using "id" column as key)
## -------------------------------------------------------

# Edge table (using gene symbols)
el_sub     <- as_edgelist(g_sub, names = TRUE)
edge_csv   <- data.frame(
  source    = to_symbol(el_sub[, 1]),
  target    = to_symbol(el_sub[, 2]),
  edge_type = E(g_sub)$edge_type,
  color_hex = all_edge_colors[E(g_sub)$edge_type],
  stringsAsFactors = FALSE
)

# Node table
node_csv <- data.frame(
  id        = names(V(g_sub)),
  symbol    = V(g_sub)$symbol,
  role      = V(g_sub)$role,
  color_hex = role_colors[V(g_sub)$role],
  stringsAsFactors = FALSE
)

write.csv(edge_csv, file.path(out_dir, "Figure1A_edges.csv"), row.names = FALSE, quote = FALSE)
write.csv(node_csv, file.path(out_dir, "Figure1A_nodes.csv"), row.names = FALSE, quote = FALSE)
cat("CSV saved: Figure1A_edges.csv / Figure1A_nodes.csv\n")

## -------------------------------------------------------
## Option C: RCy3 (control Cytoscape directly from R; requires Cytoscape running)
##   install.packages("BiocManager"); BiocManager::install("RCy3")
## -------------------------------------------------------
# library(RCy3)
# cytoscapePing()   # check if Cytoscape is running
#
# createNetworkFromIgraph(
#   igraph     = g_sub,
#   title      = "Figure1A_PFOA_network",
#   collection = "PFAS manuscript"
# )
#
# # Set node color by role
# setNodeColorMapping(
#   table.column      = "role",
#   table.column.values = names(role_colors),
#   colors            = unname(role_colors),
#   mapping.type      = "discrete"
# )
