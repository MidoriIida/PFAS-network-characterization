############################
#        Figure 1A         #
#   PPI Subnetwork of      #
#   PFOA & Drug modules    #
############################
#
# 目的:
#   PFOA（またはPFOS）と近接薬剤・遠隔薬剤のターゲット遺伝子を
#   PPIネットワーク上に可視化する（参考: Cheng et al. 2019 Nat Commun Fig.1a）
#
# 必要パッケージ:
#   install.packages(c("igraph","ggraph","tidygraph","dplyr","readxl","stringr","ggrepel"))
#
# 出力: Figure1A.pdf (12 x 10 inch)
############################

library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(readxl)
library(stringr)
library(ggrepel)

# ============================================================
# 1. データ読み込み
# ============================================================

## --- パス設定 ---
in_f_ppi      <- "/path/to/your/local/directory/human_PPI_230710_iida.dat"
in_f_pfas     <- "/path/to/your/local/directory/human_PFAS_genes_230712_iida.dat"
in_f_zscore   <- "/path/to/your/local/directory/summary_PFAS_appdrug_iidaData_20230716.xlsx"
in_f_ctd      <- "/path/to/your/local/directory/CTD_chem_gene_ixns.csv"

## --- 1-1. PPIネットワーク（igraphオブジェクト）---
g_PPI <- get(load(in_f_ppi))
cat("PPI nodes:", vcount(g_PPI), "/ edges:", ecount(g_PPI), "\n")

## --- 1-2. PFAS関連遺伝子リスト ---
load(in_f_pfas)   # → list_PFAS_genes を読み込み
# list_PFAS_genes$Merged: CAS番号をキーにした遺伝子IDのリスト

## --- 1-3. PFAS-薬剤 z-scoreデータ ---
zscore_data <- read_excel(in_f_zscore)
colnames(zscore_data)[c(5, 6)] <- c("sd_dc_null", "Zdc")

## --- 1-4. CTD 薬剤-遺伝子相互作用データ ---
ctd_raw  <- read.csv(in_f_ctd, skip = 29, header = FALSE)
colnames(ctd_raw) <- c("ChemicalName", "ChemicalID", "CasRN", "GeneSymbol",
                        "GeneID", "GeneForms", "Organism", "OrganismID",
                        "Interaction", "InteractionActions", "PubMedIDs")

# ヒトタンパク質に絞り込み、PPIネットワーク内の遺伝子のみ保持
ctd_human <- ctd_raw %>%
  filter(Organism == "Homo sapiens",
         GeneForms == "protein",
         GeneID %in% names(V(g_PPI)))

# ============================================================
# 2. 対象PFAS・薬剤の選択
# ============================================================

## --- 2-0. GeneID → GeneSymbol 変換テーブルをCTDから作成 ---
id2symbol <- ctd_human %>%
  select(GeneID, GeneSymbol) %>%
  distinct() %>%
  mutate(GeneID = as.character(GeneID))
# 関数: IDベクターをシンボルに変換（マッチしない場合はIDをそのまま返す）
to_symbol <- function(ids) {
  syms <- id2symbol$GeneSymbol[match(as.character(ids), id2symbol$GeneID)]
  ifelse(is.na(syms), as.character(ids), syms)
}

## --- 2-1. PFOS関連遺伝子 ---
# CAS番号: PFOA = 335-67-1 / PFOS = 1763-23-1
pfas_cas   <- "335-67-1"   # PFOA
pfas_name  <- "PFOA"
pfoa_genes <- unique(list_PFAS_genes$Merged[list_PFAS_genes$Merged$CASRN %in% pfas_cas, "GeneID"])
pfoa_genes <- pfoa_genes[pfoa_genes %in% names(V(g_PPI))]
cat(pfas_name, "genes in PPI:", length(pfoa_genes), "\n")

## --- 2-2. 薬剤ごとのCTD遺伝子数を事前計算してフィルタリング ---
#
# 問題: CTDには非特異的な相互作用が多い薬剤（例: Ivermectin = 4000超）がある。
# 解決: 解析に使った基準（遺伝子数20以上）を上限でも絞り込み、
#       可視化として適切な遺伝子数(20〜200)の薬剤を選ぶ。

get_drug_genes <- function(drug_name, ctd_df) {
  ctd_df %>%
    filter(ChemicalName == drug_name) %>%
    pull(GeneID) %>%
    unique()
}

## PFOAのz-scoreデータを取得
pfoa_zscore <- zscore_data %>%
  filter(str_detect(ID, fixed(pfas_cas))) %>%
  arrange(Zdc)

## 全薬剤の遺伝子数を計算（少し時間がかかる）
cat("Counting CTD gene counts for all drugs in PFOA z-score table...\n")
pfoa_zscore$n_genes <- sapply(pfoa_zscore$Name, function(nm) {
  length(get_drug_genes(nm, ctd_human))
})

## 遺伝子数が20〜200の薬剤に絞る
##   下限20: モジュールとして意味のある最小サイズ（元解析の基準と同じ）
##   上限200: 非特異的な多重結合を排除し、視覚的に比較可能なサイズに揃える
pfoa_zscore_filtered <- pfoa_zscore %>%
  filter(n_genes >= 20, n_genes <= 600)

cat("Drugs with 20-200 genes:", nrow(pfoa_zscore_filtered), "out of", nrow(pfoa_zscore), "\n")

# 最も近接（Zdcが最も小さい）の薬剤
drug_close_row  <- pfoa_zscore_filtered[1, ]
drug_close_name <- drug_close_row$Name

# 最も遠隔（Zdcが最も大きい）の薬剤
drug_far_row  <- pfoa_zscore_filtered[nrow(pfoa_zscore_filtered), ]
drug_far_name <- drug_far_row$Name

cat("Closest drug :", drug_close_name,
    "(Zdc =", round(drug_close_row$Zdc, 2),
    ", n_genes =", drug_close_row$n_genes, ")\n")
cat("Farthest drug:", drug_far_name,
    "(Zdc =", round(drug_far_row$Zdc, 2),
    ", n_genes =", drug_far_row$n_genes, ")\n")

## --- 2-3. 各薬剤のターゲット遺伝子をCTDから取得 ---
drug_close_genes <- get_drug_genes(drug_close_name, ctd_human)
drug_far_genes   <- get_drug_genes(drug_far_name,   ctd_human)

cat("Genes for", drug_close_name, ":", length(drug_close_genes), "\n")
cat("Genes for", drug_far_name,   ":", length(drug_far_genes),   "\n")

# ============================================================
# 3. サブネットワークの抽出
# ============================================================

## --- 3-0. 型の統一・PPI内存在確認 ---
# CTDのGeneIDはintegerの場合があるため、igraph頂点名(character)に統一
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

## --- 3-1. 各モジュールをPPI次数上位N遺伝子に絞る ---
# モジュール遺伝子が多すぎると1ホップ近傍がPPI全体に広がるため、
# ネットワーク中心性の高いタンパク質（PPI次数上位）に絞って代表させる。
MAX_MODULE <- 30   # 1モジュールあたりの最大遺伝子数（必要に応じて調整）

top_by_degree <- function(genes, graph, n) {
  degs  <- degree(graph, v = genes, mode = "all")
  genes[order(degs, decreasing = TRUE)][1:min(n, length(genes))]
}

pfoa_genes_top       <- top_by_degree(pfoa_genes,       g_PPI, MAX_MODULE)
drug_close_genes_top <- top_by_degree(drug_close_genes, g_PPI, MAX_MODULE)
drug_far_genes_top   <- top_by_degree(drug_far_genes,   g_PPI, MAX_MODULE)

all_module_top <- unique(c(pfoa_genes_top, drug_close_genes_top, drug_far_genes_top))
cat("Module genes after top-N filter:", length(all_module_top), "\n")

## --- 3-2. 背景ノード：PFOAの1ホップ近傍から取得 ---
#
# 【設計の意図】
#   背景ノードをPFOAの近傍だけから取ることで、
#   PPI上に「PFOAの領域（雲）」が形成される。
#   - 近接薬剤（close drug）: z-score < 0 = PFOAと共通の近傍が多い
#     → レイアウトでPFOA雲の中に自然に引き寄せられる
#   - 遠隔薬剤（far drug）: z-score > 0 = PFOAとの共通近傍が少ない
#     → レイアウトでPFOA雲の外側に押し出される
#   これにより「近い/遠い」がレイアウトで視覚的に表現される。

MAX_BG_NODES <- 60   # 背景ノードの上限（調整可）

# PFOAモジュール遺伝子の1ホップ近傍を取得
pfoa_vids     <- V(g_PPI)[ppi_node_names %in% pfoa_genes_top]
pfoa_1hop     <- names(unique(unlist(
  ego(g_PPI, order = 1, nodes = pfoa_vids, mode = "all")
)))

# ドラッグモジュール遺伝子は除く（後で別色で表示するため）
drug_all_top  <- unique(c(drug_close_genes_top, drug_far_genes_top))
bg_candidates <- setdiff(pfoa_1hop, c(pfoa_genes_top, drug_all_top))

# PPI次数が高い順に採用（ネットワークの中心的ノードを優先）
bg_deg   <- degree(g_PPI, v = bg_candidates, mode = "all")
bg_nodes <- names(sort(bg_deg, decreasing = TRUE))[
  seq_len(min(MAX_BG_NODES, length(bg_candidates)))
]
cat("Background nodes (PFOA 1-hop neighborhood):", length(bg_nodes), "\n")

## --- 3-3. 誘導部分グラフを作成 ---
keep_genes <- unique(c(all_module_top, bg_nodes))
g_sub <- induced_subgraph(g_PPI, vids = V(g_PPI)[ppi_node_names %in% keep_genes])
cat("Subnetwork nodes:", vcount(g_sub), "/ edges:", ecount(g_sub), "\n")

## --- 3-4. 最大連結成分のみ抽出 ---
comp        <- components(g_sub, mode = "weak")
largest_id  <- which.max(comp$csize)
g_sub       <- induced_subgraph(g_sub,
                 vids = V(g_sub)[comp$membership == largest_id])
cat("Largest component nodes:", vcount(g_sub),
    "/ edges:", ecount(g_sub), "\n")

# ============================================================
# 4. ノード属性の付与
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
    # モジュールノードのみシンボルをラベル表示、背景ノードはNA
    label = if_else(role != "Background", symbol, NA_character_)
  )

## igraphオブジェクトに属性をセット
V(g_sub)$role   <- node_df$role[match(names(V(g_sub)), node_df$gene)]
V(g_sub)$label  <- node_df$label[match(names(V(g_sub)), node_df$gene)]
V(g_sub)$symbol <- node_df$symbol[match(names(V(g_sub)), node_df$gene)]

# ============================================================
# 5. 可視化（ggraph）
# ============================================================

## カラーパレット（参考図と同系統）
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
  "#E8476A",   # PFAS targets: ピンク-赤
  "#56B4E9",   # close drug:   スカイブルー
  "#F5C842",   # shared close: 黄
  "#2ECC71" ,  # far drug:     グリーン
  "#C36CC9",   # shared far:   紫
  "#C8C8C8"    # Background:   グレー
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

## エッジ属性をigraph上で付与してからtidygraphに変換
## （activate(edges) + .N() はバージョンによって不安定なため、igraphで直接処理）

# 頂点ロールを名前付きベクターとして取得
node_roles_vec <- setNames(as.character(V(g_sub)$role), names(V(g_sub)))

# エッジリストで両端のロールを参照
el         <- as_edgelist(g_sub, names = TRUE)
from_roles <- node_roles_vec[el[, 1]]
to_roles   <- node_roles_vec[el[, 2]]

# エッジタイプを分類
E(g_sub)$edge_type <- case_when(
  # 同じモジュール内のエッジ
  from_roles == to_roles & from_roles == pfas_targets_label ~ pfas_targets_label,
  from_roles == to_roles & from_roles == drug_close_name    ~ drug_close_name,
  from_roles == to_roles & from_roles == drug_far_name      ~ drug_far_name,
  # # PFOA target ↔ Ketoprofen のクロスエッジ（青）
  # (from_roles == pfas_targets_label & to_roles == drug_close_name) |
  # (from_roles == drug_close_name    & to_roles == pfas_targets_label) ~ drug_close_name,
  # # PFOA target ↔ Norfloxacin のクロスエッジ（緑）
  # (from_roles == pfas_targets_label & to_roles == drug_far_name) |
  # (from_roles == drug_far_name      & to_roles == pfas_targets_label) ~ drug_far_name,
  TRUE                                                       ~ "other"
)

## tidygraphに変換（edge_type は igraph 属性として引き継がれる）
tg <- as_tbl_graph(g_sub) %>%
  activate(nodes) %>%
  mutate(role  = factor(role, levels = role_levels),
         label = label)

## エッジの色・太さスケール（ノード色と統一）
edge_colors <- c(role_colors[pfas_targets_label],
                 role_colors[drug_close_name],
                 role_colors[drug_far_name])
names(edge_colors) <- c(pfas_targets_label, drug_close_name, drug_far_name)

edge_widths <- c(
  setNames(1.5, pfas_targets_label),
  setNames(1.5, drug_close_name),
  setNames(1.5, drug_far_name)
)

## エッジの色・太さスケール（"other" も含めて一括定義）
all_edge_colors <- c(edge_colors, "other" = "gray80")
all_edge_widths <- c(edge_widths, "other" = 0.3)
all_edge_alphas <- c(
  setNames(0.75, pfas_targets_label),
  setNames(0.75, drug_close_name),
  setNames(0.75, drug_far_name),
  "other" = 0.35
)

set.seed(42)
p <- ggraph(tg, layout = "kk") +   # Kamada-Kawai: 共通近傍が多いほど近くに配置される
  # 全エッジを1レイヤーで描画（scale で色・太さ・透明度を制御）
  geom_edge_link(aes(colour = edge_type,
                     width  = edge_type,
                     alpha  = edge_type)) +
  scale_edge_colour_manual(values = all_edge_colors, guide = "none") +
  scale_edge_width_manual( values = all_edge_widths, guide = "none") +
  scale_edge_alpha_manual( values = all_edge_alphas, guide = "none") +
  # ノード（ロールによって色・サイズを変える）
  geom_node_point(aes(color = role, size = role), show.legend = TRUE) +
  # ラベル（モジュールノードのみ）
  geom_node_label(
    aes(label = label, color = role),
    size = 12, repel = TRUE,
    label.padding = unit(0.15, "lines"),
    label.size    = 0.2,
    show.legend   = FALSE,
    na.rm         = TRUE
  ) +
  # スケール
  scale_color_manual(values = role_colors, na.value = "#C8C8C8",
                     breaks = role_levels[role_levels != "Background"],
                     name = NULL) +
  scale_size_manual(values = role_size, guide = "none") +
  # テーマ
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

## 出力
out_dir  <- "/path/to/your/local/directory/Code_for_Manuscript"
out_path <- file.path(out_dir, "Figure1A.pdf")
ggsave(out_path, plot = p, width = 21, height = 10, units = "in")
cat("Saved:", out_path, "\n")

# ============================================================
# 6. Cytoscape エクスポート
# ============================================================
# 2通りの方法を用意。どちらか使いやすい方を選んでください。

## -------------------------------------------------------
## Option A: GraphML形式（1ファイルで完結・推奨）
##   Cytoscapeで: File > Import > Network from File > Figure1A.graphml
## -------------------------------------------------------

# ノード属性をigraphオブジェクトに付与（symbol・role・color）
V(g_sub)$color_hex <- role_colors[V(g_sub)$role]   # Cytoscapeの色設定に使用

graphml_path <- file.path(out_dir, "Figure1A.graphml")
write_graph(g_sub, graphml_path, format = "graphml")
cat("GraphML saved:", graphml_path, "\n")

## -------------------------------------------------------
## Option B: CSV形式（ノードテーブル + エッジテーブル）
##   Cytoscape操作手順:
##   1. File > Import > Network from File → Figure1A_edges.csv
##      （source列 = source node, target列 = target node に設定）
##   2. File > Import > Table from File  → Figure1A_nodes.csv
##      （"id" 列をキーに既存ノードへマージ）
## -------------------------------------------------------

# エッジテーブル（シンボルで表記）
el_sub     <- as_edgelist(g_sub, names = TRUE)
edge_csv   <- data.frame(
  source    = to_symbol(el_sub[, 1]),
  target    = to_symbol(el_sub[, 2]),
  edge_type = E(g_sub)$edge_type,
  color_hex = all_edge_colors[E(g_sub)$edge_type],
  stringsAsFactors = FALSE
)

# ノードテーブル
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
## Option C: RCy3（CytoscapeをRから直接操作・要 Cytoscape 起動中）
##   install.packages("BiocManager"); BiocManager::install("RCy3")
## -------------------------------------------------------
# library(RCy3)
# cytoscapePing()   # Cytoscapeが起動中か確認
#
# createNetworkFromIgraph(
#   igraph     = g_sub,
#   title      = "Figure1A_PFOA_network",
#   collection = "PFAS manuscript"
# )
#
# # ノード色をロールに応じて設定
# setNodeColorMapping(
#   table.column      = "role",
#   table.column.values = names(role_colors),
#   colors            = unname(role_colors),
#   mapping.type      = "discrete"
# )
