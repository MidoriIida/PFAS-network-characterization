####################################################
#  Table 1 + ATC + ChEMBL MOA + KEGG br08310      #
####################################################

library(readxl)
library(openxlsx)
library(dplyr)
library(stringr)
library(httr)
library(jsonlite)

# ============================================================
# 1. データ読み込み
# ============================================================

## Table 1
table1 <- read_excel(
  "/path/to/your/local/directory/Table1_PFAS_drug_proximity.xlsx"
)
colnames(table1) <- c("ID", "Name", "dc_real", "Zdc")

## ATC分類データ（複数コードは";"結合）
atc_raw <- read.table(
  "/path/to/your/local/directory/ATC_code.txt",
  sep = "\t", header = FALSE, quote = "", fill = TRUE,
  stringsAsFactors = FALSE
)
colnames(atc_raw) <- c("drug_name",
                        "atc_l5","atc_l4","atc_l4_name",
                        "atc_l3","atc_l3_name",
                        "atc_l2","atc_l2_name",
                        "atc_l1","atc_l1_name")
atc_lookup <- atc_raw %>%
  mutate(drug_name_lower = tolower(drug_name)) %>%
  group_by(drug_name_lower) %>%
  summarise(
    atc_l1      = paste(unique(atc_l1),      collapse = "; "),
    atc_l1_name = paste(unique(atc_l1_name), collapse = "; "),
    atc_l2      = paste(unique(atc_l2),      collapse = "; "),
    atc_l2_name = paste(unique(atc_l2_name), collapse = "; "),
    .groups = "drop"
  )

# ============================================================
# 2. PFAS名の付与・マッチング用キー列の作成
# ============================================================

pfas_map <- c("335-67-1" = "PFOA", "1763-23-1" = "PFOS")

table1 <- table1 %>%
  mutate(
    cas        = str_extract(ID, "^[^|]+"),
    PFAS       = if_else(is.na(pfas_map[cas]), cas, pfas_map[cas]),
    name_lower = tolower(Name),
    name_short = str_trim(str_remove(
      name_lower,
      "\\s+(hydrochloride|sulfate|acetate|sodium|chloride|mesylate|maleate|tartrate)$"
    ))
  )

# ============================================================
# 3. ATC結合（3段階マッチング）
# ============================================================

joined <- left_join(table1,
                    atc_lookup %>% select(drug_name_lower,
                                          atc_l1, atc_l1_name,
                                          atc_l2, atc_l2_name),
                    by = c("name_lower" = "drug_name_lower"))

# Step2: 塩形式除去名で再試行
na_idx <- is.na(joined$atc_l1)
if (any(na_idx)) {
  joined[na_idx, c("atc_l1","atc_l1_name","atc_l2","atc_l2_name")] <-
    left_join(table1[na_idx,] %>% select(name_short),
              atc_lookup %>% select(drug_name_lower,
                                    atc_l1, atc_l1_name,
                                    atc_l2, atc_l2_name),
              by = c("name_short" = "drug_name_lower"))[,
                c("atc_l1","atc_l1_name","atc_l2","atc_l2_name")]
}

# Step3: 手動補完
manual_atc <- data.frame(
  name_lower   = c("omacor",   "cortisone", "opicapone",  "methylergonovine"),
  atc_l1       = c("C",        "H",          "N",          "G"),
  atc_l1_name  = c("CARDIOVASCULAR SYSTEM",
                   "SYSTEMIC HORMONAL PREPARATIONS, EXCL. SEX HORMONES",
                   "NERVOUS SYSTEM",
                   "GENITO URINARY SYSTEM AND SEX HORMONES"),
  atc_l2       = c("C10",      "H02",        "N04",        "G02"),
  atc_l2_name  = c("LIPID MODIFYING AGENTS",
                   "CORTICOSTEROIDS FOR SYSTEMIC USE",
                   "ANTI-PARKINSON DRUGS",
                   "OTHER GYNECOLOGICALS"),
  stringsAsFactors = FALSE
)
na_idx2 <- is.na(joined$atc_l1)
if (any(na_idx2)) {
  joined[na_idx2, c("atc_l1","atc_l1_name","atc_l2","atc_l2_name")] <-
    left_join(table1[na_idx2,] %>% select(name_lower),
              manual_atc, by = "name_lower")[,
                c("atc_l1","atc_l1_name","atc_l2","atc_l2_name")]
}

cat("=== ATC マッチ確認 ===\n")
joined %>% select(Name, atc_l1, atc_l2, atc_l2_name) %>%
  distinct() %>% print(n = 30)

# ============================================================
# 4. ChEMBL MOA 取得
# ============================================================

get_chembl_id <- function(drug_name) {
  url  <- paste0("https://www.ebi.ac.uk/chembl/api/data/molecule?",
                 "pref_name=", URLencode(toupper(drug_name)),
                 "&format=json&limit=1")
  resp <- tryCatch(GET(url, timeout(15)), error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) return(NA_character_)
  ct   <- fromJSON(content(resp, "text", encoding = "UTF-8"))
  ids  <- ct$molecules$molecule_chembl_id
  if (length(ids) == 0 || is.null(ids)) NA_character_ else ids[1]
}

get_chembl_moa <- function(chembl_id) {
  if (is.na(chembl_id)) return(NA_character_)
  url  <- paste0("https://www.ebi.ac.uk/chembl/api/data/mechanism?",
                 "molecule_chembl_id=", chembl_id, "&format=json")
  resp <- tryCatch(GET(url, timeout(15)), error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) return(NA_character_)
  ct   <- fromJSON(content(resp, "text", encoding = "UTF-8"))
  mechs <- unique(ct$mechanisms$mechanism_of_action)
  mechs <- mechs[!is.na(mechs)]
  if (length(mechs) == 0) NA_character_ else paste(mechs, collapse = "; ")
}

unique_drugs <- unique(joined$Name)
cat("\n=== ChEMBL MOA 取得中 ===\n")

chembl_df <- data.frame(Name = unique_drugs,
                         chembl_id = NA_character_,
                         chembl_moa = NA_character_,
                         stringsAsFactors = FALSE)

for (i in seq_along(unique_drugs)) {
  drug <- unique_drugs[i]
  cat(sprintf("  [%d/%d] %s", i, length(unique_drugs), drug))

  cid <- get_chembl_id(drug)
  if (is.na(cid)) {
    cid <- get_chembl_id(table1$name_short[match(tolower(drug), table1$name_lower)])
    if (!is.na(cid)) cat(" (short)")
  }
  chembl_df$chembl_id[i]  <- cid
  chembl_df$chembl_moa[i] <- get_chembl_moa(cid)
  cat(" →", ifelse(is.na(chembl_df$chembl_moa[i]), "not found",
                   chembl_df$chembl_moa[i]), "\n")
  Sys.sleep(0.3)
}

# ============================================================
# 5. KEGG Target-based Classification (br08310) 取得・解析
# ============================================================

## --- 5-1. br08310 ダウンロード ---
cat("\n=== KEGG br08310 取得中 ===\n")
resp_brite <- GET("https://rest.kegg.jp/get/br:br08310", timeout(60))
if (status_code(resp_brite) != 200) stop("KEGG API失敗: ", status_code(resp_brite))
brite_lines <- strsplit(content(resp_brite, "text", encoding = "UTF-8"), "\n")[[1]]
cat("取得:", length(brite_lines), "行\n")

## --- 5-2. パース: Drug ID → 分類パス ---
parse_brite <- function(lines) {
  cur <- list(A = "", B = "", C = "", D = "")
  out <- list()
  for (line in lines) {
    if (grepl("^E\\s+D\\d{5}", line)) {
      kid  <- str_extract(line, "D\\d{5}")
      kname <- str_trim(str_remove(line, "^E\\s+D\\d{5}\\s*"))
      kname <- str_trim(str_remove(kname, "\\s*\\(.*"))
      out[[length(out)+1]] <- data.frame(
        kegg_id   = kid,
        kegg_name = kname,
        cat_A = str_trim(cur$A),
        cat_B = str_trim(cur$B),
        cat_C = str_trim(cur$C),
        cat_D = str_trim(cur$D),
        stringsAsFactors = FALSE
      )
    } else if (grepl("^A\\s", line)) {
      cur$A <- str_trim(str_remove(line, "^A\\s+"))
      cur$B <- ""; cur$C <- ""; cur$D <- ""
    } else if (grepl("^B\\s", line)) {
      cur$B <- str_trim(str_remove(line, "^B\\s+"))
      cur$C <- ""; cur$D <- ""
    } else if (grepl("^C\\s", line)) {
      cur$C <- str_trim(str_remove(line, "^C\\s+"))
      cur$D <- ""
    } else if (grepl("^D\\s", line)) {
      cur$D <- str_trim(str_remove(line, "^D\\s+"))
      cur$D <- str_trim(str_remove(cur$D, "^DG\\d+\\s*"))
    }
  }
  bind_rows(out)
}

brite_df <- parse_brite(brite_lines)
cat("解析:", nrow(brite_df), "Drug Entry\n")

save(brite_df, file="./brite_df.dat")

brite_lookup <- brite_df %>%
  group_by(kegg_id) %>%
  summarise(
    kegg_cat_A = paste(unique(cat_A[cat_A != ""]), collapse = "; "),
    kegg_cat_B = paste(unique(cat_B[cat_B != ""]), collapse = "; "),
    kegg_cat_C = paste(unique(cat_C[cat_C != ""]), collapse = "; "),
    kegg_cat_D = paste(unique(cat_D[cat_D != ""]), collapse = "; "),
    .groups = "drop"
  )

## --- 5-3. 薬剤名 → KEGG Drug ID ---
get_kegg_id <- function(drug_name) {
  url  <- paste0("https://rest.kegg.jp/find/drug/", URLencode(drug_name))
  resp <- tryCatch(GET(url, timeout(15)), error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) return(NA_character_)
  txt  <- content(resp, "text", encoding = "UTF-8")
  if (nchar(trimws(txt)) == 0) return(NA_character_)
  str_extract(strsplit(txt, "\n")[[1]][1], "D\\d{5}")
}

cat("\n=== KEGG Drug ID 取得中 ===\n")
kegg_id_df <- data.frame(Name = unique_drugs,
                          kegg_id = NA_character_,
                          stringsAsFactors = FALSE)

for (i in seq_along(unique_drugs)) {
  drug <- unique_drugs[i]
  cat(sprintf("  [%d/%d] %s", i, length(unique_drugs), drug))
  kid <- get_kegg_id(drug)
  if (is.na(kid)) {
    short <- table1$name_short[match(tolower(drug), table1$name_lower)]
    kid   <- get_kegg_id(short)
    if (!is.na(kid)) cat(" (short)")
  }
  kegg_id_df$kegg_id[i] <- kid
  cat(" →", ifelse(is.na(kid), "NOT FOUND", kid), "\n")
  Sys.sleep(0.3)
}

kegg_result <- kegg_id_df %>% left_join(brite_lookup, by = "kegg_id")

cat("\n=== KEGG Target-based Classification 結果 ===\n")
kegg_result %>% select(Name, kegg_id, kegg_cat_A, kegg_cat_B, kegg_cat_D) %>%
  print(n = 30)

# ============================================================
# 6. 全データ結合・整形・保存
# ============================================================

result <- joined %>%
  left_join(chembl_df, by = "Name") %>%
  left_join(kegg_result, by = "Name") %>%
  transmute(
    PFAS                          = PFAS,
    Drug                          = Name,
    Zdc                           = round(Zdc, 3),
    `ATC Level1`                  = paste0(atc_l1, ": ", str_to_title(atc_l1_name)),
    `ATC Level2`                  = paste0(atc_l2, ": ", str_to_title(atc_l2_name)),
    `ChEMBL ID`                   = chembl_id,
    `ChEMBL MOA`                  = chembl_moa,
    `KEGG Drug ID`                = kegg_id,
    `KEGG Target category (A)`    = kegg_cat_A,
    `KEGG Target subgroup (B)`    = kegg_cat_B,
    `KEGG Target receptor (C)`    = kegg_cat_C,
    `KEGG Drug action (D)`        = kegg_cat_D
  )

out_path <- "/path/to/your/local/directory/Table1_PFAS_drug_proximity_ATC.xlsx"
write.xlsx(result, file = out_path, rowNames = FALSE)
cat("\nSaved:", out_path, "\n")
