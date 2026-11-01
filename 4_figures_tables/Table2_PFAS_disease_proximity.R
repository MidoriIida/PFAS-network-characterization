####################################################
#  Table 2: PFAS-disease proximity                 #
#           + ICD-11 (disease-level)               #
#           + KEGG pathway (disease-level)         #
####################################################
#
# 各疾患についてKEGG Disease IDを取得し、
# そのエントリからICD-11コード・チャプターと
# KEGG networkパスウェイを付与する

library(readxl)
library(openxlsx)
library(dplyr)
library(stringr)
library(httr)

# ============================================================
# 1. Table4 読み込み
# ============================================================

table4 <- read_excel(
  "/path/to/your/local/directory/Table4_PFAS_disease_proximity.xlsx"
)
colnames(table4) <- c("ID", "Disease", "dc_real", "Zdc")

pfas_map <- c("335-67-1" = "PFOA", "1763-23-1" = "PFOS")
table4 <- table4 %>%
  mutate(
    cas  = str_extract(ID, "^[^|]+"),
    PFAS = if_else(is.na(pfas_map[cas]), cas, pfas_map[cas])
  )

# ============================================================
# 2. KEGG Disease ID の取得
# ============================================================

## "H00004_Chronic myeloid leukemia" → "H00004" を抽出
## 平名疾患は KEGG API で検索
get_kegg_disease_id <- function(disease_name) {
  # H0xxxx_xxx 形式 → IDを直接抽出
  hid <- str_extract(disease_name, "H\\d{5}")
  if (!is.na(hid)) return(hid)

  # 平名疾患名 → KEGG find/disease で検索
  query <- str_remove(disease_name, "\\s*\\(.*\\)$")  # 括弧以降を除去
  url   <- paste0("https://rest.kegg.jp/find/disease/", URLencode(query))
  resp  <- tryCatch(GET(url, timeout(15)), error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) return(NA_character_)
  txt   <- content(resp, "text", encoding = "UTF-8")
  if (nchar(trimws(txt)) == 0) return(NA_character_)

  lines <- strsplit(txt, "\n")[[1]]
  lines <- lines[nchar(trimws(lines)) > 0]

  # 完全一致の行を優先
  exact <- lines[tolower(str_extract(lines, "(?<=\\t).*")) == tolower(query)]
  if (length(exact) > 0) return(str_extract(exact[1], "H\\d{5}"))

  # 最初のヒット
  str_extract(lines[1], "H\\d{5}")
}

## ============================================================
## 3. KEGG Disease エントリのパース
## ============================================================

## KEGGエントリから ICD-11 と KEGG pathway を抽出
parse_kegg_disease <- function(kegg_id) {
  if (is.na(kegg_id)) {
    return(list(icd11_code = NA, icd11_chapter = NA, kegg_pathway = NA))
  }

  url  <- paste0("https://rest.kegg.jp/get/", kegg_id)
  resp <- tryCatch(GET(url, timeout(15)), error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) {
    return(list(icd11_code = NA, icd11_chapter = NA, kegg_pathway = NA))
  }
  txt   <- content(resp, "text", encoding = "UTF-8")
  lines <- strsplit(txt, "\n")[[1]]

  ## ICD-11コード: "DBLINKS     ICD-11: 2A20.0"
  icd_line <- lines[grepl("ICD-11:", lines)]
  icd11_code <- if (length(icd_line) > 0)
    paste(str_extract_all(icd_line, "[0-9A-Z]{4,}(?:\\.[0-9]+)?")[[1]], collapse = "; ")
  else NA_character_

  ## ICD-11チャプター: BRITEセクション内の "XX Neoplasms" 形式の行
  brite_start <- which(grepl("^BRITE", lines))
  brite_end   <- which(grepl("^[A-Z]", lines))
  brite_end   <- brite_end[brite_end > brite_start[1]][1] - 1
  if (length(brite_start) > 0 && !is.na(brite_end)) {
    brite_lines <- lines[brite_start:brite_end]
    chapter_line <- brite_lines[grepl("^\\s+\\d{2}\\s+", brite_lines)]
    icd11_chapter <- if (length(chapter_line) > 0)
      paste(unique(str_trim(chapter_line)), collapse = "; ")
    else NA_character_
  } else {
    icd11_chapter <- NA_character_
  }

  ## KEGG network pathway: "NETWORK     nt06276 Chronic myeloid leukemia"
  ## または BRITEの "Pathway-based classification" セクション内の "nt0xxxx" 行
  net_line <- lines[grepl("^NETWORK", lines)]
  if (length(net_line) > 0) {
    kegg_pathway <- paste(
      str_trim(str_remove_all(net_line, "^NETWORK\\s+")),
      collapse = "; "
    )
  } else {
    # BRITEのPathway-based classificationから取得
    pb_start <- which(grepl("Pathway-based classification", lines))
    if (length(pb_start) > 0) {
      pb_lines  <- lines[(pb_start+1):min(pb_start+20, length(lines))]
      nt_lines  <- pb_lines[grepl("nt\\d{5}", pb_lines)]
      kegg_pathway <- if (length(nt_lines) > 0)
        paste(unique(str_trim(nt_lines)), collapse = "; ")
      else NA_character_
    } else {
      kegg_pathway <- NA_character_
    }
  }

  list(icd11_code    = icd11_code,
       icd11_chapter = icd11_chapter,
       kegg_pathway  = kegg_pathway)
}

# ============================================================
# 4. ユニーク疾患ごとに API 照会
# ============================================================

unique_diseases <- unique(table4$Disease)
cat("照会する疾患数:", length(unique_diseases), "\n\n")

disease_df <- data.frame(
  Disease       = unique_diseases,
  kegg_id       = NA_character_,
  icd11_code    = NA_character_,
  icd11_chapter = NA_character_,
  kegg_pathway  = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(unique_diseases)) {
  disease <- unique_diseases[i]
  cat(sprintf("[%d/%d] %s\n", i, length(unique_diseases), disease))

  # KEGG Disease ID
  kid <- get_kegg_disease_id(disease)
  disease_df$kegg_id[i] <- kid
  cat("  KEGG ID:", ifelse(is.na(kid), "NOT FOUND", kid), "\n")

  if (!is.na(kid)) {
    parsed <- parse_kegg_disease(kid)
    disease_df$icd11_code[i]    <- parsed$icd11_code
    disease_df$icd11_chapter[i] <- parsed$icd11_chapter
    disease_df$kegg_pathway[i]  <- parsed$kegg_pathway
    cat("  ICD-11 code   :", ifelse(is.na(parsed$icd11_code),    "-", parsed$icd11_code),    "\n")
    cat("  ICD-11 chapter:", ifelse(is.na(parsed$icd11_chapter), "-", parsed$icd11_chapter), "\n")
    cat("  KEGG pathway  :", ifelse(is.na(parsed$kegg_pathway),  "-", parsed$kegg_pathway),  "\n")
  }
  Sys.sleep(0.3)  # API負荷軽減
}

cat("\n=== 疾患マッピング結果 ===\n")
print(disease_df %>% select(Disease, kegg_id, icd11_chapter, kegg_pathway))

# ============================================================
# 5. 結合・整形・保存
# ============================================================

result <- table4 %>%
  left_join(disease_df, by = "Disease") %>%
  transmute(
    PFAS                = PFAS,
    Disease             = Disease,
    Zdc                 = round(Zdc, 3),
    `KEGG Disease ID`   = kegg_id,
    `ICD-11 code`       = icd11_code,
    `ICD-11 chapter`    = icd11_chapter,
    `KEGG pathway`      = kegg_pathway
  )

cat("\n=== 最終テーブル ===\n")
print(result, n = 30)

out_path <- "/path/to/your/local/directory/Table2_PFAS_disease_proximity.xlsx"
write.xlsx(result, file = out_path, rowNames = FALSE)
cat("\nSaved:", out_path, "\n")
