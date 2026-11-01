##############################
#                            #
#         PFAS-genes         #
#         data make          #
#                            #
##############################
library(dplyr)
library(openxlsx)

### 1.load ####################
load("./path/to/intermediate/data/0.data_make/df_EPA_PFAS_target_ctd.dat")
load("./path/to/intermediate/data/0.data_make/df_EPA_PFAS_target_Lincs.dat")

#Takemoto lab PFAS
load("./path/to/intermediate/data/test.dat")

### 2. check ####################
x <- unique(df_EPA_PFAS_target_Lincs$CASRN) %in% unique(df_EPA_PFAS_target_ctd$CASRN)
sum(x) #LincsのPFASはCTDの中あるやつと半分同じ
unique(df_EPA_PFAS_target_Lincs$CASRN)[!x]

### 3. 必要なデータをmerge ####################
#それぞれデータを分割していらない列を削除
# df_lincs <- split(df_EPA_PFAS_target_Lincs, df_EPA_PFAS_target_Lincs$CASRN)
# df_ctds  <- split(df_EPA_PFAS_target_ctd, df_EPA_PFAS_target_ctd$CASRN)

list_lincs_genes <- tapply(df_EPA_PFAS_target_Lincs$genes_xy, df_EPA_PFAS_target_Lincs$CASRN, function(x){list(sort(unique(as.integer(x))))})
list_ctd_genes   <- tapply(df_EPA_PFAS_target_ctd$GeneID, df_EPA_PFAS_target_ctd$CASRN, function(x){list(sort(unique(as.integer(x))))})

#同じ名前なので、重複を避けながらmergeする
list_PFAS_genes <- lapply(unique(c(names(list_lincs_genes), names(list_ctd_genes))), 
                          function(x){sort(unique(c(list_lincs_genes[[x]], list_ctd_genes[[x]])))})
names(list_PFAS_genes) <- unique(c(names(list_lincs_genes), names(list_ctd_genes)))

### 4. nameを入れる ####################
#Ref_name
ref_name_Lincs <- unique(df_EPA_PFAS_target_Lincs$PREFERRED.NAME)
ref_name_ctd   <- unique(df_EPA_PFAS_target_ctd$PREFERRED.NAME)
names(ref_name_Lincs) <- unique(df_EPA_PFAS_target_Lincs$CASRN)
names(ref_name_ctd) <- unique(df_EPA_PFAS_target_ctd$CASRN)

#CTD_name
name_ctd        <- unique(df_EPA_PFAS_target_ctd$ChemicalName)
names(name_ctd) <- unique(df_EPA_PFAS_target_ctd$CASRN)

#CTD_chemical ID
ChemID_ctd        <- unique(df_EPA_PFAS_target_ctd$ChemicalID)
names(ChemID_ctd) <- unique(df_EPA_PFAS_target_ctd$CASRN)

#合体
df_PFAS_genes <- bind_rows(lapply(list_PFAS_genes, function(GeneID){data.frame(GeneID)})) %>% 
  mutate(CASRN          = rep(names(list_PFAS_genes), sapply(list_PFAS_genes, length)),
         PREFERRED.NAME = rep(c(ref_name_Lincs, ref_name_ctd)[names(list_PFAS_genes)], sapply(list_PFAS_genes, length)),
         ChemicalName   = rep(name_ctd[names(list_PFAS_genes)], sapply(list_PFAS_genes, length)),
         ChemicalID     = rep(ChemID_ctd[names(list_PFAS_genes)], sapply(list_PFAS_genes, length))
         ) %>% 
  select(CASRN, everything())

#確認
head(df_PFAS_genes)
tail(df_PFAS_genes)

### 5. 保存 ###############################
list_PFAS_genes <- list(Merged = df_PFAS_genes,
                        CTD    = df_EPA_PFAS_target_ctd,
                        Lincs  = df_EPA_PFAS_target_Lincs)

save(list_PFAS_genes, file="./path/to/intermediate/data/0.data_make/human_PFAS_genes/human_PFAS_genes_230712_iida.dat")
write.xlsx(list_PFAS_genes, file="./path/to/intermediate/data/0.data_make/human_PFAS_genes/human_PFAS_genes_230712_iida.xlsx")

### 6. 確認 ################################
#竹本LabでとってたPFAS分子種との比較
y <- unique(df_EPA_PFAS_target_ctd[df_EPA_PFAS_target_ctd$CASRN %in% names(list_PFAS_genes), "ChemicalName"])

