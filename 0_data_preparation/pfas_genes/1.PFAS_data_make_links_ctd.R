#########################
#                       #
#    Check the PFAS     #
#     LINCS & CTD       #
#########################
#PFASいくつ入ってる？
library(tidyverse)
library(tictoc)
library(fs)
library(vroom)
library(data.table)
library(openxlsx)

### 1. load ######################################
#CTD_chems
in_f2 <- "/Users/iida/Documents/Lab_Iida/Orignal_data/CTD/CTD_chemicals.csv"
CTD_chems <- read.csv(in_f2, skip=29, header = F)
colnames(CTD_chems) <- c("ChemicalName","ChemicalID","CasRN","Definition","ParentIDs","TreeNumbers","ParentTreeNumbers","Synonyms")

#CTD Chem-genes
load(file = "./Analysis_230701/0.data_make/ctd_hPPI.dat")

#Lincks PFAS-genes
load(file = "./Analysis_230701/0.data_make//df_EPA_PFAS_target_Lincs.dat")

#CTD PFAS-genes
load(file="./Analysis_230701/0.data_make//df_EPA_PFAS_target_ctd.dat")

### 2. Check ######################################
x <- unique(df_EPA_PFAS_target_Lincs$CASRN) %in% unique(df_EPA_PFAS_target_ctd$CASRN)
sum(x) #LincsのPFASはCTDの中あるやつと半分同じ

#これらのデータはdf_EPA_PFAS_target_ctdにはない
unique(df_EPA_PFAS_target_Lincs$CASRN)[!x]
#[1] "293753-05-6" "847925-91-1" "103055-07-8" "467459-31-0"
unique(df_EPA_PFAS_target_Lincs[df_EPA_PFAS_target_Lincs$CASRN %in% unique(df_EPA_PFAS_target_Lincs$CASRN)[!x],"PREFERRED.NAME"])
# [1] "N-[4-(1,1,1,3,3,3-Hexafluoro-2-hydroxy-2-propanyl)phenyl]-2-thiophenesulfonamide"                                 
# [2] "N-((7S)-6,7-Dihydro-6-oxo-5H-dibenz(b,d)azepin-7-yl)-2,2-dimethyl-N'-(2,2,3,3,3-pentafluoropropyl)-propanediamide"
# [3] "N-[2,5-dichloro-4-(1,1,2,3,3,3-hexafluoropropoxy)-phenyl-aminocarbonyl]-2,6-difluorobenzamide"                    
# [4] "2-(6-Fluoro-1H-indol-3-yl)-N-[[3-(2,2,3,3-tetrafluoropropoxy)phenyl]methyl]ethanamine" 
 
### 3. CTDのデータには入っている?######################################
#名前で検索
sum(unique(ctd_hPPI$ChemicalName) %in% unique(df_EPA_PFAS_target_Lincs[df_EPA_PFAS_target_Lincs$CASRN %in% unique(df_EPA_PFAS_target_Lincs$CASRN)[!x],"PREFERRED.NAME"]))

#CASで検索
sum(unique(ctd_hPPI$CasRN) %in% unique(df_EPA_PFAS_target_Lincs[df_EPA_PFAS_target_Lincs$CASRN %in% unique(df_EPA_PFAS_target_Lincs$CASRN)[!x],"CASRN"]))

#手動で検索
# [1] "293753-05-6" #データベース自体にはあるがhumanProtein情報はない(SR 3335)
# [2] "847925-91-1" #データベース自体にはあるがhumanProtein情報はない(2,2-dimethyl-N-(6-oxo-6,7-dihydro-5H-dibenzo(b,d)azepin-7-yl)-N'-(2,2,3,3,3-pentafluoropropyl)malonamide)
# [3] "103055-07-8" #データベース自体にはあるがhumanProtein情報はない(fluphenacur)
# [4] "467459-31-0" #データベース自体にはあるがhumanProtein情報はない((2-(6-fluoro-1H-indol-3-yl)-ethyl)-(3-(2,2,3,3-tetrafluoropropoxy)benzyl)amine)





















