###################
#  data check     #
###################

### 1. load ###############
PFAS_drugs    <- read.table("/path/to/your/local/directory/summary_PFAS_appdrug_new_20230703.txt", header = F, sep = "\t", comment.char = "", quote = '"')

#drug-protein
in_f <- "/path/to/your/local/directory/CTD_chem_gene_ixns.csv"
ctd <- read.csv(paste0(in_f), skip=29, header = F) #2502394 interactions
colnames(ctd) <- c("ChemicalName",	"ChemicalID",	"CasRN",	"GeneSymbol",
                   "GeneID", "GeneForms",	"Organism",	"OrganismID",
                   "Interaction",	"InteractionActions",	"PubMedIDs")

in_f <- "/path/to/your/local/directory/drug_module.dat"
load(in_f)

in_f <- "/path/to/your/local/directory/all_chem_module.dat"
load(in_f)

### 2. check ###############
ctd_human <- ctd[ctd$OrganismID == "9606",]

sum(unique(ctd$ChemicalName) %in% unique(PFAS_drugs$V2)) #全部入っている
sum(unique(ctd_human$ChemicalName) %in% unique(PFAS_drugs$V2)) #全部入っている

### 3. extraction ###############
df_ctd_drugs <- ctd_human[(ctd_human$ChemicalName) %in% unique(PFAS_drugs$V2),]

dim(df_ctd_drugs) #[1] 271419     11　interaction　数
length(unique(df_ctd_drugs$ChemicalID)) #900
length(unique(df_ctd_drugs$GeneID)) #21441 protein 数