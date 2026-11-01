###########################
#                         #
#  disease-protein data   #
#          check          #
#                         #  
###########################


### 1. load ####################
load(file="./path/to/intermediate/data/0.make_data/list_disease_module_causalgenes.dat")
load(file="./path/to/intermediate/data/0.data_make/causal_genes/list_Disease_GeneIDs.dat")

### 2. load ####################
names(list_disease_module)
names(list_Disease_GeneIDs)

KEGGID_module <- 
c("H00003",  "H00004"  , "H00018",  "H00020",  "H00031", "H00056" ,"H00079", "H00409" ,"H01633" ,"H01358")

KEGGID_GeneIDs <- c("H00021", "H00024", "H00032", "H00048", "H00082", "H01635")

MESHID <-
c("D015470", "D015464", "D013274", "D015179", "D001943", "D000544", "D001249", "D003924", "D006973", "D003876")

KEGG_Disease_name_module <-
c("Acute myeloid leukemia",     "Chronic myeloid leukemia",     "Gastric cancer", "Colorectal cancer", "Breast cancer", 
  "Alzheimer disease", "Asthma", "Type 2 diabetes mellitus", "High blood pressure (Hypertension)", "Atopic dermatitis")

KEGG_Disease_name_GeneIDs <-
  c("Renal cell cancer", "Prostate cancer",     "Thyroid cancer", "Hepatocellular carcinoma", "Graves disease", "Dyslipidemia")

MESH_Disease_name <-
  c("Leukemia, Myeloid, Acute",     "Leukemia, Myelogenous, Chronic, BCR-ABL Positive",    
    "Stomach Neoplasms", "Colorectal Neoplasms", "Breast Neoplasms", 
    "Alzheimer Disease", "Asthma", "Diabetes Mellitus, Type 2", "Hypertension", "Dermatitis, Atopic")

### 3. Merge ####################
list_disease_EntrezIDs <- c(list_disease_module, list_Disease_GeneIDs)
names(list_disease_EntrezIDs) <- paste0(names(list_disease_EntrezIDs), "_", c(KEGG_Disease_name_module, KEGG_Disease_name_GeneIDs))

save(list_disease_EntrezIDs, file = "./path/to/intermediate/data/0.data_make/causal_genes/list_disease_EntrezIDs.dat")
