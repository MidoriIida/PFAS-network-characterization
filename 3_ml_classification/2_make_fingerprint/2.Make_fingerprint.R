############################
#    fingerprintの生成     #
############################
library(rcdk)
library(fingerprint)

### 1. load #################################
df = read.csv("./Analysis_250412/cas_to_smiles.csv")

### 2. get fingerprints #####################
mols <- parse.smiles(df$SMILES)

## get ECF6 (morganfinger print)
morgan <- lapply(mols, get.fingerprint, type = "circular", size = 2048)

## get MACCS keys for all the molecules
maccs <- lapply(mols, get.fingerprint, type='maccs')

## get feature, count fingerprinter
featurecounts <- lapply(mols, get.fingerprint, type='signature', fp.mode='raw')

## get Substructure fingerprint for functional group fragments
Substructure <- lapply(mols, get.fingerprint, type='substructure')

## get Pubchem fingerprint 
pubchem <- lapply(mols, get.fingerprint, type='pubchem')
