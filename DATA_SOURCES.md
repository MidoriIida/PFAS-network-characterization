# Data Sources

This document lists the publicly available databases and resources used to construct the input data for this study. All data were downloaded as of 2023–2024.

---

## Human Protein-Protein Interaction (PPI) Network

The human interactome was constructed by integrating protein-protein interactions from the following databases:

| Database | URL | Reference |
|---|---|---|
| APID | http://apid.dep.usal.es/ | Alonso-López et al., 2016 |
| BioGRID | https://thebiogrid.org/ | Oughtred et al., 2021 |
| CORUM | https://mips.helmholtz-muenchen.de/corum/ | Giurgiu et al., 2019 |
| HuRI | http://www.interactome-atlas.org/ | Luck et al., 2020 |
| InnateDB | https://www.innatedb.com/ | Breuer et al., 2013 |
| IntAct | https://www.ebi.ac.uk/intact/ | Kerrien et al., 2012 |
| MINT | https://mint.bio.uniroma2.it/ | Licata et al., 2012 |
| Reactome | https://reactome.org/ | Jassal et al., 2020 |

---

## PFAS-Associated Genes

| Database | URL | Data used |
|---|---|---|
| CTD (Comparative Toxicogenomics Database) | https://ctdbase.org/ | Chemical–gene interactions for PFAS |
| LINCS L1000 | https://lincsproject.org/ | Gene expression profiles for PFAS |

---

## Drug-Associated Genes

| Database | URL | Data used |
|---|---|---|
| CTD (Comparative Toxicogenomics Database) | https://ctdbase.org/ | Chemical–gene interactions for approved drugs |
| KEGG Drug | https://www.genome.jp/kegg/drug/ | ATC codes for approved drugs |
| DrugBank | https://go.drugbank.com/ | Drug target information |

---

## Disease-Associated Genes

| Database | URL | Data used |
|---|---|---|
| OMIM | https://www.omim.org/ | Disease–gene associations |
| ClinVar | https://www.ncbi.nlm.nih.gov/clinvar/ | Clinically relevant variants |
| GWAS Catalog | https://www.ebi.ac.uk/gwas/ | GWAS-identified gene associations |
| PheWAS Catalog | https://phewascatalog.org/ | Phenome-wide association data |

---

## PFAS Chemical Structures

| Database | URL | Data used |
|---|---|---|
| EPA CompTox Chemicals Dashboard | https://comptox.epa.gov/dashboard/ | SMILES structures for EPA PFAS Master List |

---

## Disease Classifications

| Database | URL | Data used |
|---|---|---|
| ICD-11 | https://icd.who.int/en | Disease classification for enrichment analysis |
| KEGG Disease | https://www.genome.jp/kegg/disease/ | KEGG pathway-disease associations |

---

## Notes

- Input data files used by scripts in `0_data_preparation/` must be downloaded from the databases listed above and placed in your local directory.
- Intermediate computed data files (`.dat` files: network proximity scores, UMAP coordinates, clustering results) are not included in this repository due to file size. These are available from the corresponding author upon request.
- Supplementary data files (`.txt`) in the `data/` folder are the final analysis outputs provided as supplementary materials to the manuscript.
