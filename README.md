# PFAS Network Characterization

This repository contains the analysis code for the manuscript:

**"A Network-based Framework for Characterizing Biological Activities and Potential Health Outcomes of Per- and Polyfluoroalkyl Substances (PFAS)"**

Midori Iida, Keito Takahashi, Kazuhiro Takemoto  
Kyushu Institute of Technology, Fukuoka, Japan

---

## Repository Structure

```
PFAS-network-characterization/
├── function.R                          # Shared R functions (network proximity calculation)
├── function.py                         # Shared Python utilities
│
├── data/                               # Input data files (Supplementary Data 1–9, KEGG BRITE)
│
├── 0_data_preparation/                 # Input data construction (run in order 1→4; 2–4 can run in parallel)
│   ├── 1.PPI_data_load.R               # [Step 1] Human PPI network → g_comp.dat
│   ├── 2.PFAS_genes_load.R             # [Step 2] PFAS-associated genes → list_PFAS_genes.dat
│   ├── 3.drug_data_load.R              # [Step 3] Drug-associated genes → df_ctd_drugs.dat
│   └── 4.disease_data_load.R           # [Step 4] Disease-associated genes → list_disease_genes.dat
│
├── 1_drug_space/                       # Drug-space analysis
│   ├── 1.compute_network_proximity_between_PFAS_drug.R   # Network proximity → list_dist_mat_drug_PFAS.dat
│   ├── 2.Get_proximity_of_45PFAS_drug.R                  # z-score matrix + UMAP → list_umap_PFAS_drug.dat
│   ├── 3.Clustering.R                                    # Hierarchical clustering → list_cluster.dat / list_cluster.csv
│   └── enrichment_analysis/
│       ├── 4.Enrichment_analysis_ATC.R                   # ATC drug class enrichment
│       └── 5.Enrichment_analysis_target.R                # Drug target enrichment
│
├── 2_disease_space/                    # Disease-space analysis
│   ├── 1.compute_network_proximity_between_PFAS_disease.R  # Network proximity → list_dist_mat_disease_PFAS.dat
│   ├── 2.Get_proximity_of_45PFAS_disease.R                 # z-score matrix + UMAP → list_umap_45PFAS_148disease.dat
│   ├── 3.Clustering.R                                      # Hierarchical clustering → list_cluster_disease.dat
│   └── enrichment_analysis/
│       ├── 4.Enrichment_anaysis_ICD-11.R                   # ICD-11 disease enrichment
│       └── 5.Enrichment_anaysis_KEGG_pathway.R             # KEGG pathway enrichment
│
├── 3_ml_classification_drug/           # ML classification of PFAS by drug-space cluster (Random Forest)
│   ├── 0.Get_SMILES.py                 # Retrieve SMILES from PubChem (uses list_cluster.csv)
│   ├── 1.Get_structure.py              # Draw chemical structures by cluster
│   ├── 2.Onehotencoding_Morgan_fingerprint_2048.py   # Morgan fingerprints (radius=2, 2048-bit) → pkl + CSV
│   ├── 3.Classification_RF_Morgan.py                 # Random Forest classifier (all features)
│   ├── 4_feature_selection/
│   │   ├── 4.1.Feature_Selection_cor_Morgan.py       # Remove zero-variance and correlated features
│   │   ├── 4.2.Feature_Selection_RFE_Morgan.py       # Recursive Feature Elimination (RFE)
│   │   └── 4.3.Feature_Selection_RFE_Morgan_accuracy.py  # 5-fold CV accuracy by feature count
│   ├── 5.SHAP_Morgan_classification.py               # SHAP analysis on selected features
│   ├── 6.Visualization_Morgan.py                     # Highlight important substructures
│   ├── 7.Get_Morgan_finger_print_PFASmasterlist.py   # Fingerprints for EPA PFAS Master List
│   └── 8.Prediction_Morgan_fingerprint_PFASmasterlist.py  # Large-scale prediction
│
├── 3_ml_classification_disease/        # ML classification of PFAS by disease-space cluster (XGBoost)
│   ├── 0.Get_SMILES.py                 # Retrieve SMILES from PubChem (uses list_cluster_disease.csv)
│   ├── 1.Get_Morgan_fingerprint_chem.py  # Morgan fingerprints + 23 toxicity descriptors → pkl
│   ├── 2.Classification_XGBoost_Morgan_chems.py  # XGBoost classifier (CV=0.772) → model pkl
│   ├── 3.SHAP_XGBoost_Morgan_chems.py            # SHAP analysis on all features
│   ├── 4.Visualization_Morgan.py                 # Highlight important substructures → Figs S7, S8
│   ├── 5.Get_Morgan_fingerprint_PFASmasterlist.py  # Fingerprints + descriptors for EPA Master List
│   └── 6.Prediction_XGBoost_PFASmasterlist.py     # Large-scale prediction
│
└── 4_figures_tables/                   # Figure and table generation scripts
    ├── figures/                        # Main figures (Figure 1–5)
    ├── tables/                         # Main tables (Table 1–2)
    ├── supplementary_figures/          # Supplementary figures
    ├── supplementary_tables/           # Supplementary tables
    └── supplementary_data/             # Supplementary data scripts
```

---

## Requirements

### R packages
- `ggplot2`, `dplyr`, `umap`, `Rtsne`, `reshape2`, `cowplot`, `lemon`, `stringr`, `igraph`, `openxlsx`, `jsonlite`, `randomForest`, `caret`

### Python packages
- `rdkit`, `scikit-learn`, `xgboost`, `shap`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `pubchempy`

---

## Data

Input data used in this study are obtained from publicly available sources:

| Data | Source |
|---|---|
| Human PPI network | APID, BioGRID, CORUM, HuRI, InnateDB, and others (see Methods) |
| PFAS-associated genes | CTD, LINCS |
| Drug-associated genes | CTD |
| Disease-associated genes | OMIM, ClinVar, GWAS Catalog, PheWAS Catalog |
| PFAS structures (45 compounds) | EPA CompTox Chemicals Dashboard |
| Drug ATC codes / targets | KEGG, DrugBank |
| Disease classifications | ICD-11, KEGG Disease (br08402, br08403) |

### EPA PFAS Master List (required for `3_ml_classification_drug/7.`, `8.` and `3_ml_classification_disease/5.`, `6.`)

Download the PFAS Master List from the EPA CompTox Chemicals Dashboard:

1. Go to: https://comptox.epa.gov/dashboard/chemical-lists/PFASMASTER
2. Click **Export** → **CSV**
3. Save the file and update the path in the relevant scripts:
   ```python
   masterlist = pd.read_csv('/path/to/your/local/directory/Chemical_List_pfasmaster.csv')
   ```

The merged classification results (drug-space + disease-space, 4 groups) for all EPA PFAS compounds are provided as `data/Supplementary_Data_9_Classification_results_for_all_EPA_PFAS.txt`.
The R script to reproduce this file is `4_figures_tables/supplementary_data/Supplementary_Data_7_Classification_results_for_all_EPA_PFAS.R`.

Intermediate data files (network proximity scores, UMAP coordinates, clustering results) required to run the visualization scripts are available from the corresponding author upon request.

---

## Execution Order

Run scripts in the following order to reproduce the full analysis:

1. `0_data_preparation/` (1→4) — build PPI network and gene lists
2. `1_drug_space/` (1→3, then enrichment_analysis/) — drug proximity and clustering
3. `2_disease_space/` (1→3, then enrichment_analysis/) — disease proximity and clustering
4. `3_ml_classification_drug/` (0→8) — Random Forest classifier using drug-space clusters
5. `3_ml_classification_disease/` (0→6) — XGBoost classifier using disease-space clusters
6. `4_figures_tables/` — generate figures and tables

Intermediate `.dat` files are saved to a local directory of your choice; update `./path/to/intermediate/data/` placeholders in each script accordingly.

---

## Contact

Midori Iida, Ph.D.  
Kyushu Institute of Technology  
Email: iida.midori517@mail.kyutech.jp
