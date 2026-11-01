# PFAS Network Characterization

This repository contains the analysis code for the manuscript:

**"A Network-based Framework for Characterizing Biological Activities and Potential Health Outcomes of Per- and Polyfluoroalkyl Substances (PFAS)"**

Midori Iida, Keito Takahashi, Kazuhiro Takemoto  
Kyushu Institute of Technology, Fukuoka, Japan

---

## Repository Structure

```
PFAS-network-characterization/
├── function.R                        # Shared R functions (network proximity calculation)
├── function.py                       # Shared Python utilities
│
├── 0_data_preparation/               # Input data construction (run in order 1→5; 2–4 can run in parallel)
│   ├── 1_ppi/                        # [Step 1] Human PPI network construction (reads data/g_PPI_ext.dat)
│   ├── 2_pfas_genes/                 # [Step 2] PFAS-associated gene retrieval (CTD, LINCS)
│   ├── 3_disease_genes/              # [Step 3] Disease-associated gene retrieval (OMIM, ClinVar, GWAS, PheWAS)
│   ├── 4_drug_genes/                 # [Step 4] Drug-associated gene data (CTD)
│   └── 5_modules/                    # [Step 5] Gene module construction (requires steps 1–3)
│
├── 1_drug_space/                     # Drug-space analysis
│   ├── 2.Calculate_45PFAS_on_900drug_space.R   # Network proximity: PFAS vs. drugs
│   └── enrichment_analysis/          # ATC and target enrichment
│
├── 2_disease_space/                  # Disease-space analysis
│   ├── 1.compute_network_proximity_between_PFAS_disease.R
│   ├── 2.Get_proximity_of_45PFAS.R
│   ├── 3.Get_similarity_45PFAS_148diseases.R
│   └── enrichment_analysis/          # ICD-11 and KEGG enrichment
│
├── 3_ml_classification/              # Machine learning classification
│   ├── 0.Get_SMILES.py               # Retrieve SMILES from EPA CompTox
│   ├── 1_get_structure/              # Structure retrieval
│   ├── 2_make_fingerprint/           # Morgan fingerprint generation (radius=2, 2048-bit)
│   ├── 3_full_model/                 # Random forest / XGBoost classifiers
│   ├── 4_feature_selection/          # Recursive feature elimination (RFE)
│   ├── 5_shap/                       # SHAP analysis
│   ├── 6_visualization/              # Visualization
│   ├── 7.Get_Morgan_finger_print_PFASmasterlist.py   # Fingerprints for EPA PFAS Master List
│   └── 8.Prediction_Morgan_fingerprint_PFASmasterlist.py  # Large-scale prediction
│
└── 4_figures_tables/                 # Figure and table generation scripts
    ├── figures/                      # Main figures (Figure 1–5)
    ├── supplementary_figures/        # Supplementary figures (S1–S10)
    ├── supplementary_tables/         # Supplementary tables (S1–S4)
    └── supplementary_data/           # Supplementary data scripts (Data S1–S6, S9)
```

---

## Requirements

### R packages
- `ggplot2`, `dplyr`, `umap`, `reshape2`, `cowplot`, `lemon`, `stringr`, `openxlsx`, `igraph`

### Python packages
- `rdkit`, `scikit-learn`, `xgboost`, `shap`, `pandas`, `numpy`, `matplotlib`, `umap-learn`, `requests`

---

## Data

Input data used in this study are obtained from publicly available sources:

| Data | Source |
|---|---|
| Human PPI network | APID, BioGRID, CORUM, HuRI, InnateDB, and others (see Methods) |
| PFAS-associated genes | CTD, LINCS |
| Drug-associated genes | CTD |
| Disease-associated genes | OMIM, ClinVar, GWAS Catalog, PheWAS Catalog |
| PFAS structures | EPA CompTox Chemicals Dashboard |
| Drug ATC codes / targets | KEGG, DrugBank |
| Disease classifications | ICD-11, KEGG Disease |

Intermediate data files (network proximity scores, UMAP coordinates, clustering results) required to run the visualization scripts are available from the corresponding author upon request.

---

## Notes on File Paths

Scripts in `4_figures_tables/` load pre-computed intermediate data (`.dat` files) from local paths.  
To reproduce the full analysis pipeline, run scripts in the following order:

1. `0_data_preparation/` → PPI construction, gene retrieval, module building
2. `1_drug_space/` → network proximity calculation (PFAS vs. drugs)
3. `2_disease_space/` → network proximity calculation (PFAS vs. diseases)
4. `3_ml_classification/` → fingerprints → classification → SHAP
5. `4_figures_tables/` → generate figures and tables

---

## Contact

Midori Iida, Ph.D.  
Kyushu Institute of Technology  
Email: iida.midori517@mail.kyutech.jp
