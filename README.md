# Immune Cell Deconvolution in Canine Mammary Tumors (GSE136197)

This repository contains the complete analysis pipeline for immune cell deconvolution of canine mammary tumor RNA-seq data (GSE136197), using two complementary algorithms: **EPIC** and **CIBERSORTx**. The goal of this project is to profile immune infiltration patterns, with a specific focus on **CD4⁺ T cell enrichment** across histological subtypes of canine mammary tumors.

---

## 🧬 Overview

- **Dataset**: [GSE136197](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136197) – RNA-seq from normal and tumor canine mammary tissues  
- **Tools Used**:
  - [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/) – bulk immune deconvolution
  - [CIBERSORTx](https://cibersortx.stanford.edu/) – subset-level immune cell profiling
- **Main Goal**: Evaluate immune remodeling (esp. CD4⁺ T cells) across tumor progression
- **Study Status**: Submitted for publication (Herath, 2025)

---

## 📁 Repository Structure

├── data/
│ ├── GSE136197_TPM_humanOrthologs.tsv # Normalized expression matrix with human gene symbols
│ ├── CIBERSORTx_Output.tsv # Output matrix from CIBERSORTx
│ └── metadata.tsv # Sample metadata with histology and phenotype
│
├── results/
│ ├── EPIC_Immune_Fractions_Boxplots.pdf # Final boxplots (EPIC)
│ ├── EPIC_Heatmap_Scaled.pdf # Heatmap of z-score scaled EPIC fractions
│ ├── EPIC_PCA.pdf # PCA of immune profiles (EPIC)
│ └── CIBERSORTx_Boxplots.pdf # CD4+ subsets by histology (CIBERSORTx)
│
├── scripts/
│ ├── 01_Preprocessing_EPIC.R # Data formatting and TPM conversion
│ ├── 02_Visualization_EPIC.R # ggplot2, pheatmap, PCA
│ ├── 03_Statistics_EPIC.R # Kruskal-Wallis, Dunn test, Cohen’s d
│ └── 04_Visualization_CIBERSORTx.R # CD4⁺ subset plotting
│
├── LICENSE
└── README.md

yaml
Copy
Edit

---

## 📊 Main Visualizations

- 📌 **Boxplots**: Immune fractions (EPIC & CIBERSORTx) by histological subtype  
- 📌 **Heatmap**: Z-score scaled immune landscape across 48 samples  
- 📌 **PCA**: Dimensional reduction based on immune profiles  
- 📌 **Subset Analysis**: CD4⁺ memory, activated, Tregs, and Tfh (CIBERSORTx)

---

## 📦 Required R Packages

```r
install.packages(c("tidyverse", "ggpubr", "ggplot2", "effectsize", 
                   "pheatmap", "RColorBrewer", "factoextra", "FSA"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("Biobase")
▶️ How to Reproduce
  -Prepare expression matrix
  -Convert raw counts to TPM
  -Map canine Ensembl IDs to human orthologs
  -Format as EPIC/CIBERSORTx input

* Run EPIC deconvolution
  -Use local R script or EPIC Shiny App
  -Save output as .tsv

* Run CIBERSORTx
  -Upload TPM matrix to CIBERSORTx
  -Download immune subset fractions

* Visualize results
  -Run 02_Visualization_EPIC.R and 04_Visualization_CIBERSORTx.R
  -Export plots in .pdf format for publication

📜 Citation
If you use this repository or build upon it, please cite:

Herath, R.M.H.H. (2025). Dual Deconvolution of Canine Mammary Tumors Reveals Robust CD4⁺ T Cell Enrichment via EPIC and CIBERSORTx. Submitted to Genomics & Informatics.

📘 License
This project is licensed under the MIT License. See the LICENSE file for details.

💡 Acknowledgments
The Gfeller Lab and EPIC developers

The CIBERSORTx team at Stanford University

The authors of the GSE136197 dataset (Graim et al., 2021)

SLU Veterinary Faculty for supervision and support

🤝 Contributions
Feel free to fork, reuse, or cite this repository. For major suggestions or collaboration requests, please contact:

📧 ruwiniherath92@gmail.com
🔗 ORCID: 0009-0000-7982-0690

