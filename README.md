# EPIC-Immune-Deconvolution-GSE136197
Immune cell deconvolution of canine mammary tumor RNA-seq data (GSE136197) using the EPIC algorithm. Includes full R code for preprocessing, EPIC analysis, statistical comparisons, and figure generation for publication.
# Immune Cell Deconvolution Analysis in Canine Mammary Tumors Using EPIC

This repository contains R scripts, data, and visualizations used to perform immune cell deconvolution of canine mammary tumor RNA-seq data using the EPIC algorithm, with a comparative analysis between normal and tumor samples.

The pipeline includes immune cell fraction quantification, statistical testing, effect size calculation, and publication-ready visualizations including boxplots, PCA, and heatmaps.

---

## 🧬 Overview

- **Data Source**: GSE136197 (Canine mammary tumor RNA-seq data)
- **Deconvolution Tool**: [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/)
- **Cell Types**: B cells, CD4/CD8 T cells, NK cells, Macrophages, CAFs, Endothelial cells, Other immune cells
- **Groups Compared**: Tumor vs Normal
- **Visualization**: ggplot2 boxplots with Wilcoxon test p-values and Cohen’s d effect sizes

---

## 📁 Repository Structure
├── data/
│ └── EPIC_input_matrix.tsv # TPM matrix or preprocessed expression data
├── results/
│ ├── EPIC_Immune_Heatmap_Scaled.pdf # Heatmap of scaled immune fractions
│ ├── Immune_Cell_Fractions_EPIC.pdf # Final publication-ready boxplots
│ └── EPIC_PCA.pdf # PCA of immune cell fractions
├── scripts/
│ ├── 01_EPIC_processing.R # Data cleaning and formatting
│ ├── 02_EPIC_visualization.R # Main figure generation
│ └── 03_EPIC_statistics.R # Effect sizes and statistical analysis
├── README.md


---

## 📊 Main Visualizations

- 📌 **Boxplots** (Tumor vs Normal) of each immune cell type with Wilcoxon test p-values
- 📌 **Heatmap** of scaled immune fractions across all samples
- 📌 **PCA plot** of immune landscapes (color-coded by phenotype)

---

## 📦 R Package Requirements

install.packages(c("tidyverse", "ggpubr", "ggplot2", "effectsize", "pheatmap", "RColorBrewer", "factoextra"))

install Bioconductor packages if needed:

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("Biobase")

▶️ How to Run
Prepare the expression matrix
Format the matrix as required by EPIC and load it into scripts/01_EPIC_processing.R

Run immune deconvolution externally
Submit the expression matrix to EPIC and download the results.

Generate plots and statistics
Execute 02_EPIC_visualization.R and 03_EPIC_statistics.R to create all plots.

Export publication figures
Set ggsave() or pdf() output to save at dpi = 300 as PDF/PNG.

📜 Citation
If you use this code or adapt it for your research, please cite:

Herath, R. (2025). Immune Cell Landscape of Canine Mammary Tumors Inferred by EPIC Deconvolution. Submitted for publication.

📘 License
This repository is licensed under the MIT License. See LICENSE for details.

🤝 Contributions
You’re welcome to fork, cite, or adapt this work. For major contributions, please open an issue or contact the original author.

💡 Acknowledgments
EPIC developers and Gfeller Lab for providing the algorithm

The original authors of GSE136197 dataset

SLU lab team for mentorship and support during the project


---

Let me know if you'd like:

- A downloadable `README.md` file
- A ready-to-upload `.zip` of this full repo structure
- Help writing your `LICENSE` file or a `citation.cff` for GitHub's citation button

I'm happy to help!
