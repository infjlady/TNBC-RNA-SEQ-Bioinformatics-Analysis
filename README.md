# 🧬 Transcriptomic Analysis of Carboplatin Resistance in TNBC

> Differential gene expression analysis workflow for investigating transcriptomic alterations associated with carboplatin resistance in Triple-Negative Breast Cancer (TNBC).

![R Version](https://img.shields.io/badge/R-4.6.0-blue.svg)
![DESeq2](https://img.shields.io/badge/Bioconductor-DESeq2-green.svg)
![License](https://img.shields.io/badge/License-MIT-purple.svg)

## 📌 Project Overview

This repository contains a complete pipeline for analyzing bulk RNA-seq data to uncover transcriptomic remodeling associated with chemotherapy resistance. Using DESeq2 in R, we analyzed the **GEO dataset GSE309617** to identify differentially expressed genes between carboplatin-sensitive and carboplatin-resistant TNBC tumors.

The workflow encompasses:
- Data preprocessing & metadata generation
- Differential expression analysis (DESeq2)
- Advanced visualizations (Volcano plots, Heatmaps, PCA)
- Preparation for downstream systems biology analysis

---

## 📊 Dataset Information

- **GEO Accession:** [GSE309617](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309617)
- **Title:** Therapeutic Synergy Overcomes Carboplatin Resistance in Triple-Negative Breast Cancer
- **Organism:** *Homo sapiens*
- **Platform:** Illumina NovaSeq X Plus

---

## 🔬 Methodology & Workflow

Our RNA-seq analysis pipeline follows these core steps:

1. **Data Ingestion:** Import count matrix and generate metadata
2. **Modeling:** Construct DESeq2 object and perform differential expression analysis
3. **Extraction:** Filter and extract significant Differentially Expressed Genes (DEGs)
4. **Visualization:** Generate high-quality, publication-ready plots (Volcano, PCA, Heatmap)
5. **Downstream Prep:** Format data for subsequent pathway enrichment

---

## 📂 Repository Structure

```text
├── data/                    # Raw count matrices and metadata
├── scripts/
│   └── deseq2_analysis.R    # Main DESeq2 pipeline script
├── results/                 # Output CSVs and lists of significant DEGs
│   ├── volcano_plot.png     # Example volcano plot
│   ├── heatmap.png          # Example clustering heatmap
│   └── PCA_PLOT.png         # Example PCA separation plot
├── requirements.txt         # R dependencies
├── results_summary.txt      # High-level summary of analysis findings
└── README.md                # Project documentation
```

*(Note: Data and result files are generated via the `scripts/deseq2_analysis.R` pipeline.)*

---

## 🚀 Key Findings

- **Distinct Clustering:** Resistant and sensitive TNBC tumors exhibited clearly distinct transcriptomic clustering profiles.
- **Significant Alterations:** Multiple DEGs associated with therapeutic adaptation and tumor progression were identified.
- **Pathway Implications:** Preliminary analyses point toward major remodeling in cancer survival signaling pathways.

For a more detailed breakdown, please see the `results_summary.txt` file.

---

## 🛠️ Requirements & Installation

Analysis was conducted using **R version 4.6.0**. 

Key packages include:
- `DESeq2` (Bioconductor)
- `GEOquery`
- `tidyverse` (dplyr, etc.)
- `pheatmap`
- `EnhancedVolcano`

A complete list of dependencies can be found in `requirements.txt`.

---

## 🔮 Future Directions

- **Pathway Enrichment:** GO Biological Process and KEGG pathway analyses.
- **Network Systems Biology:** Hub gene identification using STRING and Cytoscape.
- **Multi-omics Integration:** Integration with epigenetic (DNA methylation) datasets.
- **Structural Biology:** Molecular docking and MD simulation studies on target proteins.

---

## ✍️ Author
**Sakshi Swami**  
*Bioinformatics | Transcriptomics | Computational Biology*

## 📜 License
This project is licensed under the [MIT License](LICENSE).
