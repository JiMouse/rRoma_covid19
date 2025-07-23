# rRoma_covid19

## 🧬 Pseudobulk Seurat pipeline for COVID-19 single-cell datasets

This repository provides a reproducible, automated pipeline to build **Seurat v5** objects from publicly available **COVID-19 scRNA-seq datasets**, perform **cell type annotation**, and extract **pseudobulk counts** for downstream analysis. The initial objective was to apply the **ROMA** framework to infer regulatory module activity on aggregated expression profiles. However, the pipeline is **modular and generalizable** to any pseudobulk-based analysis or differential expression approach.

---

## 🚀 Features

- 📦 **Seurat v5 object construction** from raw matrix or HDF5 input (several example, 1.xxx scripts)
- 🧬 **Cell type annotation** using SingleR and curated references (3. singleR_annotation.R)
- 🧮 **Pseudobulk matrix generation** (gene × sample × cell type) (both raw and normalised counts, 4. seurat_pseudobulk_extraction.R)
- 🔁 Optional **batch-wise processing** and merging for multi-sample studies (5. merge batch pseudobulk.R)
- 🧰 Compatible with downstream pseudobulk tools (e.g., DESeq2, edgeR, ROMA) (6. rRoma_analysis.R as example)

---

## 📂 Included Processed Datasets

We provide preprocessed Seurat v5 objects and pseudobulk matrices for **five publicly available COVID-19 datasets**:

| Dataset | GEO ID | Tissue | Format | Cells |
|--------|--------|--------|--------|-------|
| Liao et al. 2020 | GSE145926 | BALF | Seurat v5 | ~43,000 |
| Wilk et al. 2020 | GSE150728 | PBMC | Seurat v5 | ~60,000 |
| Lee et al. 2020 | GSE157344 | Nasopharyngeal | Seurat v5 | ~31,000 |
| Stephenson et al. 2021 | GSE171524 | Whole blood | Seurat v5 | ~95,000 |
| Unterman et al. 2022 | GSE158055 | BALF + PBMC | Seurat v5 | ~47,000 |

> 🔎 Note: The datasets vary in format and annotation quality. See individual folders for dataset-specific notes.

---

## 📦 Folder Structure

```
rRoma_covid19/
├── scripts/                  # Core R scripts for import, processing, pseudobulk
├── data/                    # Raw and processed data (Seurat objects)
├── metadata/                # Cell type labels, references, annotations
├── pseudobulk_outputs/      # Output matrices (pseudocounts per cell type/sample)
├── results/                 # UMAPs, cluster stats, QC summaries
└── README.md
```

---

## 🛠️ Dependencies

- R (≥ 4.3)
- [Seurat v5](https://satijalab.org/seurat/)
- SingleR
- edgeR / DESeq2 (optional)
- Bioconductor packages: `scRNAseq`, `celldex`, `SummarizedExperiment`

Install required R packages:

```r
install.packages("Seurat")
BiocManager::install(c("SingleR", "celldex", "scRNAseq", "SummarizedExperiment"))
install.packages('devtools')
devtools::install_url("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.7.tar.gz")

```

---

## ⚙️ How to Use

### 1. Build Seurat Objects

Use `seurat_construction.R` and ` 2. seurat_features_annotation.R`  to construct a Seurat object from raw count matrices (MTX or H5).

### 2. Annotate Cell Types

3. singleR_annotation.R

### 3. Generate Pseudobulk Matrices

4. seurat_pseudobulk_extraction.R
---

## 🔬 ROMA or Custom Downstream Analysis

The output pseudobulk matrices (gene × sample × cell type) can be used as input to:

- **ROMA** (for regulatory module activity inference)
- **DESeq2/edgeR** (for pseudobulk DE analysis)
- **WGCNA**, pathway scoring, or other bulk-like techniques

Scripts to facilitate this are under development.

---

## 📜 Citation

If you use this pipeline or processed datasets, please cite this repository:

> Jean-Marie Ravel, *rRoma_covid19: A reproducible pipeline for pseudobulk analysis of COVID-19 scRNA-seq datasets*, 2025. [GitHub](https://github.com/JiMouse/rRoma_covid19)

---

## 💡 Contributing / Feedback

This is a work in progress. Contributions, pull requests, or bug reports are very welcome!
