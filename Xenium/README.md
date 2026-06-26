# Xenium Spatial Transcriptomics — CUP/IC Analysis Pipeline

This repository contains the analysis pipeline for spatial transcriptomics profiling of Cancer of Unknown Primary (CUP) and Identified Cancer (IC) samples using the **10x Genomics Xenium** platform. The cohort comprises **16 samples** (11 CUP + 5 IC) processed with the **Xenium 5K panel**. The five IC samples additionally have MultiTissue Xenium panel (541 genes) data available.

Core analytical steps include: BPCells-backed single-cell processing, CancerFinder3 (CF3) tumour/stroma classification, negative-control-probe and GMM-based expression thresholding, Spatial H-Score computation, tissue-of-origin scoring, PROGENy pseudobulk pathway analysis, UCDBase immune deconvolution, and SingleR cell-type annotation.

---

## Scripts

| Script | Purpose | Key inputs | Key outputs |
|--------|---------|------------|-------------|
| `01_extract_cohort_data.R` | Load raw Xenium data via BPCells for 16 samples; QC filter; log-normalise; export MTX → h5ad; run CancerFinder3 (CF3) tumour/stroma inference; sketch PCA/UMAP; save per-sample slim data frames and cohort RDS | `sample_sheet.csv`, raw Xenium output folders, CF3 model + scripts | `cohort_slim.rds`, per-sample `*_slim.rds` |
| `02_cohort_visualisation.qmd` | Cohort-level Quarto report: violin plots, dot plots, heatmaps and Wilcoxon tests for actionable targets and ICI markers split by Tumour/Stroma compartment | `cohort_slim_v2.rds` | HTML/PDF report, PNG figures |
| `03_threshold_classification.R` | Three-phase cell-level thresholding (Negative Control Probes, Stromal Floor 95th percentile, GMM via mclust); Spatial H-Score per sample; K-means (gap statistic) classification into Negative/Low/High | `cohort_slim_v2.rds`, BPCells on-disk matrices | `classification_results.csv`, H-Score and K-means tables |
| `04_classification_methods_report.qmd` | Quarto report documenting and visualising the classification methodology (NegCtrl, Stromal Floor, GMM, H-Score, K-means) | Outputs from script 03 | HTML/PDF report |
| `05_tissue_origin_immune_classification.R` | (A) Weighted tissue-of-origin scoring for 12 candidate primary sites using Spatial H-Score; (B) Immune signature scoring (Effector Cell Traffic, Antitumor Cytokines, ICI) and composite Immune Activity Score; combined 37-marker classification | H-Score tables from script 03 | `tissue_origin_scores.csv`, `immune_classification.csv` |
| `06_PI_report.qmd` | Comprehensive PI-facing Quarto report: benchmarks Xenium tissue-of-origin and ONC-CUP/IMI-CUP classification against VisiumST, dkfz, OncoNPC, and EPICUP; reads clinical data from Excel | Script 03/05 outputs, PROGENy scores, clinical Excel file | HTML/PDF report |
| `07_progeny_pathway_analysis.R` | Pseudobulk PROGENy pathway activity for 7 pathways (PI3K, EGFR, MAPK, VEGF, JAK-STAT, NFkB, TNFa); CPM-normalised pseudobulk from BPCells; hierarchical clustering into ONC-CUP vs IMI-CUP | `cohort_slim_v2.rds`, BPCells on-disk matrices | `progeny_scores.csv`, classification heatmap |
| `08_spatial_pathway_maps.R` | Per-sample spatial maps of PROGENy pathway activity using weighted gene scoring from BPCells (no full matrix materialisation); 7-panel PNGs with dark background | BPCells on-disk matrices, PROGENy weights | Per-sample PNG spatial maps |
| `09_immune_infiltration_markers.R` | Characterises immune infiltration in the tumour compartment; extracts 13 immune genes (B cell, DC, macrophage, NK, neutrophil markers) from BPCells; H-Score per population; enrichment classification | BPCells on-disk matrices, `cohort_slim_v2.rds` | `immune_infiltration.csv`, spatial and summary plots |
| `10_ucdbase_deconvolution.py` | UCDBase single-cell deconvolution on 10 k sketch cells per sample; authenticates with UCDBase API token; extracts 8 immune population fractions; aggregates to sample level; spatial immune maps | Sketch h5ad files from script 01 | `ucdbase_fractions.csv`, spatial proportion plots |
| `11_singleR_copykat_annotation.R` | Per-cell annotation: tumour cells labelled by tissue origin from script 05; stroma cells annotated by SingleR (HumanPrimaryCellAtlasData, coarse labels); spatial annotation plots | `cohort_slim_v2.rds`, BPCells on-disk matrices | `cell_annotations.csv`, spatial annotation PNGs |
| `11_Classification_Comparison.qmd` | Quarto report benchmarking Xenium origin predictions vs VisiumST/dkfz/OncoNPC/EPICUP; ONC/IMI PROGENy v3 (mixed-compartment) vs clinical comparison; reads Excel clinical data | Script 03/05/07 outputs, clinical Excel file | HTML/PDF report |
| `12_emt_stemness_IC.R` | IC-only analysis: EMT and stemness scoring from MultiTissue Xenium panel; adds ALDH1A3 (T1 progenitor) and LY6D (T2 inflammatory); z-scored EMT-T1, EMT-T2, and Stemness scores across 5 IC samples; spatial plots with H&E thumbnails | MultiTissue Xenium output folders, H&E TIF images | `emt_stemness_IC.csv`, spatial score plots |
| `12b_tumor_origin_IC.R` | Tumour-of-origin prediction for IC samples: CK7/CK20 primary classification + 13 lineage-specific markers; 5K panel (20 markers in tumour) + MultiTissue panel (KRT7, HMGCS2); binary positivity + likelihood scoring against 13 cancer type profiles | MultiTissue Xenium output folders, H&E TIF images, script 03 thresholds | `tumor_origin_IC.csv`, origin likelihood heatmaps |

---

## Prerequisites

### R packages

Install from CRAN unless noted otherwise.

```r
# CRAN
install.packages(c(
  "Seurat",        # v5 required
  "BPCells",       # on-disk matrix backend
  "Matrix",
  "dplyr", "tidyr", "tibble",
  "ggplot2", "patchwork", "scales", "ggrepel",
  "ggbeeswarm", "rstatix", "ggpubr",
  "knitr", "kableExtra",
  "mclust",        # GMM thresholding
  "cluster",       # K-means gap statistic
  "pheatmap", "RColorBrewer",
  "readxl",
  "png", "grid", "gridExtra"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "SingleR", "celldex",
  "BiocParallel",
  "AnnotationDbi", "org.Hs.eg.db",
  "edgeR"
))

# GitHub / specialised
remotes::install_github("saezlab/progeny")    # PROGENy pathway analysis
remotes::install_github("saezlab/decoupleR") # decoupleR spatial scoring
```

### Python packages

Tested with Python 3.9+ inside the `scf` conda environment (required for CancerFinder3):

```bash
conda activate scf
pip install numpy pandas scanpy matplotlib seaborn ucdeconvolve
```

---

## Data

### Xenium output folders

Each sample must have a standard 10x Xenium output directory containing at least:

```
<sample_name>/
  cell_feature_matrix/          # MTX sparse matrix
  cells.csv.gz                  # cell centroids and metadata
  cell_boundaries.parquet       # (optional, for polygon overlays)
```

Provide a **`sample_sheet.csv`** in the working directory (`ANALYSIS_ROOT`) with the following columns:

| Column | Description |
|--------|-------------|
| `sample_id` | Unique sample identifier |
| `xenium_dir` | Path to the Xenium output directory |
| `group` | `CUP` or `IC` |
| `...` | Any additional clinical metadata columns |

### CancerFinder3 (CF3)

Download from <https://github.com/JiangBioLab/CancerFinder> and update the three placeholder paths in `01_extract_cohort_data.R`:

```r
PY_MTX2H5AD <- "path/to/CancerFinder3/mtx_to_h5ad.py"
PY_INFER    <- "path/to/CancerFinder3/infer.py"
CF3_MODEL   <- "path/to/CancerFinder3/models/st_pretrain_article.pkl"
```

### IC MultiTissue panel (scripts 12 and 12b only)

Set `MT_ROOT` to the parent directory containing one Xenium output subfolder per IC sample and a subdirectory `HE_POST_XENIUM_20X/` holding 20x H&E TIF images:

```r
MT_ROOT <- "path/to/IC_MultiTissue"
# Expected layout:
#   MT_ROOT/output-<run_id_IC_002>/
#   MT_ROOT/output-<run_id_IC_004>/
#   ...
#   MT_ROOT/HE_POST_XENIUM_20X/STR1337.tif
#   MT_ROOT/HE_POST_XENIUM_20X/STR1338.tif
#   ...
```

---

## Usage

All scripts are designed to be run from the repository root (`ANALYSIS_ROOT <- "."`). Run them in order.

```bash
# R analysis scripts
Rscript 01_extract_cohort_data.R
Rscript 03_threshold_classification.R
Rscript 05_tissue_origin_immune_classification.R
Rscript 07_progeny_pathway_analysis.R
Rscript 08_spatial_pathway_maps.R
Rscript 09_immune_infiltration_markers.R
Rscript 11_singleR_copykat_annotation.R
Rscript 12_emt_stemness_IC.R       # IC samples only
Rscript 12b_tumor_origin_IC.R      # IC samples only

# Quarto reports (requires quarto >= 1.3)
quarto render 02_cohort_visualisation.qmd
quarto render 04_classification_methods_report.qmd
quarto render 06_PI_report.qmd
quarto render 11_Classification_Comparison.qmd

# Python deconvolution (activate scf conda env first)
conda activate scf
python 10_ucdbase_deconvolution.py --analysis_root .
```
