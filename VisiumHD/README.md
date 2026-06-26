# CUP Visium HD Pipeline

Analysis pipeline for 3 Cancer of Unknown Primary (CUP) samples profiled with Visium HD
(16 µm bins). The pipeline goes from raw Space Ranger output through QC, tumour/stroma
compartment annotation via UniCell Deconvolve (UCD), and IHC marker classification with
tumour-of-origin prediction.

---

## Samples

- CUP_062
- CUP_138
- CUP_204_1O

---

## Prerequisites

### R packages

| Package | Purpose |
|---------|---------|
| `Seurat` (≥ 5.0) | Spatial data loading, subsetting, normalisation |
| `Matrix` | Sparse matrix operations |
| `dplyr` / `tidyr` / `tibble` / `stringr` | Data manipulation |
| `mclust` | GMM-based expression thresholds (script 04) |

### Python environment (UCD)

Requires the `ucdeconvolve` package and a UCD API token:

```bash
pip install ucdeconvolve anndata scanpy pandas
export UCD_TOKEN="your_token_here"
```

Register at <https://celldeconvolve.com> to obtain a token.

---

## Usage

Run scripts in order. Each step reads the output from the previous one.

```bash
# Step 1 — QC filtering and normalisation
Rscript 01_QC.R

# Step 2a — Export UCD input (full tissue)
# Run Stage 1 only (comment out Stage 2 loop)
Rscript 02_tumor_annotation.R

# Step 2b — Run full-tissue UCD deconvolution
export UCD_TOKEN="your_token_here"
python 03_UCD.py  # point to UCD_input_{sample}/ folders

# Step 2c — Annotate bins using UCD output
# Run Stage 2 (uncomment Stage 2 loop)
Rscript 02_tumor_annotation.R

# Step 3 — Export tumour-only UCD input
Rscript 03_UCD.R

# Step 4 — Run tumour-only UCD deconvolution
python 03_UCD.py  # point to UCD_input_{sample}_TumorOnly/ folders

# Step 5 — Marker classification and tumour-origin prediction
Rscript 04_marker_classification.R
```

### Expected directory layout

```
0.visiumHD/
├── objects/
│   ├── L02_16um_seu_qc.rds
│   ├── L03_16um_seu_qc.rds
│   └── L04_1O_16um_seu_qc.rds
└── 5.UCD/
    ├── UCD_input_{sample}/           # full-tissue UCD input (Stage 1)
    ├── UCD_output_{sample}/          # full-tissue UCD output
    ├── UCD_annotated_objects_16um/   # annotated Seurat objects
    ├── UCD_input_{sample}_TumorOnly/ # tumour-only UCD input (Step 3)
    └── UCD_output_{sample}_TumorOnly/# tumour-only UCD output
```

### QC parameters

| Parameter | Value |
|-----------|-------|
| Bin size | 16 µm |
| MT threshold | < 15 % |
| Min genes/bin | Sample-specific (see `01_QC.R`) |
| Upper UMI outlier | 99.5th percentile |
| Min gene UMI | ≥ 5 total across all bins |
| Tumour ratio cutoff | ≥ 0.60 |
