# CUP Visium Spatial Transcriptomics Pipeline

Analysis pipeline for 20 Cancer of Unknown Primary (CUP) Visium spatial transcriptomics
samples. The pipeline proceeds from raw SpaceRanger output through QC, spatial enhancement,
tumour/stroma classification, immune deconvolution, pathway scoring, niche definition, and
ligand–receptor interaction analysis.

---

## Scripts

| Script | Purpose | Key Inputs | Key Outputs |
|--------|---------|------------|-------------|
| `V01_qc_normalisation.qmd` | Per-sample QC (UMI, gene count, mitochondrial %) filtering and SCTransform normalisation. Produces a merged QC UMAP for inspection. | SpaceRanger output directories (`<sample>/filtered_feature_bc_matrix.h5`, `<sample>/spatial/`) | `results/01_qc/visium_seurat_list.rds` |
| `V02_bayesspace_enhancement.qmd` | BayesSpace spatial domain clustering and sub-spot enhancement (`nsubspots.per.edge = 4`, ~25 µm effective spacing). Runs all 20 samples in parallel via `mclapply`; per-sample RDS caches written to `bayesspace_cache/`. | `results/01_qc/visium_seurat_list.rds` | `visium_enhanced_sce_list.rds`, `bayesspace_cache/<sample>_enhanced.rds` |
| `V03_cancerfinder3.qmd` | Tumour/Stroma compartment annotation on BayesSpace sub-spots using CancerFinder3 (threshold: cancer_prob ≥ 0.85). Exports pseudo-count MTX for each sample, runs CF3 via `conda run`, and attaches `compartment` labels to each sub-spot. | `visium_enhanced_sce_list.rds`, CancerFinder3 model + infer.py (external) | `visium_annotated_sce_list.rds`, `cf3_results/<sample>_predictions.csv` |
| `V04_cell2location.qmd` | Fine-grained immune deconvolution (10 purified PBMC cell types) using cell2location. Builds a `RegressionModel` reference from `scVI purified_pbmc_dataset` then runs `Cell2location` per sample. Tumour sub-spots receive immune infiltration scores; Stroma sub-spots are labelled by dominant immune type. | `visium_annotated_sce_list.rds` | `visium_celltypes_sce_list.rds`, `c2l_results/<sample>_abundances.csv` |
| `V04_ucdbase.qmd` | Coarse immune deconvolution using UCD_Base (UniCell Deconvolve cloud API, no reference required). Maps 840+ cell types to 7 broad immune categories. Provides a complementary annotation to V04_cell2location. | `visium_annotated_sce_list.rds` | `visium_celltypes_ucd_sce_list.rds`, `ucd_results/<sample>_fractions.csv` |
| `V05_progeny_classification.qmd` | Pseudobulk PROGENy pathway scoring (7 pathways, 500 top genes). ONC score = mean(PI3K, EGFR, MAPK, VEGF) from Tumour pseudobulk; IMI score = mean(JAK-STAT, NFkB, TNFa) from Stroma pseudobulk. Samples classified ONC (δ > 0) or IMI (δ ≤ 0). | `visium_celltypes_sce_list.rds` | `visium_progeny_classification.csv` |
| `V06_niche_definition.qmd` | Assigns each sub-spot to a spatial niche zone: **Tumour**, **Peritumoral** (< 100 µm from nearest Tumour sub-spot), or **Distal** (≥ 100 µm). Uses pixel coordinates converted to µm via per-sample `scalefactors_json.json`. Performs ONC/IMI-stratified immune enrichment analysis using `metadata/ONC_IMI_phenotypes.csv`. | `visium_celltypes_ucd_sce_list.rds`, `metadata/ONC_IMI_phenotypes.csv` | `visium_niches_ucd_sce_list.rds`, `Plots/` (PNG + XLSX) |
| `V07_cellphonedb.qmd` | Ligand–receptor interaction analysis (LIANA v1.7.1, CellPhoneDB method) restricted to the tumour–immune interface (Peritumoral sub-spots + Tumour sub-spots per sample). Produces per-sample and consolidated dot plots, lollipop plots of top ONC/IMI-differential interactions, and spatial LR maps. | `visium_niches_ucd_sce_list.rds`, `metadata/ONC_IMI_phenotypes.csv` | `cpdb_results/<sample>/liana_results.csv`, `cpdb_results/LIANA_Visium.csv`, `cpdb_ucd_results/` (PNG + XLSX) |

---

## Prerequisites / Dependencies

### R packages

| Package | Version tested | Purpose |
|---------|---------------|---------|
| `Seurat` | ≥ 5.0 | Spatial data loading, SCTransform normalisation |
| `BayesSpace` | ≥ 1.12 | Spatial domain clustering, sub-spot enhancement |
| `SingleCellExperiment` | ≥ 1.24 | SCE data structure throughout pipeline |
| `Matrix` | — | Sparse matrix operations |
| `parallel` | — | `mclapply` parallelisation (V02) |
| `RhpcBLASctl` | optional | Per-worker BLAS thread control (V02) |
| `ggplot2` | — | Plotting |
| `patchwork` | — | Multi-panel figure assembly |
| `dplyr` / `tidyr` / `tibble` | — | Data manipulation |
| `scales` | — | Axis formatting |
| `jsonlite` | — | Reading `scalefactors_json.json` (V06) |
| `RANN` | — | Fast nearest-neighbour search for niche distances (V06) |
| `progeny` | — | PROGENy pathway scoring (V05) |
| `decoupleR` | — | Footprint methods (imported alongside progeny in V05) |
| `pheatmap` | — | Pathway and LR heatmaps |
| `ggrepel` | — | Non-overlapping text labels (V05 scatter) |
| `openxlsx` | — | Excel export of plot data (V06, V07) |
| `ggnewscale` | — | Dual colour scales in spatial LR maps (V07) |

### Python environments (conda)

| Environment | Tools | Used in |
|-------------|-------|---------|
| `scf` | `anndata`, `scipy`, `pandas`, `numpy`, CancerFinder3 | V03 |
| `cell2loc_env` | `cell2location`, `scvi-tools`, `scanpy`, `torch` | V04 |
| `ucdenv` | `ucdeconvolve`, `anndata`, `scanpy` | V04b |
| `squidpy` | `liana`, `anndata`, `pandas`, `scipy` | V07 |

---

## Data

### Visium SpaceRanger output

Each sample directory must contain:

```
<sample_id>/
    filtered_feature_bc_matrix.h5
    spatial/
        tissue_positions.csv
        tissue_lowres_image.png
        scalefactors_json.json
```

Place all 20 sample directories at the repository root (i.e., next to the `.qmd` scripts).
Sample IDs used in the pipeline: `114B`, `271`, `274`, `297`, `351`, `353`, `355`, `368`,
`369`, `372`, `373`, `374`, `375`, `376`, `381`, `383`, `385`, `386`, `387`, `394`.

### Metadata file

`metadata/ONC_IMI_phenotypes.csv` — two-column CSV mapping sample IDs (column 1,
optionally prefixed `CUP-`) to ONC/IMI phenotype (column 2). Used in V06 and V07 for
ONC/IMI-stratified analyses.

### CancerFinder3 (external tool)

Download CancerFinder3 from <https://github.com/JiangBioLab/CancerFinder> and set the
following variables at the top of `V03_cancerfinder3.qmd`:

```r
CF3_MODEL <- "path/to/CancerFinder3/models/st_pretrain_article.pkl"
CF3_INFER <- "path/to/CancerFinder3/infer.py"
```

### UCD_Base API token

UCD_Base requires a cloud API token. Set it via the environment variable `UCD_TOKEN`
or directly in `V04_ucdbase.qmd`:

```r
UCD_TOKEN <- Sys.getenv("UCD_TOKEN", "your-ucd-token-here")
```

Register at <https://celldeconvolve.com> to obtain a token.

---

## Usage

Run scripts **in order** from the repository root. Each script reads checkpoint `.rds`
files written by the previous one. Render with Quarto:

```bash
# Step 1 — QC and normalisation
quarto render V01_qc_normalisation.qmd

# Step 2 — BayesSpace spatial enhancement (~60–90 min, 20 samples)
quarto render V02_bayesspace_enhancement.qmd

# Step 3 — CancerFinder3 tumour/stroma annotation (requires scf conda env)
quarto render V03_cancerfinder3.qmd

# Step 4a — Cell2location immune deconvolution (requires cell2loc_env conda env + GPU recommended)
quarto render V04_cell2location.qmd

# Step 4b — UCD_Base coarse immune deconvolution (requires ucdenv + UCD_TOKEN)
quarto render V04_ucdbase.qmd

# Step 5 — PROGENy ONC/IMI classification
quarto render V05_progeny_classification.qmd

# Step 6 — Spatial niche definition and immune enrichment
quarto render V06_niche_definition.qmd

# Step 7 — LIANA CellPhoneDB ligand-receptor analysis (requires squidpy conda env)
quarto render V07_cellphonedb.qmd
```

Each script writes a self-contained HTML report alongside the output `.rds` / `.csv` files.
Per-sample caches (BayesSpace in `bayesspace_cache/`, cell2location in `c2l_results/`,
UCD in `ucd_results/`, CancerFinder3 in `cf3_results/`) are skipped on re-runs, so it is
safe to re-render a partially completed pipeline.

### Expected directory layout after a full run

```
.
├── <sample_id>/          # SpaceRanger output per sample (20 total)
├── metadata/
│   └── ONC_IMI_phenotypes.csv
├── results/01_qc/
│   └── visium_seurat_list.rds
├── bayesspace_cache/
├── cf3_input/ / cf3_results/
├── c2l_results/
├── ucd_results/
├── cpdb_input/ / cpdb_results/ / cpdb_ucd_results/
├── Plots/
├── visium_enhanced_sce_list.rds
├── visium_annotated_sce_list.rds
├── visium_celltypes_sce_list.rds
├── visium_celltypes_ucd_sce_list.rds
├── visium_niches_ucd_sce_list.rds
└── visium_progeny_classification.csv
```
