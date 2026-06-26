#!/usr/bin/env python3
# =============================================================================
# 10b_ucdbase_singlecell.py
# -----------------------------------------------------------------------------
# PURPOSE : Run UCDBase deep-learning deconvolution at single-cell resolution
#           using a representative sketch of 10k cells per sample.
#
#   Cells are selected from the pre-computed sketch (UMAP_1 not NA in
#   cohort_slim_v2.rds, ~10% of each sample), then further subsampled to
#   MAX_SKETCH_CELLS (default 10,000) so every sample uses the same scale.
#   QC filters (nCount, nFeature) are applied before subsampling.
#
# METHOD  : UniCell Deconvolve Base (UCDBase) — Charytonowicz et al.,
#           Nature Communications 2023. Pre-trained on 28M annotated single
#           cells spanning 840 cell types. No reference data needed.
#
# INPUTS  : per_sample/<SAMPLE>/input/<SAMPLE>.h5ad   (cells × genes, raw counts)
#           cohort_slim_v2.rds → exported as cohort_sketch_metadata.csv at startup
#
# OUTPUTS : ImmuneInfiltration/UCDBase/singlecell/
#   ├── cohort_sketch_metadata.csv                     (R export, cached)
#   ├── per_cell/
#   │   └── <SAMPLE>_immune_fractions.csv              (per sketch-cell, 8 populations)
#   ├── adata_cache/
#   │   └── <SAMPLE>_ucdbase.h5ad                      (full results, checkpoint)
#   ├── ucdbase_sc_immune_fractions.csv                (sample-level aggregated)
#   ├── ucdbase_sc_immune_fractions_norm.csv
#   └── PLOTs/10_ucdbase/spatial/
#       └── <SAMPLE>_spatial_immune.png                (compartment + 8 immune panels)
#
# USAGE   : conda run -n ucdenv python 10b_ucdbase_singlecell.py --token $UCD_TOKEN
# =============================================================================

import os
import sys
import argparse
import logging
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

# ── 0. Argument parsing ───────────────────────────────────────────────────────

parser = argparse.ArgumentParser()
parser.add_argument("--token",         type=str,  default=os.environ.get("UCD_TOKEN", ""))
parser.add_argument("--analysis_root", type=str,  default=".")
parser.add_argument("--samples",       nargs="+", default=None)
parser.add_argument("--force",         action="store_true")
parser.add_argument("--max_sketch",    type=int,  default=10000,
                    help="Max sketch cells per sample after QC (default: 10000)")
parser.add_argument("--min_counts",    type=int,  default=50)
parser.add_argument("--max_counts",    type=int,  default=4000)
parser.add_argument("--min_features",  type=int,  default=20)
parser.add_argument("--seed",          type=int,  default=42)
args = parser.parse_args()

ANALYSIS_ROOT = args.analysis_root
OUT_DIR       = os.path.join(ANALYSIS_ROOT, "ImmuneInfiltration", "UCDBase", "singlecell")
PERCELL_DIR   = os.path.join(OUT_DIR, "per_cell")
CACHE_DIR     = os.path.join(OUT_DIR, "adata_cache")
PLOT_DIR      = os.path.join(ANALYSIS_ROOT, "PLOTs", "10_ucdbase", "spatial")
META_CSV      = os.path.join(OUT_DIR, "cohort_sketch_metadata.csv")
COHORT_RDS    = os.path.join(ANALYSIS_ROOT, "cohort_slim_v2.rds")

for d in [OUT_DIR, PERCELL_DIR, CACHE_DIR, PLOT_DIR,
          os.path.join(ANALYSIS_ROOT, "logs")]:
    os.makedirs(d, exist_ok=True)

# ── 1. Logging ────────────────────────────────────────────────────────────────

log_path = os.path.join(OUT_DIR, "ucdbase_singlecell.log")
logging.basicConfig(
    level    = logging.INFO,
    format   = "[%(asctime)s] %(message)s",
    datefmt  = "%H:%M:%S",
    handlers = [
        logging.FileHandler(log_path),
        logging.StreamHandler(sys.stdout)
    ]
)
log = logging.getLogger()
log.info("Script 10b — UCDBase Sketch Single-Cell Deconvolution — started")
log.info(f"Max sketch cells per sample: {args.max_sketch}")

if not args.token:
    log.error("No API token provided. Use --token or set UCD_TOKEN.")
    sys.exit(1)

# ── 2. Load UCDBase ───────────────────────────────────────────────────────────

try:
    import ucdeconvolve as ucd
    log.info(f"ucdeconvolve version: {ucd.__version__}")
except ImportError:
    log.error("ucdeconvolve not installed.")
    sys.exit(1)

log.info("Authenticating...")
ucd.api.authenticate(args.token)
log.info("Authentication successful.")

# ── 3. Export sketch metadata from R (cached) ─────────────────────────────────
# Sketch cells: UMAP_1 is not NA in cohort_slim_v2.rds

if not os.path.exists(META_CSV):
    log.info("Exporting sketch metadata from R...")
    r_script = f"""
suppressPackageStartupMessages(library(dplyr))
rds <- readRDS("{COHORT_RDS}")
df  <- rds |>
  filter(!is.na(UMAP_1)) |>
  dplyr::select(cell_id, sample_name, Compartment_binary, x, y)
write.csv(df, "{META_CSV}", row.names = FALSE)
cat(sprintf("Exported %d sketch cells from %d samples\\n",
    nrow(df), length(unique(df$sample_name))))
"""
    result = subprocess.run(["Rscript", "--vanilla", "-e", r_script],
                            capture_output=True, text=True)
    if result.returncode != 0:
        log.error(f"R export failed:\n{result.stderr}")
        sys.exit(1)
    log.info(result.stdout.strip())
else:
    log.info(f"Using cached sketch metadata: {META_CSV}")

log.info("Loading sketch metadata...")
meta_all = pd.read_csv(META_CSV)
log.info(f"  {len(meta_all):,} sketch cells | {meta_all['sample_name'].nunique()} samples")
log.info("  Sketch cells per sample:")
for s, n in meta_all.groupby("sample_name").size().sort_index().items():
    log.info(f"    {s}: {n:,}")

all_samples    = sorted(meta_all["sample_name"].unique().tolist())
samples_to_run = args.samples if args.samples else all_samples

# ── 4. Immune population keywords ─────────────────────────────────────────────

IMMUNE_KEYWORDS = {
    "T_CD8":          ["cd8", "cytotoxic t", "t cell cd8"],
    "T_CD4_Treg":     ["cd4", "t cell cd4", "treg", "regulatory t",
                       "helper t", "th1", "th2", "th17"],
    "B_cells":        ["b cell", "b-cell", "naive b", "memory b",
                       "germinal", "marginal zone b"],
    "Plasma_cells":   ["plasma cell", "plasmablast", "plasma b"],
    "DC":             ["dendritic", "cdc", "pdc"],
    "Macrophage_TAM": ["macrophage", "monocyte", "kupffer", "myeloid leukocyte",
                       "myeloid cell", "tissue-resident macrophage"],
    "NK_cells":       ["natural killer", "nk cell", "nk-cell"],
    "MDSC_Neutro":    ["mdsc", "neutrophil", "granulocyte",
                       "polymorphonuclear", "myeloid suppressor"],
}

IMMUNE_COLORS = {
    "T_CD8":          "#E63946",
    "T_CD4_Treg":     "#FF9F1C",
    "B_cells":        "#2EC4B6",
    "Plasma_cells":   "#4CC9F0",
    "DC":             "#F72585",
    "Macrophage_TAM": "#7B2FBE",
    "NK_cells":       "#06D6A0",
    "MDSC_Neutro":    "#FFB703",
}

COMP_COLORS = {"Tumor": "#FF6B6B", "Stroma": "#4CC9F0"}

# ── 5. Helper functions ───────────────────────────────────────────────────────

def resolve_colnames(adata, key="ucdbase_primary"):
    stem = key.split("ucdbase_")[-1]
    headers = adata.uns.get("ucdbase", {}).get("headers", {})
    if stem in headers:
        return list(headers[stem])
    n = adata.obsm[key].shape[1]
    log.warning(f"Cell type names not found for {key} — using generic names")
    return [f"ct_{i}" for i in range(n)]


def match_immune_populations(colnames):
    matched = {}
    for pop, kws in IMMUNE_KEYWORDS.items():
        hits = [c for c in colnames if any(kw in c.lower() for kw in kws)]
        matched[pop] = hits
        log.info(f"  {pop}: {len(hits)} cell types matched")
    return matched


# ── 6. Per-sample UCDBase loop ────────────────────────────────────────────────

log.info("=" * 60)
rng = np.random.default_rng(args.seed)

for sname in samples_to_run:
    h5ad_path  = os.path.join(ANALYSIS_ROOT, "per_sample", sname, "input",
                              f"{sname}.h5ad")
    cache_h5ad = os.path.join(CACHE_DIR, f"{sname}_ucdbase.h5ad")
    out_csv    = os.path.join(PERCELL_DIR, f"{sname}_immune_fractions.csv")

    log.info(f"\n{'─' * 50}")
    log.info(f"Sample: {sname}")

    if not os.path.exists(h5ad_path):
        log.warning(f"  [SKIP] h5ad not found: {h5ad_path}")
        continue

    if os.path.exists(out_csv) and not args.force:
        log.info(f"  [CACHED] — skipping API call")
        continue

    # ── Select sketch cells for this sample ───────────────────────────────────
    meta_s = meta_all[meta_all["sample_name"] == sname].copy()
    log.info(f"  Sketch pool: {len(meta_s):,} cells")

    # ── Load h5ad and apply QC filter on sketch cells ─────────────────────────
    adata_full = sc.read_h5ad(h5ad_path)
    sketch_ids = meta_s["cell_id"].values
    shared     = np.intersect1d(sketch_ids, adata_full.obs_names)
    adata      = adata_full[shared].copy()
    del adata_full

    n_before = adata.n_obs
    qc_pass = (
        (adata.obs["nCount_Xenium"]   >= args.min_counts)  &
        (adata.obs["nCount_Xenium"]   <= args.max_counts)  &
        (adata.obs["nFeature_Xenium"] >= args.min_features)
    )
    if "nCount_BlankCodeword" in adata.obs.columns:
        blank_rate = (adata.obs["nCount_BlankCodeword"] /
                      (adata.obs["nCount_Xenium"] + 1e-9))
        qc_pass = qc_pass & (blank_rate < 0.01)
    adata = adata[qc_pass].copy()
    log.info(f"  After QC: {adata.n_obs:,} / {n_before:,} cells "
             f"({100 * adata.n_obs / n_before:.1f}% kept)")

    # ── Subsample to MAX_SKETCH_CELLS ─────────────────────────────────────────
    if adata.n_obs > args.max_sketch:
        idx = rng.choice(adata.n_obs, size=args.max_sketch, replace=False)
        idx = np.sort(idx)
        adata = adata[idx].copy()
        log.info(f"  Subsampled to {adata.n_obs:,} cells")
    else:
        log.info(f"  Using all {adata.n_obs:,} cells (below max_sketch)")

    # Attach spatial / compartment metadata
    meta_idx = meta_s.set_index("cell_id")
    shared2  = np.intersect1d(adata.obs_names, meta_idx.index)
    adata    = adata[shared2].copy()
    adata.obs["Compartment_binary"] = meta_idx.loc[shared2, "Compartment_binary"].values
    adata.obs["x"]                  = meta_idx.loc[shared2, "x"].values
    adata.obs["y"]                  = meta_idx.loc[shared2, "y"].values

    # ── Normalise ─────────────────────────────────────────────────────────────
    adata.raw = adata
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    log.info(f"  Normalised {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # ── Run UCDBase ───────────────────────────────────────────────────────────
    log.info(f"  Running ucd.tl.base()...")
    try:
        ucd.tl.base(adata, split=True, sort=True)
        log.info("  UCDBase complete.")
    except Exception as e:
        log.error(f"  UCDBase failed: {e}")
        continue

    # ── Cache ──────────────────────────────────────────────────────────────────
    try:
        adata.write_h5ad(cache_h5ad)
        log.info(f"  Cached: {cache_h5ad}")
    except Exception as e:
        log.warning(f"  Cache write failed: {e}")

    # ── Extract immune fractions ──────────────────────────────────────────────
    if "ucdbase_primary" not in adata.obsm:
        log.error("  ucdbase_primary not in obsm — skipping")
        continue

    colnames = resolve_colnames(adata, "ucdbase_primary")
    matched  = match_immune_populations(colnames)

    primary_df = pd.DataFrame(
        adata.obsm["ucdbase_primary"],
        index   = adata.obs_names,
        columns = colnames
    )

    immune_df = pd.DataFrame(index=adata.obs_names)
    for pop, cts in matched.items():
        immune_df[pop] = primary_df[cts].sum(axis=1) if cts else 0.0

    immune_df["Compartment_binary"] = adata.obs["Compartment_binary"].values
    immune_df["x"]                  = adata.obs["x"].values
    immune_df["y"]                  = adata.obs["y"].values
    immune_df["sample_name"]        = sname

    immune_df.to_csv(out_csv)
    log.info(f"  Per-cell fractions saved: {out_csv}")


# ── 7. Aggregate to sample level ─────────────────────────────────────────────

log.info("\n" + "=" * 60)
log.info("Aggregating to sample level")

pop_names = list(IMMUNE_KEYWORDS.keys())
agg_rows  = []

for sname in samples_to_run:
    csv_path = os.path.join(PERCELL_DIR, f"{sname}_immune_fractions.csv")
    if not os.path.exists(csv_path):
        continue
    df = pd.read_csv(csv_path, index_col=0)
    avail = [p for p in pop_names if p in df.columns]
    row = {"sample_name": sname}
    for pop in avail:
        row[pop] = df[pop].mean()
    agg_rows.append(row)

if not agg_rows:
    log.warning("No results to aggregate — check API token and re-run.")
    sys.exit(0)

agg_df = pd.DataFrame(agg_rows).set_index("sample_name")
agg_df.to_csv(os.path.join(OUT_DIR, "ucdbase_sc_immune_fractions.csv"))
agg_norm = agg_df.div(agg_df.sum(axis=1), axis=0).fillna(0)
agg_norm.to_csv(os.path.join(OUT_DIR, "ucdbase_sc_immune_fractions_norm.csv"))

log.info(f"Aggregated {len(agg_df)} samples")
log.info("\nMean immune fractions:")
log.info(agg_df.mean().sort_values(ascending=False).round(4).to_string())


# ── 8. Summary plots ──────────────────────────────────────────────────────────

# Heatmap
fig, ax = plt.subplots(figsize=(max(10, len(agg_df) * 0.65), 6))
sns.heatmap(
    agg_df.T, cmap="YlOrRd", linewidths=0.5, linecolor="white",
    annot=True, fmt=".3f", annot_kws={"size": 7}, ax=ax,
    cbar_kws={"label": "Mean predicted fraction", "shrink": 0.6}
)
ax.set_title(f"UCDBase Single-Cell Immune Deconvolution\n"
             f"(sketch {args.max_sketch:,} cells/sample, mean fraction)",
             fontsize=12, fontweight="bold", pad=10)
ax.tick_params(axis="x", rotation=45, labelsize=8)
ax.tick_params(axis="y", rotation=0,  labelsize=8)
plt.tight_layout()
hm_path = os.path.join(ANALYSIS_ROOT, "PLOTs", "10_ucdbase",
                       "ucdbase_sc_immune_heatmap.png")
plt.savefig(hm_path, dpi=200, bbox_inches="tight", facecolor="white")
plt.close()
log.info(f"Heatmap saved: {hm_path}")

# Stacked bar (normalised)
fig, ax = plt.subplots(figsize=(max(12, len(agg_norm) * 0.75), 6))
order      = agg_norm.sort_values("Macrophage_TAM").index.tolist()
agg_sorted = agg_norm.loc[order]
bottom     = np.zeros(len(order))
x          = np.arange(len(order))
for pop in agg_sorted.columns:
    ax.bar(x, agg_sorted[pop].values, bottom=bottom,
           color=IMMUNE_COLORS.get(pop, "grey"),
           label=pop.replace("_", " "),
           width=0.75, edgecolor="white", linewidth=0.5)
    bottom += agg_sorted[pop].values
ax.set_xticks(x)
ax.set_xticklabels(order, rotation=45, ha="right", fontsize=9, fontweight="bold")
ax.set_ylabel("Normalised immune fraction", fontsize=10)
ax.set_title(f"UCDBase Single-Cell Immune Composition  "
             f"(sketch {args.max_sketch:,} cells/sample, sorted by TAM fraction)",
             fontsize=11, fontweight="bold")
ax.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=8, frameon=False)
ax.set_ylim(0, 1.05)
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
bar_path = os.path.join(ANALYSIS_ROOT, "PLOTs", "10_ucdbase",
                        "ucdbase_sc_immune_barplot.png")
plt.savefig(bar_path, dpi=200, bbox_inches="tight", facecolor="white")
plt.close()
log.info(f"Bar chart saved: {bar_path}")


# ── 9. Per-sample spatial plots ───────────────────────────────────────────────

log.info("\n" + "=" * 60)
log.info("Generating per-sample spatial immune maps")

for sname in samples_to_run:
    csv_path = os.path.join(PERCELL_DIR, f"{sname}_immune_fractions.csv")
    if not os.path.exists(csv_path):
        continue

    df = pd.read_csv(csv_path, index_col=0).dropna(subset=["x", "y"])
    if len(df) == 0:
        continue

    avail_pops = [p for p in pop_names if p in df.columns]
    n_panels   = 1 + len(avail_pops)
    ncols      = 3
    nrows      = int(np.ceil(n_panels / ncols))

    fig = plt.figure(figsize=(ncols * 5, nrows * 5), facecolor="black")
    gs  = gridspec.GridSpec(nrows, ncols, figure=fig, hspace=0.35, wspace=0.25)

    # Adaptive point size — sketch cells are sparse so use larger points
    pt = max(0.5, min(4.0, 100000 / len(df)))

    # Panel 0: compartment
    ax = fig.add_subplot(gs[0, 0])
    ax.set_facecolor("black")
    for comp, col in COMP_COLORS.items():
        mask = df["Compartment_binary"] == comp
        ax.scatter(df.loc[mask, "x"], df.loc[mask, "y"],
                   c=col, s=pt, alpha=0.8, linewidths=0, rasterized=True,
                   label=comp)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Tumour / Stroma", color="white", fontsize=9, fontweight="bold")
    ax.legend(loc="lower right", frameon=False, labelcolor="white",
              fontsize=7, markerscale=4, handletextpad=0.3)

    # Immune population panels
    for idx, pop in enumerate(avail_pops):
        r = (idx + 1) // ncols
        c = (idx + 1) %  ncols
        ax = fig.add_subplot(gs[r, c])
        ax.set_facecolor("black")

        vals     = df[pop].values
        pos_vals = vals[vals > 0]
        q99      = np.quantile(pos_vals, 0.99) if len(pos_vals) > 0 else 1e-6

        sc_plot = ax.scatter(
            df["x"], df["y"],
            c=np.clip(vals, 0, q99), s=pt, alpha=0.85,
            linewidths=0, cmap="magma", vmin=0, vmax=q99, rasterized=True
        )
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(pop.replace("_", " "), color="white",
                     fontsize=9, fontweight="bold")

        cbar = plt.colorbar(sc_plot, ax=ax, shrink=0.6, pad=0.02,
                            orientation="horizontal")
        cbar.ax.tick_params(labelsize=6, colors="white")
        cbar.outline.set_edgecolor("white")
        cbar.set_label("Predicted fraction", color="white", fontsize=6)

    # Hide unused panels
    for idx in range(n_panels, nrows * ncols):
        fig.add_subplot(gs[idx // ncols, idx % ncols]).axis("off")

    fig.suptitle(f"{sname} — UCDBase Immune Deconvolution  "
                 f"(sketch {len(df):,} cells)",
                 color="white", fontsize=11, fontweight="bold", y=1.01)

    out_png = os.path.join(PLOT_DIR, f"{sname}_spatial_immune.png")
    plt.savefig(out_png, dpi=150, bbox_inches="tight",
                facecolor="black", edgecolor="none")
    plt.close()
    log.info(f"  {sname}: spatial map saved ({len(df):,} cells plotted)")


# ── Final summary ─────────────────────────────────────────────────────────────

log.info("=" * 60)
log.info("SCRIPT 10b COMPLETE")
log.info(f"Outputs: {OUT_DIR}")
log.info("=" * 60)
