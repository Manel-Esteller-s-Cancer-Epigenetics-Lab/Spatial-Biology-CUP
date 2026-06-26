#!/usr/bin/env python3
# ============================================================
# Visium HD 16 µm — UniCell Deconvolve (UCD) Run
# ============================================================
# Runs UCD deconvolution on tumor-only spatial bins for all
# three CUP samples. Produces raw, primary, and cancer-level
# deconvolution scores used for tumor-origin prediction.
#
# Prerequisite: tumor-only counts.csv and coord.csv must exist
# in each sample's input directory. These are produced by
# 03_UCD.R.
#
# Authentication: set the UCD_TOKEN environment variable
# before running:
#   export UCD_TOKEN="your_token_here"
#   python 03_UCD.py
#
# Input:  UCD_input_*_TumorOnly/counts.csv + coord.csv
# Output: UCD_output_*_TumorOnly/ucdbase_{raw,primary,cancer}.csv
#         UCD_output_*_TumorOnly/adata_vis.h5ad
# ============================================================

import os
import sys
from pathlib import Path

import pandas as pd
from anndata import AnnData
import ucdeconvolve as ucd

# ── Authentication ────────────────────────────────────────────
token = os.environ.get("UCD_TOKEN")
if not token:
    sys.exit(
        "Error: UCD_TOKEN environment variable is not set.\n"
        "Set it with: export UCD_TOKEN='your_token_here'"
    )

# ── Paths ─────────────────────────────────────────────────────
base_dir = Path(__file__).resolve().parent

samples = {
    "CUP_062": {
        "sample_id":  "L02_16um",
        "input_dir":  base_dir / "UCD_input_L02_16um_TumorOnly",
        "output_dir": base_dir / "UCD_output_L02_16um_TumorOnly",
    },
    "CUP_138": {
        "sample_id":  "L03_16um",
        "input_dir":  base_dir / "UCD_input_L03_16um_TumorOnly",
        "output_dir": base_dir / "UCD_output_L03_16um_TumorOnly",
    },
    "CUP_204_1O": {
        "sample_id":  "L04_1O_16um",
        "input_dir":  base_dir / "UCD_input_L04_1O_16um_TumorOnly",
        "output_dir": base_dir / "UCD_output_L04_1O_16um_TumorOnly",
    },
}

# ── Authenticate once ─────────────────────────────────────────
ucd.api.authenticate(token)

# ── Process each sample ───────────────────────────────────────
for sample_name, s in samples.items():

    print(f"\n=== {sample_name} ({s['sample_id']}) ===")

    counts_path = s["input_dir"] / "counts.csv"
    coord_path  = s["input_dir"] / "coord.csv"

    if not counts_path.exists():
        raise FileNotFoundError(f"Counts file not found: {counts_path}")
    if not coord_path.exists():
        raise FileNotFoundError(f"Coordinate file not found: {coord_path}")

    # Load input files
    counts = pd.read_csv(counts_path, index_col=0)
    coord  = pd.read_csv(coord_path,  index_col=0)

    # Align barcodes
    common_bins = counts.columns.intersection(coord.index)
    if len(common_bins) == 0:
        raise ValueError(
            f"No matching barcodes between counts and coordinates for {sample_name}."
        )
    counts = counts.loc[:, common_bins]
    coord  = coord.loc[common_bins, :]

    print(f"  Counts : {counts.shape[0]} genes x {counts.shape[1]} bins")
    print(f"  Matched: {len(common_bins)} bins")

    # Build AnnData (bins x genes)
    adata = AnnData(counts.T)
    adata.var["SYMBOL"]      = adata.var_names
    adata.obs["array_row"]   = coord.loc[adata.obs_names, "y"].values
    adata.obs["array_col"]   = coord.loc[adata.obs_names, "x"].values
    adata.obsm["spatial"]    = coord.loc[adata.obs_names, ["x", "y"]].to_numpy()
    adata.obs["sample"]      = s["sample_id"]

    s["output_dir"].mkdir(parents=True, exist_ok=True)

    # Run UCD deconvolution
    try:
        ucd.tl.base(adata)

        for category in ["raw", "primary", "cancer"]:
            result = ucd.utils.read_results(adata, category=category)
            out_path = s["output_dir"] / f"ucdbase_{category}.csv"
            result.to_csv(out_path, sep="\t")
            print(f"  Saved {category}: {out_path}")

        adata.write(s["output_dir"] / "adata_vis.h5ad")
        print(f"  Saved AnnData: {s['output_dir'] / 'adata_vis.h5ad'}")
        print(f"  UCD complete for {sample_name}.")

    except Exception as e:
        print(f"\n  UCD failed for {sample_name}.")
        print(f"  {type(e).__name__}: {e}")
        try:
            adata.write(s["output_dir"] / "adata_vis_partial.h5ad")
            print("  Partial AnnData saved.")
        except Exception:
            pass

print("\nUCD complete for all samples.")
