# Visium HD 16 µm — Quality Control
# Samples: CUP_062 (L02), CUP_138 (L03), CUP_204_1O (L04_1O)

library(Seurat)
library(Matrix)

# Parameters
assay_name            <- "Spatial.016um"
image_name            <- "slice1.016um"
bin_size              <- 16
percent_mt_threshold  <- 15      # bins with > 15 % MT reads excluded
min_gene_umi          <- 5       # genes with < 5 total UMI excluded
ncount_upper_quantile <- 0.995   # upper UMI outlier removal threshold

# nFeature_threshold: sample-specific, set by inspecting retention curves during QC.
# roi_filter: coordinate boundaries to isolate each section (two samples share a slide).

samples <- list(

  CUP_062 = list(
    data_dir           = "/home/irsp_carles/phd/data/CUP/SpaceRanger/L02_ind/outs",
    sample_id          = "L02",
    nFeature_threshold = 10,
    # L02 shares the capture area with another section; keep only bins
    # with y_plot (= -y) strictly above -15 000 and at most 0.
    roi_filter = function(coords) {
      coords$y_plot > -15000 & coords$y_plot <= 0
    }
  ),

  CUP_138 = list(
    data_dir           = "/home/irsp_carles/phd/data/CUP/SpaceRanger/A1_095/outs",
    sample_id          = "L03",
    nFeature_threshold = 60,
    # L03 shares the capture area with another section; keep only bins
    # with y_plot at or below -19 000 and within the tissue y-range.
    roi_filter = function(coords) {
      coords$y_plot <= -19000 &
      coords$y_plot >= -25500 &
      coords$y_plot <=  19000
    }
  ),

  CUP_204_1O = list(
    data_dir           = "/home/irsp_carles/phd/data/CUP/SpaceRanger/D1_1125/outs",
    sample_id          = "L04_1O",
    nFeature_threshold = 5,
    # L04 capture area contains two sections (1O and 1P).
    # Primary split: bins with y_plot > -55 300 belong to 1O.
    # A boundary strip (x 8 000-19 500, y_plot -56 000 to -50 000)
    # is reassigned to 1P to clean the border between sections.
    roi_filter = function(coords) {
      region <- ifelse(coords$y_plot > -55300, "1O", "1P")
      boundary <- coords$x      >= 8000  & coords$x      <= 19500 &
                  coords$y_plot >= -56000 & coords$y_plot <= -50000
      region[boundary] <- "1P"
      region == "1O"
    }
  )
)

# Output directory
out_dir <- "./0.visiumHD/objects"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (sample_name in names(samples)) {

  s <- samples[[sample_name]]
  message("\n=== ", sample_name, " (", s$sample_id, ") ===")

  # 1. Load Space Ranger output
  seu <- Load10X_Spatial(
    data.dir      = s$data_dir,
    bin.size      = bin_size,
    filter.matrix = TRUE)

  DefaultAssay(seu) <- assay_name

  # 2. Spatial ROI subsetting
  coords           <- GetTissueCoordinates(seu[[image_name]])
  coords$y_plot    <- -coords[["y"]]
  coords$cell_id   <- rownames(coords)
  cells_keep       <- coords$cell_id[s$roi_filter(coords)]
  seu              <- subset(seu, cells = cells_keep)

  # 3. QC metrics: % mitochondrial reads
  gene_names  <- rownames(seu)
  mt_pattern  <- if      (any(grepl("^MT-", gene_names))) "^MT-"
                 else if (any(grepl("^mt-", gene_names))) "^mt-"
                 else NULL

  if (!is.null(mt_pattern)) {
    seu[["percent.mt"]] <- PercentageFeatureSet(
      seu, pattern = mt_pattern, assay = assay_name)
  } else {
    seu[["percent.mt"]] <- 0
  }

  ncount_col   <- paste0("nCount_",   assay_name)
  nfeature_col <- paste0("nFeature_", assay_name)

  # 4. Remove upper UMI outliers
  ncount_upper <- quantile(
    seu@meta.data[[ncount_col]], ncount_upper_quantile, na.rm = TRUE)
  seu <- subset(seu, cells = colnames(seu)[
    seu@meta.data[[ncount_col]] <= ncount_upper])

  # 5. Filter bins: sample-specific nFeature threshold + MT < 15 %
  keep_bins <- colnames(seu)[
    seu@meta.data[[nfeature_col]] >  s$nFeature_threshold &
    seu@meta.data[["percent.mt"]] <  percent_mt_threshold]
  seu <- subset(seu, cells = keep_bins)

  # 6. Filter genes: fewer than 5 total UMI across all bins
  counts_mat <- GetAssayData(seu, assay = assay_name, layer = "counts")
  seu        <- seu[rowSums(counts_mat) >= min_gene_umi, ]

  # 7. Normalize: LogNormalize, scale factor 10 000
  seu <- NormalizeData(
    seu,
    assay                = assay_name,
    normalization.method = "LogNormalize",
    scale.factor         = 10000,
    verbose              = FALSE)

  message(sprintf(
    "  Bins retained: %d | Genes retained: %d",
    ncol(seu), nrow(seu)))

  # 8. Save
  out_file <- file.path(out_dir, paste0(s$sample_id, "_16um_seu_qc.rds"))
  saveRDS(seu, file = out_file)
  message("  Saved: ", out_file)
}

message("\nQC complete for all samples.")
