# Visium HD 16 µm — Tumor-Only UCD Input Preparation
# Exports counts + coordinates for tumor bins only (tumor-of-origin input).
# Requires annotated objects from 02_tumor_annotation.R. Then run 03_UCD.py.

library(Seurat)
library(Matrix)

# Parameters
assay_name <- "Spatial.016um"
image_name <- "slice1.016um"

# Sample paths
samples <- list(

  CUP_062 = list(
    sample_id      = "L02",
    annotated_rds  = "./0.visiumHD/5.UCD/UCD_annotated_objects_16um/L02_16um_seu_ucd_annotated.rds",
    ucd_input_dir  = "./0.visiumHD/5.UCD/UCD_input_L02_16um_TumorOnly"
  ),

  CUP_138 = list(
    sample_id      = "L03",
    annotated_rds  = "./0.visiumHD/5.UCD/UCD_annotated_objects_16um/L03_16um_seu_ucd_annotated.rds",
    ucd_input_dir  = "./0.visiumHD/5.UCD/UCD_input_L03_16um_TumorOnly"
  ),

  CUP_204_1O = list(
    sample_id      = "L04_1O",
    annotated_rds  = "./0.visiumHD/5.UCD/UCD_annotated_objects_16um/L04_1O_16um_seu_ucd_annotated.rds",
    ucd_input_dir  = "./0.visiumHD/5.UCD/UCD_input_L04_1O_16um_TumorOnly"
  )
)

# Export tumor-only UCD inputs
for (sample_name in names(samples)) {

  s <- samples[[sample_name]]
  message("\n=== ", sample_name, " (", s$sample_id, ") ===")

  # Load annotated Seurat object
  seu <- readRDS(s$annotated_rds)
  DefaultAssay(seu) <- assay_name

  # Subset to Tumor bins only
  tumor_bins <- colnames(seu)[
    !is.na(seu$ucd_compartment) & seu$ucd_compartment == "Tumor"]

  if (length(tumor_bins) == 0) {
    warning("No Tumor bins found for ", sample_name, " — skipping.")
    next
  }

  seu_tumor <- subset(seu, cells = tumor_bins)

  # Counts matrix (genes x bins)
  counts <- GetAssayData(seu_tumor, assay = assay_name, layer = "counts")

  # Spatial coordinates
  coords <- GetTissueCoordinates(seu_tumor[[image_name]])
  coords <- coords[colnames(seu_tumor), , drop = FALSE]
  coord_export <- data.frame(x = coords$x, y = coords$y, row.names = rownames(coords))

  # Write outputs
  dir.create(s$ucd_input_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(as.matrix(counts),  file.path(s$ucd_input_dir, "counts.csv"))
  write.csv(coord_export,        file.path(s$ucd_input_dir, "coord.csv"))

  message(sprintf(
    "  Tumor bins: %d / %d total (%.1f%%)",
    length(tumor_bins), ncol(seu), 100 * length(tumor_bins) / ncol(seu)))
  message(sprintf(
    "  Exported: %d genes x %d bins", nrow(counts), ncol(counts)))
  message("  Written to: ", s$ucd_input_dir)
}

message("\nTumor-only UCD inputs ready. Run 03_UCD.py to execute UCD.")
