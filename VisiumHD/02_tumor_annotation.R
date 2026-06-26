# Visium HD 16 µm — Tumor / Stroma Annotation via UCD
# Stage 1: export counts + coordinates → run 03_UCD.py
# Stage 2: load UCD output, assign compartments, save

library(Seurat)
library(Matrix)
library(dplyr)
library(stringr)

# Parameters
assay_name <- "Spatial.016um"
image_name <- "slice1.016um"

# Bins with UCD tumor-like ratio >= 0.60 are classified as Tumor
tumor_ratio_threshold <- 0.60

# Sample paths
samples <- list(

  CUP_062 = list(
    sample_id         = "L02",
    input_rds         = "./0.visiumHD/objects/L02_16um_seu_qc.rds",
    ucd_input_dir     = "./0.visiumHD/5.UCD/UCD_input_L02_16um",
    ucd_output_dir    = "./0.visiumHD/5.UCD/UCD_output_L02_16um",
    annotated_out     = "./0.visiumHD/5.UCD/UCD_annotated_objects_16um/L02_16um_seu_ucd_annotated.rds"
  ),

  CUP_138 = list(
    sample_id         = "L03",
    input_rds         = "./0.visiumHD/objects/L03_16um_seu_qc.rds",
    ucd_input_dir     = "./0.visiumHD/5.UCD/UCD_input_L03_16um",
    ucd_output_dir    = "./0.visiumHD/5.UCD/UCD_output_L03_16um",
    annotated_out     = "./0.visiumHD/5.UCD/UCD_annotated_objects_16um/L03_16um_seu_ucd_annotated.rds"
  ),

  CUP_204_1O = list(
    sample_id         = "L04_1O",
    input_rds         = "./0.visiumHD/objects/L04_1O_16um_seu_qc.rds",
    ucd_input_dir     = "./0.visiumHD/5.UCD/UCD_input_L04_1O_16um",
    ucd_output_dir    = "./0.visiumHD/5.UCD/UCD_output_L04_1O_16um",
    annotated_out     = "./0.visiumHD/5.UCD/UCD_annotated_objects_16um/L04_1O_16um_seu_ucd_annotated.rds"
  )
)

# Keyword patterns to classify UCD cell-type columns (case-insensitive)

tumor_like_patterns <- c(
  "epithelial", "malignant", "carcinoma", "adenocarcinoma",
  "duct", "glandular", "secretory", "basal cell", "luminal",
  "urothelial", "cholangiocyte", "hepatocyte", "enterocyte",
  "goblet", "paneth", "intestinal", "stomach", "pancreatic ductal",
  "pneumocyte", "lung secretory", "lung goblet", "bronchial epithelial",
  "respiratory basal", "prostate epithelium", "ovarian surface epithelial",
  "mammary gland glandular", "colon epithelial"
)

cancer_type_patterns <- c(
  "luad", "lusc", "coad", "stad", "pdac", "chol", "lihc",
  "brca", "ov", "hgsoc", "ovca", "prad", "rcc", "blca",
  "skcm", "hnsc", "gbm", "uvm", "net"
)

cancer_line_patterns <- c(
  "ncih", "caov", "mcf", "hcc", "hct116", "sw480", "rko",
  "panc", "aspc", "pk", "snu", "ovcar", "igrov", "a549",
  "calu", "huh", "caki", "mkn", "t47d", "bt474", "mdamb",
  "pc3", "rt4", "tccsup", "ht1376"
)

stroma_patterns <- c(
  "fibroblast", "stromal", "mesenchymal", "myofibroblast",
  "smooth muscle", "pericyte", "mural cell", "adipocyte", "fat cell",
  "endothelial", "vascular", "stellate cell", "connective tissue"
)

immune_patterns <- c(
  "t cell", "b cell", "plasma", "plasmablast", "lymphocyte",
  "nk", "monocyte", "macrophage", "myeloid", "dendritic",
  "neutrophil", "mast", "eosinophil", "basophil", "leukocyte",
  "hematopoietic", "granulocyte"
)

# Helper: classify UCD columns
classify_ucd_columns <- function(celltypes) {
  all_tumor <- c(tumor_like_patterns, cancer_type_patterns, cancer_line_patterns)

  tumor_cols <- celltypes[str_detect(
    celltypes, regex(paste(all_tumor, collapse = "|"), ignore_case = TRUE))]

  stroma_cols <- celltypes[str_detect(
    celltypes, regex(paste(stroma_patterns, collapse = "|"), ignore_case = TRUE))]

  immune_cols <- celltypes[str_detect(
    celltypes, regex(paste(immune_patterns, collapse = "|"), ignore_case = TRUE))]

  non_tumor_cols <- unique(c(stroma_cols, immune_cols))

  # Remove any overlap (cell types matching both tumor and non-tumor) from both
  overlap        <- intersect(tumor_cols, non_tumor_cols)
  tumor_cols     <- setdiff(tumor_cols, overlap)
  non_tumor_cols <- setdiff(non_tumor_cols, overlap)

  list(tumor = tumor_cols, non_tumor = non_tumor_cols)
}

# STAGE 1 — Export UCD input (counts + coordinates)
# Run 03_UCD.py on each sample's input folder before Stage 2.

for (sample_name in names(samples)) {

  s <- samples[[sample_name]]
  message("\n=== STAGE 1 | ", sample_name, " ===")

  seu <- readRDS(s$input_rds)
  DefaultAssay(seu) <- assay_name

  dir.create(s$ucd_input_dir, showWarnings = FALSE, recursive = TRUE)

  # Counts matrix (genes x bins)
  counts <- GetAssayData(seu, assay = assay_name, layer = "counts")
  write.csv(as.matrix(counts), file.path(s$ucd_input_dir, "counts.csv"))

  # Spatial coordinates
  coords <- GetTissueCoordinates(seu[[image_name]])
  coords <- coords[colnames(seu), , drop = FALSE]
  coord_export <- data.frame(x = coords$x, y = coords$y, row.names = rownames(coords))
  write.csv(coord_export, file.path(s$ucd_input_dir, "coord.csv"))

  message(sprintf(
    "  Exported: %d genes x %d bins | coords: %d bins",
    nrow(counts), ncol(counts), nrow(coord_export)))
  message("  UCD input written to: ", s$ucd_input_dir)
}

message("\n>>> Run 03_UCD.py on each sample's input folder before continuing to Stage 2.")

# STAGE 2 — Load UCD output, annotate, and save
# Requires ucdbase_raw.csv in each sample's UCD output folder.

for (sample_name in names(samples)) {

  s <- samples[[sample_name]]
  message("\n=== STAGE 2 | ", sample_name, " ===")

  # Load QC Seurat object
  seu <- readRDS(s$input_rds)
  DefaultAssay(seu) <- assay_name

  # Load UCD raw deconvolution output
  ucd_raw_file <- file.path(s$ucd_output_dir, "ucdbase_raw.csv")
  if (!file.exists(ucd_raw_file)) {
    stop("UCD output not found: ", ucd_raw_file,
         "\nRun 03_UCD.py first.")
  }
  ucd_raw <- read.delim(ucd_raw_file, row.names = 1, check.names = FALSE)

  # Classify UCD cell-type columns as tumor-like vs non-tumor
  cols       <- classify_ucd_columns(colnames(ucd_raw))
  common_bins <- intersect(colnames(seu), rownames(ucd_raw))
  ucd_use    <- ucd_raw[common_bins, , drop = FALSE]

  # Compute tumor-like ratio per bin
  tumor_score <- if (length(cols$tumor) > 0)
    rowSums(ucd_use[, cols$tumor,     drop = FALSE], na.rm = TRUE) else 0
  non_tumor_score <- if (length(cols$non_tumor) > 0)
    rowSums(ucd_use[, cols$non_tumor, drop = FALSE], na.rm = TRUE) else 0

  total_score  <- tumor_score + non_tumor_score
  tumor_ratio  <- ifelse(total_score > 0, tumor_score / total_score, NA_real_)

  # Assign compartments: Tumor >= 0.60, else Stroma; NA → Unclassified
  compartment <- case_when(
    is.na(tumor_ratio)                    ~ "Unclassified",
    tumor_ratio >= tumor_ratio_threshold  ~ "Tumor",
    TRUE                                  ~ "Stroma"
  )

  # Add UCD-derived metadata to Seurat object
  seu$ucd_tumor_score     <- NA_real_
  seu$ucd_non_tumor_score <- NA_real_
  seu$ucd_tumor_ratio     <- NA_real_
  seu$ucd_compartment     <- NA_character_

  seu@meta.data[common_bins, "ucd_tumor_score"]     <- tumor_score
  seu@meta.data[common_bins, "ucd_non_tumor_score"] <- non_tumor_score
  seu@meta.data[common_bins, "ucd_tumor_ratio"]     <- tumor_ratio
  seu@meta.data[common_bins, "ucd_compartment"]     <- compartment

  message(table(seu$ucd_compartment))

  # Save annotated object
  dir.create(dirname(s$annotated_out), showWarnings = FALSE, recursive = TRUE)
  saveRDS(seu, file = s$annotated_out)
  message("  Saved: ", s$annotated_out)
}

message("\nTumor annotation complete for all samples.")
