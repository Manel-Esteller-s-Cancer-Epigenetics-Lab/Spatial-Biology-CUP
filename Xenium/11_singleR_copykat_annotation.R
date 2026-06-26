#!/usr/bin/env Rscript
# =============================================================================
# 11_singleR_copykat_annotation.R
# -----------------------------------------------------------------------------
# PURPOSE : Per-cell annotation for ALL cells using:
#
#   Seurat Compartment_binary — cancer vs non-cancer (Tumor / Stroma labels
#                               already computed from marker expression in
#                               Script 03; used directly as the cancer detector)
#
#   SingleR (HPCA)            — coarse cell type annotation on non-cancer cells
#                               Reference: HumanPrimaryCellAtlasData (celldex)
#                               Fine labels used to resolve CD4 vs CD8 T cells
#
# CLASSIFICATION:
#   Compartment_binary == "Tumor"  → Cancer cell  → label = sample tissue origin
#   Compartment_binary == "Stroma" → Non-cancer   → label = SingleR coarse type
#
# INPUTS  : per_sample/<S>/bpcells/counts/       BPCells on-disk count matrices
#           cohort_slim_v2.rds                   cell metadata + Compartment_binary
#           Classification/origin_classification_v2.csv
#
# OUTPUTS : CellAnnotation/
#   ├── per_cell/<S>_cell_annotation.csv
#   ├── cohort_cell_annotation.csv
#   ├── sample_cancer_fractions.csv
#   └── sample_celltype_fractions.csv
#   PLOTs/11_cellannotation/<S>_spatial_annotation.png
#
# USAGE   : Rscript 11_singleR_copykat_annotation.R
# =============================================================================

suppressPackageStartupMessages({
  library(BPCells)
  library(SingleR)
  library(celldex)
  library(BiocParallel)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ── 0. Configuration ──────────────────────────────────────────────────────────

ANALYSIS_ROOT <- "."
COHORT_RDS    <- file.path(ANALYSIS_ROOT, "cohort_slim_v2.rds")
ORIGIN_CSV    <- file.path(ANALYSIS_ROOT, "Classification", "origin_classification_v2.csv")

OUT_DIR     <- file.path(ANALYSIS_ROOT, "CellAnnotation")
PERCELL_DIR <- file.path(OUT_DIR, "per_cell")
PLOT_DIR    <- file.path(ANALYSIS_ROOT, "PLOTs", "11_cellannotation")

for (d in c(OUT_DIR, PERCELL_DIR, PLOT_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

SAMPLES_EXCLUDE <- c()
N_CORES         <- 18

# IC samples first, then CUP alphabetically
all_samples  <- setdiff(unique(readRDS(COHORT_RDS)$sample_name), SAMPLES_EXCLUDE)

# ── 1. Load cohort metadata ───────────────────────────────────────────────────

message("Loading cohort metadata...")
cohort_meta <- readRDS(COHORT_RDS) |>
  dplyr::select(cell_id, sample_name, x, y, Compartment_binary)

sample_names <- c(sort(grep("^IC_",  all_samples, value = TRUE)),
                  sort(grep("^IC_",  all_samples, value = TRUE, invert = TRUE)))

message(sprintf("  %d samples | %s cells",
                length(sample_names), format(nrow(cohort_meta), big.mark = ",")))
message("  Order: ", paste(sample_names, collapse = ", "))

# Sample-level tissue origin (strip score: "Lung_Adeno (237.6)" → "Lung_Adeno")
origin_df <- read.csv(ORIGIN_CSV) |>
  mutate(cancer_label = sub(" \\(.*", "", top1)) |>
  dplyr::select(sample_name, cancer_label)

# ── 2. Load SingleR reference ─────────────────────────────────────────────────

message("\nLoading HumanPrimaryCellAtlasData reference...")
hpca_ref <- celldex::HumanPrimaryCellAtlasData()

# Convert Ensembl rownames → gene symbols if needed (celldex >= 1.0)
if (grepl("^ENSG", rownames(hpca_ref)[1])) {
  message("  Converting Ensembl IDs to gene symbols...")
  sym  <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys      = rownames(hpca_ref),
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  keep <- !is.na(sym) & !duplicated(sym)
  hpca_ref <- hpca_ref[keep, ]
  rownames(hpca_ref) <- sym[keep]
}
message(sprintf("  %d genes | %d cell-type profiles",
                nrow(hpca_ref), ncol(hpca_ref)))

# ── 3. HPCA label → coarse display name ───────────────────────────────────────

HPCA_COARSE <- c(
  "B_cell"               = "B_cell",
  "T_cells"              = "T_cell",        # overridden per CD4/CD8 below
  "NK_cell"              = "NK",
  "DC"                   = "DC",
  "Monocyte"             = "Monocyte",
  "Macrophage"           = "Macrophage_TAM",
  "Neutrophils"          = "Neutrophil_MDSC",
  "Myelocyte"            = "Neutrophil_MDSC",
  "Pro-Myelocyte"        = "Neutrophil_MDSC",
  "Plasma_cells"         = "Plasma",
  "Fibroblasts"          = "Fibroblast_CAF",
  "MSC"                  = "Fibroblast_CAF",
  "Pre-adipocyte"        = "Fibroblast_CAF",
  "Tissue_stem_cells"    = "Fibroblast_CAF",
  "Endothelial_cells"    = "Endothelial",
  "Smooth_muscle_cells"  = "Pericyte_SMC",
  "Epithelial_cells"     = "Epithelial",
  "Keratinocytes"        = "Epithelial",
  "Hepatocytes"          = "Epithelial",
  "Neuroepithelial_cell" = "Neural",
  "Neurons"              = "Neural",
  "Astrocyte"            = "Neural",
  "Chondrocytes"         = "Other",
  "Osteoblasts"          = "Other",
  "MEP"                  = "Other",
  "CMP"                  = "Other",
  "GMP"                  = "Other",
  "HSC_-G-CSF"           = "Other",
  "HSC_CD34+"            = "Other",
  "Pro-B_cell_CD34-"     = "B_cell",
  "Pro-B_cell_CD34+"     = "B_cell",
  "Pre-B_cell_CD34-"     = "B_cell",
  "BM"                   = "Other",
  "BM & Prog."           = "Other",
  "Platelets"            = "Other",
  "Erythroblast"         = "Other",
  "Gametocytes"          = "Other",
  "iPS_cells"            = "Other",
  "Embryonic_stem_cells" = "Other"
)

# Fine T-cell label → CD4 / CD8 / generic T cell
tcell_fine_coarse <- function(fine) {
  l <- tolower(fine)
  dplyr::case_when(
    grepl("cd8",                  l) ~ "T_CD8",
    grepl("cd4|regulatory|treg",  l) ~ "T_CD4_Treg",
    TRUE                             ~ "T_cell"
  )
}

map_to_coarse <- function(main_label, fine_label) {
  coarse        <- HPCA_COARSE[main_label]
  coarse[is.na(coarse)] <- "Other"
  is_t          <- !is.na(main_label) & main_label == "T_cells"
  coarse[is_t]  <- tcell_fine_coarse(fine_label[is_t])
  coarse
}

# fine → main lookup built from reference metadata
fine_main_map <- unique(data.frame(
  fine = hpca_ref$label.fine,
  main = hpca_ref$label.main,
  stringsAsFactors = FALSE
))

# ── 4. Colour scheme ──────────────────────────────────────────────────────────

CELL_COLORS <- c(
  T_CD8           = "#e67639",
  T_CD4_Treg      = "#FF9F1C",
  T_cell          = "#b01f22",
  B_cell          = "#0b776c",
  Plasma          = "#4CC9F0",
  NK              = "#06edb0",
  Macrophage_TAM  = "#7B2FBE",
  Monocyte        = "#eeb7f1",
  DC              = "#f725ed",
  Neutrophil_MDSC = "#FFB703",
  Fibroblast_CAF  = "#3A86FF",
  Endothelial     = "#003580",
  Pericyte_SMC    = "#73c5b9",
  Epithelial      = "#699910",
  Neural          = "#786519",
  Other           = "#666666"
)

COMP_COLORS <- c(Tumor = "#FF6B6B", Stroma = "#4CC9F0")
CANCER_GREY    <- "#2A2A2A"
NONCANCER_GREY <- "#333333"

ORIGIN_COLORS <- c(
  Melanoma           = "#FBBF24",
  Bladder_Urothelial = "#BBF7D0",
  Breast             = "#F472B6",
  Lung_Adeno         = "#E63946",
  Thyroid            = "#A78BFA",
  SCLC               = "#C1121F",
  Renal              = "#60A5FA",
  Liver_HCC          = "#F59E0B",
  Rectum             = "#34D399",
  Colorectal         = "#6EE7B7",
  Upper_GI_Gastric   = "#FCA5A5",
  Endometrial_Ovary  = "#DDD6FE",
  Prostate           = "#FEF08A",
  Other              = "#9CA3AF",
  Unknown            = "#9CA3AF"
)

dark_theme <- theme_void(base_size = 9) +
  theme(
    panel.background  = element_rect(fill = "black", colour = NA),
    plot.background   = element_rect(fill = "black", colour = NA),
    legend.background = element_rect(fill = "black", colour = NA),
    legend.text       = element_text(colour = "white", size = 7),
    legend.title      = element_text(colour = "white", face = "bold", size = 8),
    legend.key        = element_rect(fill = "black", colour = NA),
    plot.title        = element_text(colour = "white", face = "bold",
                                     hjust = 0.5, size = 9)
  )

# ── 5. Per-sample loop ────────────────────────────────────────────────────────

all_ann <- list()

for (sname in sample_names) {
  message(sprintf("\n%s\nSample: %s", strrep("─", 50), sname))

  ann_csv <- file.path(PERCELL_DIR, sprintf("%s_cell_annotation.csv", sname))

  # ── Checkpoint ───────────────────────────────────────────────────────────────
  if (file.exists(ann_csv)) {
    message("  [CACHED] Loading annotation CSV")
    all_ann[[sname]] <- read.csv(ann_csv)
    next
  }

  # ── 5a. Load count matrix (BPCells, genes × cells) ───────────────────────────
  bpcells_dir <- file.path(ANALYSIS_ROOT, "per_sample", sname, "bpcells", "counts")
  if (!dir.exists(bpcells_dir)) {
    message(sprintf("  [SKIP] BPCells not found: %s", bpcells_dir)); next
  }

  message("  Loading BPCells matrix...")
  bp_mat  <- open_matrix_dir(bpcells_dir) |> convert_matrix_type("uint32_t")
  meta_s  <- cohort_meta |> filter(sample_name == sname)
  cells_s <- intersect(colnames(bp_mat), meta_s$cell_id)
  message(sprintf("  %s genes × %s cells | %s shared with metadata",
                  format(nrow(bp_mat),    big.mark = ","),
                  format(ncol(bp_mat),    big.mark = ","),
                  format(length(cells_s), big.mark = ",")))

  if (length(cells_s) == 0) { message("  [SKIP] No shared cells"); next }

  # Subset to shared cells — already genes × cells, no transpose needed
  mat_gc <- bp_mat[, cells_s] |> as("dgCMatrix")
  rm(bp_mat); gc()

  # ── 5b. Cancer vs non-cancer (Seurat Compartment_binary) ─────────────────────
  meta_cells <- meta_s |>
    filter(cell_id %in% cells_s) |>
    mutate(is_cancer = Compartment_binary == "Tumor")

  n_cancer <- sum(meta_cells$is_cancer, na.rm = TRUE)
  message(sprintf("  Cancer (Tumor): %s | Non-cancer (Stroma): %s",
                  format(n_cancer,                         big.mark = ","),
                  format(length(cells_s) - n_cancer,       big.mark = ",")))

  # ── 5c. SingleR on ALL cells ──────────────────────────────────────────────────
  message("  Running SingleR (HPCA)...")

  col_sums          <- Matrix::colSums(mat_gc)
  col_sums[col_sums == 0] <- 1
  mat_norm          <- Matrix::t(Matrix::t(mat_gc) / col_sums * 1e4)
  mat_norm@x        <- log1p(mat_norm@x)

  common_genes <- intersect(rownames(mat_norm), rownames(hpca_ref))
  message(sprintf("  Common genes with HPCA: %s / %s",
                  format(length(common_genes), big.mark = ","),
                  format(nrow(hpca_ref),       big.mark = ",")))

  sr_result <- SingleR::SingleR(
    test      = mat_norm[common_genes, ],
    ref       = hpca_ref[common_genes, ],
    labels    = hpca_ref$label.fine,
    BPPARAM   = BiocParallel::MulticoreParam(N_CORES),
    fine.tune = FALSE
  )

  main_labels   <- fine_main_map$main[match(sr_result$labels, fine_main_map$fine)]
  coarse_labels <- map_to_coarse(main_labels, sr_result$labels)

  sr_df <- data.frame(
    cell_id        = cells_s,
    singler_fine   = sr_result$labels,
    singler_main   = main_labels,
    singler_coarse = coarse_labels,
    stringsAsFactors = FALSE
  )

  message("  SingleR coarse label distribution:")
  print(sort(table(sr_df$singler_coarse), decreasing = TRUE))

  # ── 5d. Build final annotation ────────────────────────────────────────────────
  sample_origin <- origin_df |> filter(sample_name == sname) |> pull(cancer_label)
  if (length(sample_origin) == 0) sample_origin <- "Unknown"

  ann_df <- meta_cells |>
    left_join(sr_df |> dplyr::select(cell_id, singler_coarse), by = "cell_id") |>
    mutate(
      annotation = dplyr::if_else(
        is_cancer,
        sample_origin,
        dplyr::coalesce(singler_coarse, "Other")
      )
    )

  write.csv(ann_df, ann_csv, row.names = FALSE)
  message(sprintf("  Saved: %s cells", format(nrow(ann_df), big.mark = ",")))
  all_ann[[sname]] <- ann_df

  rm(mat_gc, mat_norm, sr_result, sr_df); gc()
}

# ── 6. Cohort summaries ───────────────────────────────────────────────────────

if (length(all_ann) == 0) stop("No annotations produced.")

cohort_ann <- dplyr::bind_rows(all_ann)
write.csv(cohort_ann,
          file.path(OUT_DIR, "cohort_cell_annotation.csv"), row.names = FALSE)

cancer_frac <- cohort_ann |>
  group_by(sample_name) |>
  summarise(n_cells         = n(),
            n_cancer        = sum(is_cancer, na.rm = TRUE),
            cancer_fraction = n_cancer / n_cells,
            .groups         = "drop")
write.csv(cancer_frac,
          file.path(OUT_DIR, "sample_cancer_fractions.csv"), row.names = FALSE)

ct_frac <- cohort_ann |>
  group_by(sample_name, annotation) |>
  tally(name = "n") |>
  group_by(sample_name) |>
  mutate(fraction = n / sum(n)) |>
  ungroup()
write.csv(ct_frac,
          file.path(OUT_DIR, "sample_celltype_fractions.csv"), row.names = FALSE)

message(sprintf("\nCohort: %s cells | %d samples",
                format(nrow(cohort_ann), big.mark = ","),
                length(unique(cohort_ann$sample_name))))

# ── 7. Per-sample spatial plots ───────────────────────────────────────────────

message("\nGenerating spatial plots...")

for (sname in sample_names) {
  ann_csv <- file.path(PERCELL_DIR, sprintf("%s_cell_annotation.csv", sname))
  if (!file.exists(ann_csv)) next

  df <- read.csv(ann_csv) |> filter(!is.na(x), !is.na(y))
  if (nrow(df) == 0) next

  pt <- max(0.05, min(1.2, 30000 / nrow(df)))

  c_df  <- df |> filter(is_cancer == TRUE)
  nc_df <- df |> filter(is_cancer == FALSE | is.na(is_cancer))

  sample_origin <- unique(na.omit(c_df$annotation))[1]
  if (is.na(sample_origin) || length(sample_origin) == 0) sample_origin <- "Unknown"
  cancer_col <- dplyr::coalesce(ORIGIN_COLORS[sample_origin], ORIGIN_COLORS["Other"])

  nc_counts <- nc_df |>
    count(annotation, name = "n") |>
    arrange(desc(n)) |>
    mutate(
      col     = dplyr::coalesce(CELL_COLORS[annotation], CELL_COLORS["Other"]),
      leg_lbl = sprintf("%s  (%s)", gsub("_", " ", annotation),
                        format(n, big.mark = ","))
    )
  nc_df <- nc_df |>
    left_join(nc_counts |> dplyr::select(annotation, col, leg_lbl), by = "annotation") |>
    mutate(col     = dplyr::coalesce(col,     CELL_COLORS["Other"]),
           leg_lbl = dplyr::coalesce(leg_lbl, "Other"))

  # Panel 1: Seurat Compartment
  p1 <- ggplot(df, aes(x, y, colour = Compartment_binary)) +
    geom_point(size = pt, stroke = 0, alpha = 0.8) +
    scale_colour_manual(values = COMP_COLORS, name = "Compartment") +
    coord_fixed() + dark_theme +
    ggtitle("Compartment (Seurat)") +
    guides(colour = guide_legend(override.aes = list(size = 3),
                                 title.position = "top"))

  # Panel 2: Cancer cells coloured by tissue origin
  p2 <- ggplot() +
    geom_point(data = nc_df, aes(x, y),
               colour = NONCANCER_GREY, size = pt * 0.4, stroke = 0, alpha = 0.25) +
    { if (nrow(c_df) > 0)
        geom_point(data = c_df, aes(x, y),
                   colour = as.character(cancer_col),
                   size = pt, stroke = 0, alpha = 0.9) } +
    coord_fixed() + dark_theme +
    ggtitle(sprintf("Cancer  —  %s  (%s · %d%%)",
                    sample_origin,
                    format(nrow(c_df), big.mark = ","),
                    round(100 * nrow(c_df) / nrow(df)))) +
    scale_colour_manual(
      values = setNames(as.character(cancer_col), sample_origin),
      name   = "Origin"
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3),
                                 title.position = "top"))

  # Panel 3: Non-cancer cells coloured by SingleR coarse type
  p3 <- ggplot() +
    geom_point(data = c_df, aes(x, y),
               colour = CANCER_GREY, size = pt * 0.4, stroke = 0, alpha = 0.25) +
    geom_point(data = nc_df, aes(x, y, colour = leg_lbl),
               size = pt, stroke = 0, alpha = 0.9) +
    scale_colour_manual(
      values = setNames(nc_df$col, nc_df$leg_lbl),
      breaks = nc_counts$leg_lbl,
      name   = "Cell type"
    ) +
    coord_fixed() + dark_theme +
    ggtitle(sprintf("Non-cancer  (%s · %d%%)",
                    format(nrow(nc_df), big.mark = ","),
                    round(100 * nrow(nc_df) / nrow(df)))) +
    guides(colour = guide_legend(override.aes = list(size = 3), ncol = 1,
                                 title.position = "top"))

  p_all <- (p1 | p2 | p3) +
    plot_annotation(
      title = sprintf("%s  —  Cell Annotation  (Seurat + SingleR/HPCA · %s cells)",
                      sname, format(nrow(df), big.mark = ",")),
      theme = theme(
        plot.background = element_rect(fill = "black", colour = NA),
        plot.title      = element_text(colour = "white", face = "bold",
                                       hjust = 0.5, size = 11)
      )
    )

  out_png <- file.path(PLOT_DIR, sprintf("%s_spatial_annotation.png", sname))
  ggsave(out_png, p_all, width = 22, height = 8, dpi = 150, bg = "black")
  message(sprintf("  %s: saved (%s cells)", sname, format(nrow(df), big.mark = ",")))
}

# ── 8. Cohort summary plots ───────────────────────────────────────────────────

# Non-cancer composition stacked bar
ct_frac_nc <- ct_frac |>
  filter(!annotation %in% names(ORIGIN_COLORS)) |>
  group_by(sample_name) |>
  mutate(frac_nc = n / sum(n)) |>
  ungroup()

sample_order <- ct_frac_nc |>
  filter(annotation == "Macrophage_TAM") |>
  arrange(frac_nc) |>
  pull(sample_name)

nc_bar <- ggplot(ct_frac_nc,
                 aes(factor(sample_name, levels = sample_order),
                     frac_nc, fill = annotation)) +
  geom_bar(stat = "identity", width = 0.8, colour = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = CELL_COLORS,
    labels = function(x) gsub("_", " ", x),
    name   = "Cell type"
  ) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Fraction of non-cancer cells",
       title = "Non-cancer cell composition  (SingleR / HPCA)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(PLOT_DIR, "cohort_noncancer_composition.png"),
       nc_bar, width = 14, height = 6, dpi = 200, bg = "white")

# Cancer fraction bar
cfrac_bar <- ggplot(cancer_frac,
                    aes(reorder(sample_name, cancer_fraction), cancer_fraction)) +
  geom_bar(stat = "identity", fill = "#E63946", colour = "white", linewidth = 0.3) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Cancer cell fraction",
       title = "Cancer cell fraction per sample  (Seurat Compartment_binary)") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(PLOT_DIR, "cohort_cancer_fractions.png"),
       cfrac_bar, width = 12, height = 5, dpi = 200, bg = "white")

message("\n", strrep("=", 60))
message("SCRIPT 11 COMPLETE")
message(sprintf("Outputs : %s", OUT_DIR))
message(sprintf("Plots   : %s", PLOT_DIR))
message(strrep("=", 60))

