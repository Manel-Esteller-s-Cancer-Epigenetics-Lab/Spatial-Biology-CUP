# =============================================================================
# 12c_tumor_origin_IC.R
# -----------------------------------------------------------------------------
# PURPOSE : Tumor-of-origin prediction for IC samples using markers from
#           Table S4 (CK7/CK20 primary + lineage-specific additional markers).
#
#   Marker availability:
#     5K panel (tumour compartment):
#       AMACR, ARG1, CDX2, ESR1, FOLH1, GATA3, KRT20, MME, NKX2-1, NKX3-1,
#       PAX2, PAX8, PGR, SATB2, SCGB1D2, SOX10, SYP, TG, TP63, WT1  (20)
#     MultiTissue panel (all cells — no compartment labels):
#       KRT7, HMGCS2  (+2)
#     Unavailable (absent from both panels):
#       KRT19, NAPSA, SCGB2A2
#
#   Classification strategy:
#     1. Binary positivity call per gene: pct_positive > POS_PCT_THRESHOLD
#     2. Likelihood score per cancer type: (#positive_expected) / (#expected)
#        minus penalty for exclusive negative markers being positive
#     3. Top-scoring type = prediction; second score = alternative
#
# INPUT  : cohort_slim_v2.rds
#          per_sample/<IC>/bpcells/counts  (5K panel)
#          IC_MultiTissue/<STR>/cell_feature_matrix  (MT panel)
#
# OUTPUT : IC_cases/
#   ├── tumor_origin_hscores_IC.csv
#   ├── tumor_origin_classification_IC.csv
#   └── PLOTs/12c_tumor_origin_IC/
#       ├── marker_heatmap_IC.png
#       ├── prediction_bars_IC.png
#       ├── classification_table_IC.png
#       └── spatial/<sample>_spatial_origin.png
# =============================================================================

# ── 0. Configuration ──────────────────────────────────────────────────────────

ANALYSIS_ROOT   <- "."
COHORT_V2_RDS   <- file.path(ANALYSIS_ROOT, "cohort_slim_v2.rds")
OUT_DIR         <- file.path(ANALYSIS_ROOT, "IC_cases")
PLOT_DIR        <- file.path(ANALYSIS_ROOT, "PLOTs", "12c_tumor_origin_IC")
LOG_PATH        <- file.path(ANALYSIS_ROOT, "logs", "12c_tumor_origin_IC.log")
MT_ROOT         <- "path/to/IC_MultiTissue"  # directory containing per-sample Xenium output for IC samples

dir.create(OUT_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ANALYSIS_ROOT, "logs"), recursive = TRUE, showWarnings = FALSE)

IC_SAMPLES <- c("IC_002T1", "IC_004T1", "IC_005T1", "IC_017T2", "IC_029T2")

IC_MT_MAP <- list(
  IC_002T1 = file.path(MT_ROOT, "output-XETG00289__0104567__STR1337__20260514__101648_IC_002-T1"),
  IC_004T1 = file.path(MT_ROOT, "output-XETG00289__0104567__STR1338__20260514__101648_IC_004-T1"),
  IC_005T1 = file.path(MT_ROOT, "output-XETG00289__0104567__STR1339__20260514__101648"),
  IC_017T2 = file.path(MT_ROOT, "output-XETG00289__0104567__STR1340__20260514__101648"),
  IC_029T2 = file.path(MT_ROOT, "output-XETG00289__0104567__STR1341__20260514__101648")
)

# Positivity threshold: ≥5% of tumour cells above raw-count threshold
POS_PCT_THRESHOLD <- 5.0

NEGCTRL_RDS <- file.path(ANALYSIS_ROOT, "Classification", "negcontrol_expression.rds")
if (file.exists(NEGCTRL_RDS)) {
  negctrl_df  <- readRDS(NEGCTRL_RDS)
  negctrl_df$neg_per_probe <- negctrl_df$neg_sum_raw / negctrl_df$n_probes
  nc_p99      <- quantile(negctrl_df$neg_per_probe, 0.99, na.rm = TRUE)
  if (is.na(nc_p99) || nc_p99 == 0)
    nc_p99 <- mean(negctrl_df$neg_per_probe, na.rm = TRUE) +
              2 * sd(negctrl_df$neg_per_probe, na.rm = TRUE)
  H_THRESHOLD_RAW <- nc_p99; rm(negctrl_df)
} else {
  H_THRESHOLD_RAW <- 1
}


# ── Marker sets (gene symbols) ────────────────────────────────────────────────

MARKERS_5K <- c(
  "KRT20",              # CK20 — primary classifier
  "NKX2-1",            # TTF1  — Lung (NSCLC adeno / SCLC), Thyroid
  "SYP",               # Synaptophysin — SCLC
  "TG",                # Thyroglobulin — Thyroid
  "PAX8",              # Thyroid / Endometrial-Ovary / Renal
  "GATA3",             # Breast / Bladder
  "SOX10",             # Breast (triple-negative)
  "ESR1",              # Breast / Endometrial-Ovary
  "PGR",               # Breast / Endometrial-Ovary
  "SCGB1D2",          # Mammaglobin-B — Breast
  "CDX2",             # Upper GI / Colorectal / Gastric
  "WT1",              # Ovary (serous) / Endometrial
  "PAX2",             # Renal / Gynaecological
  "AMACR",            # Renal (papillary) / Prostate
  "MME",              # CD10 — Renal / Lymphoma
  "SATB2",            # Colorectal / Rectum
  "ARG1",             # Hepatocellular
  "CPS1",             # HepPar1 — Hepatocellular (carbamoyl-phosphate synthase 1)
  "FOLH1",            # PSMA — Prostate
  "NKX3-1",          # Prostate
  "TP63"              # Bladder / Squamous
)

MARKERS_MT <- c(
  "KRT7",             # CK7 — primary classifier (from MT panel)
  "HMGCS2"            # HepPar1 equivalent — Hepatocellular
)

ALL_MARKERS <- c(MARKERS_5K, MARKERS_MT)

# ── Cancer type marker profiles (Table S4) ────────────────────────────────────
# positive:  must/should be positive  (scored +1 per marker / total expected)
# negative:  must be negative         (penalty -1 if found positive)
# support:   adds confidence if +ve   (scored +0.5 per marker / total support)
# Note: unavailable markers (KRT19, NAPSA, SCGB2A2) are omitted; CPS1 added from 5K panel

CANCER_PROFILES <- list(
  Lung_Adeno = list(
    label    = "Lung (NSCLC adeno)",
    positive = c("KRT7", "NKX2-1"),
    negative = c("KRT20"),
    support  = c()                   # NAPSA absent
  ),
  Lung_SCLC = list(
    label    = "Lung (SCLC)",
    positive = c("SYP", "NKX2-1"),
    negative = c(),
    support  = c("KRT7")
  ),
  Thyroid = list(
    label    = "Thyroid",
    positive = c("KRT7", "PAX8"),
    negative = c("KRT20"),
    support  = c("TG", "NKX2-1")
  ),
  Breast = list(
    label    = "Breast",
    positive = c("KRT7", "GATA3"),
    negative = c("KRT20"),
    support  = c("SOX10", "ESR1", "PGR", "SCGB1D2")  # SCGB2A2 absent
  ),
  UpperGI_Pancreatobiliary = list(
    label    = "Upper GI / Pancreatobiliary",
    positive = c("KRT7", "CDX2"),
    negative = c(),
    support  = c("KRT20")            # KRT19 absent
  ),
  Endometrial_Ovary = list(
    label    = "Endometrial / Endocervical / Ovary",
    positive = c("KRT7", "PAX8"),
    negative = c("KRT20"),
    support  = c("ESR1", "PGR", "WT1")
  ),
  Renal_Papillary = list(
    label    = "Renal (papillary)",
    positive = c("KRT7", "PAX8", "PAX2"),
    negative = c("KRT20"),
    support  = c("AMACR", "MME")
  ),
  Renal_ClearCell = list(
    label    = "Renal (clear cell)",
    positive = c("PAX8", "PAX2"),
    negative = c("KRT7", "KRT20"),
    support  = c("MME")
  ),
  Bladder = list(
    label    = "Bladder / Urothelial",
    positive = c("KRT7", "GATA3", "TP63"),
    negative = c(),
    support  = c("KRT20")
  ),
  Colorectal = list(
    label    = "Colorectal / Rectum",
    positive = c("KRT20", "CDX2", "SATB2"),
    negative = c("KRT7"),
    support  = c()
  ),
  Hepatocellular = list(
    label    = "Hepatocellular",
    positive = c("ARG1", "CPS1"),
    negative = c("KRT7", "KRT20"),
    support  = c("HMGCS2")
  ),
  Prostate = list(
    label    = "Prostate",
    positive = c("FOLH1", "NKX3-1"),
    negative = c("KRT7", "KRT20"),
    support  = c("AMACR")
  ),
  Gastric = list(
    label    = "Gastric",
    positive = c("CDX2"),
    negative = c("KRT7", "KRT20"),
    support  = c()
  )
)


# ── 1. Libraries ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(BPCells)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(pheatmap)
  library(RColorBrewer)
  library(png)
  library(grid)
  library(gridExtra)
})


# ── 2. Logging ────────────────────────────────────────────────────────────────

.log_con <- file(LOG_PATH, open = "wt")
log_msg <- function(...) {
  txt <- paste0(format(Sys.time(), "[%H:%M:%S] "), paste(..., sep = ""))
  message(txt); tryCatch(writeLines(txt, con = .log_con), error = function(e) NULL)
}
log_cat <- function(...) {
  txt <- paste(..., sep = "")
  cat(txt); tryCatch(writeLines(trimws(txt, "right"), con = .log_con), error = function(e) NULL)
}

log_msg("Script 12c started — Tumor-of-Origin prediction (IC samples)")
log_cat(sprintf("H-Score threshold: %.5f | Positivity call: pct_positive > %.1f%%\n",
                H_THRESHOLD_RAW, POS_PCT_THRESHOLD))


# ── 3. Load metadata ──────────────────────────────────────────────────────────

log_msg("Loading cohort metadata...")
cohort_meta <- readRDS(COHORT_V2_RDS) |>
  dplyr::select(cell_id, sample_name, Compartment_binary, x, y) |>
  filter(sample_name %in% IC_SAMPLES)


# ── 4a. H-Scores: 5K panel (tumour compartment) ───────────────────────────────

log_msg("Extracting 5K-panel H-Scores (tumour compartment)...")

compute_hscores_5k <- function(sname) {
  bpc_dir <- file.path(ANALYSIS_ROOT, "per_sample", sname, "bpcells", "counts")
  if (!dir.exists(bpc_dir)) return(NULL)
  bpc         <- open_matrix_dir(bpc_dir)
  genes_avail <- intersect(MARKERS_5K, rownames(bpc))
  meta_s      <- cohort_meta |> filter(sample_name == sname, Compartment_binary == "Tumor")
  cells_use   <- intersect(meta_s$cell_id, colnames(bpc))
  if (length(cells_use) == 0 || length(genes_avail) == 0) return(NULL)

  expr <- t(as.matrix(bpc[genes_avail, cells_use, drop = FALSE]))

  lapply(genes_avail, function(g) {
    vals <- expr[, g]; n_total <- length(vals)
    n_pos    <- sum(vals > H_THRESHOLD_RAW)
    pct_pos  <- 100 * n_pos / n_total
    mean_int <- if (n_pos > 0) mean(vals[vals > H_THRESHOLD_RAW]) else 0
    data.frame(sample_name = sname, gene = g, panel = "5K",
               n_cells = n_total, n_positive = n_pos, pct_pos = pct_pos,
               mean_int = mean_int, H_Score = pct_pos * mean_int,
               stringsAsFactors = FALSE)
  }) |> dplyr::bind_rows()
}

hscores_5k <- lapply(IC_SAMPLES, function(s) {
  log_msg(sprintf("  5K: %s", s)); compute_hscores_5k(s)
}) |> dplyr::bind_rows()


# ── 4b. H-Scores: MultiTissue panel (all cells, MT-specific threshold) ────────
# Threshold derived from MT NegControl probes per sample (mean + 2×SD).
# For these integer-count data the value falls below 1.0 in all 5 samples,
# making it functionally equivalent to "count ≥ 1" (same as the 5K threshold).
# Computed per-sample so IC_017T2's higher background (thr ≈ 0.016) is handled.

log_msg("Computing per-sample MT NegControl thresholds...")

compute_mt_threshold <- function(sname) {
  mtx_dir  <- file.path(IC_MT_MAP[[sname]], "cell_feature_matrix")
  if (!dir.exists(mtx_dir)) return(H_THRESHOLD_RAW)
  features <- read.table(gzfile(file.path(mtx_dir, "features.tsv.gz")),
                         sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  mat      <- as(Matrix::readMM(gzfile(file.path(mtx_dir, "matrix.mtx.gz"))), "CsparseMatrix")
  neg_idx  <- which(features$V3 == "Negative Control Probe")
  if (length(neg_idx) == 0) return(H_THRESHOLD_RAW)
  per_cell <- as.numeric(Matrix::colSums(mat[neg_idx, , drop = FALSE]) / length(neg_idx))
  mu <- mean(per_cell); sg <- sd(per_cell)
  thr <- mu + 2 * sg
  log_msg(sprintf("  MT threshold %s: mean=%.5f  sd=%.5f  mean+2SD=%.5f  (neg probes=%d)",
                  sname, mu, sg, thr, length(neg_idx)))
  thr
}

MT_THRESHOLDS <- setNames(
  lapply(IC_SAMPLES, compute_mt_threshold),
  IC_SAMPLES
)

log_msg("Extracting MT-panel H-Scores (all cells, no compartment filter)...")

compute_hscores_mt <- function(sname) {
  mtx_dir <- file.path(IC_MT_MAP[[sname]], "cell_feature_matrix")
  if (!dir.exists(mtx_dir)) return(NULL)
  features <- read.table(gzfile(file.path(mtx_dir, "features.tsv.gz")),
                         sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  mat <- as(Matrix::readMM(gzfile(file.path(mtx_dir, "matrix.mtx.gz"))), "CsparseMatrix")
  thr <- MT_THRESHOLDS[[sname]]

  lapply(MARKERS_MT, function(g) {
    gi <- which(features$V2 == g)
    if (length(gi) == 0) return(NULL)
    vals <- as.numeric(mat[gi[1], ]); n_total <- length(vals)
    n_pos    <- sum(vals > thr)
    pct_pos  <- 100 * n_pos / n_total
    mean_int <- if (n_pos > 0) mean(vals[vals > thr]) else 0
    data.frame(sample_name = sname, gene = g, panel = "MT",
               mt_threshold = round(thr, 5),
               n_cells = n_total, n_positive = n_pos, pct_pos = pct_pos,
               mean_int = mean_int, H_Score = pct_pos * mean_int,
               stringsAsFactors = FALSE)
  }) |> dplyr::bind_rows()
}

hscores_mt <- lapply(IC_SAMPLES, function(s) {
  log_msg(sprintf("  MT: %s", s)); compute_hscores_mt(s)
}) |> dplyr::bind_rows()


# ── 4c. Combine and save ──────────────────────────────────────────────────────

marker_scores <- dplyr::bind_rows(hscores_5k, hscores_mt) |>
  mutate(H_Score = round(H_Score, 4), pct_pos = round(pct_pos, 2),
         is_positive = pct_pos > POS_PCT_THRESHOLD,
         mt_threshold = ifelse(is.na(mt_threshold), NA_real_, mt_threshold))

write.csv(marker_scores, file.path(OUT_DIR, "tumor_origin_hscores_IC.csv"), row.names = FALSE)
log_msg("Marker H-Scores saved.")

# Summary of positive calls
pos_summary <- marker_scores |>
  filter(is_positive) |>
  group_by(sample_name) |>
  summarise(positive_markers = paste(sort(gene), collapse = ", "), .groups = "drop")
log_cat("\nPositive markers per sample:\n")
log_cat(capture.output(print(pos_summary, width = 120)), "\n\n")


# ── 5. Classification ─────────────────────────────────────────────────────────

log_msg("Classifying tumor of origin...")

# Build per-sample positivity matrix
pos_matrix <- marker_scores |>
  dplyr::select(sample_name, gene, is_positive) |>
  pivot_wider(names_from = gene, values_from = is_positive, values_fill = FALSE)

score_sample <- function(sname, profile) {
  pos_row <- pos_matrix |> filter(sample_name == sname)
  if (nrow(pos_row) == 0) return(0)

  pos_genes     <- intersect(profile$positive, colnames(pos_row))
  neg_genes     <- intersect(profile$negative, colnames(pos_row))
  support_genes <- intersect(profile$support,  colnames(pos_row))

  n_pos_avail <- length(pos_genes)
  if (n_pos_avail == 0) return(0)

  pos_hits     <- if (n_pos_avail > 0)
    sum(as.logical(pos_row[, pos_genes, drop = TRUE])) / n_pos_avail else 0
  neg_penalty  <- if (length(neg_genes) > 0)
    sum(as.logical(pos_row[, neg_genes, drop = TRUE])) * 0.5 else 0
  sup_hits     <- if (length(support_genes) > 0)
    0.3 * sum(as.logical(pos_row[, support_genes, drop = TRUE])) / length(support_genes) else 0

  score <- pos_hits + sup_hits - neg_penalty
  max(score, 0)
}

classification <- lapply(IC_SAMPLES, function(sname) {
  scores_ct <- sapply(CANCER_PROFILES, function(prof) score_sample(sname, prof))
  names(scores_ct) <- sapply(CANCER_PROFILES, `[[`, "label")
  scores_sorted <- sort(scores_ct, decreasing = TRUE)
  data.frame(
    sample_name   = sname,
    prediction_1  = names(scores_sorted)[1],
    score_1       = round(scores_sorted[1], 3),
    prediction_2  = names(scores_sorted)[2],
    score_2       = round(scores_sorted[2], 3),
    prediction_3  = names(scores_sorted)[3],
    score_3       = round(scores_sorted[3], 3),
    stringsAsFactors = FALSE
  )
}) |> dplyr::bind_rows()

# Add the primary CK7/CK20 pattern
ck_wide <- marker_scores |>
  filter(gene %in% c("KRT7", "KRT20")) |>
  dplyr::select(sample_name, gene, is_positive) |>
  pivot_wider(names_from = gene, values_from = is_positive, values_fill = FALSE)
if (!"KRT7"  %in% colnames(ck_wide)) ck_wide$KRT7  <- FALSE
if (!"KRT20" %in% colnames(ck_wide)) ck_wide$KRT20 <- FALSE

ck_pattern <- ck_wide |>
  mutate(
    CK_pattern = dplyr::case_when(
       KRT7 &  KRT20 ~ "CK7+/CK20+",
       KRT7 & !KRT20 ~ "CK7+/CK20-",
      !KRT7 &  KRT20 ~ "CK7-/CK20+",
      TRUE            ~ "CK7-/CK20-"
    )
  ) |>
  dplyr::select(sample_name, CK_pattern)

classification <- classification |> left_join(ck_pattern, by = "sample_name") |>
  dplyr::select(sample_name, CK_pattern, everything())

write.csv(classification, file.path(OUT_DIR, "tumor_origin_classification_IC.csv"), row.names = FALSE)

log_cat("\nTumor-of-Origin Classification:\n")
log_cat(capture.output(print(classification |> dplyr::select(sample_name, CK_pattern, prediction_1, score_1, prediction_2, score_2))), "\n\n")


# ── 6. Plots ──────────────────────────────────────────────────────────────────

log_msg("Generating plots...")

# Colour palette for cancer types
CANCER_COLS <- c(
  "Lung (NSCLC adeno)"                     = "#E41A1C",
  "Lung (SCLC)"                            = "#FF7F00",
  "Thyroid"                                = "#984EA3",
  "Breast"                                 = "#F781BF",
  "Upper GI / Pancreatobiliary"            = "#A65628",
  "Endometrial / Endocervical / Ovary"     = "#FFFF33",
  "Renal (papillary)"                      = "#4DAF4A",
  "Renal (clear cell)"                     = "#377EB8",
  "Bladder / Urothelial"                   = "#17BECF",
  "Colorectal / Rectum"                    = "#8C564B",
  "Hepatocellular"                         = "#BCBD22",
  "Prostate"                               = "#1F77B4",
  "Gastric"                                = "#9467BD"
)

# ── Plot A: Marker heatmap (z-scored H-Score across IC samples) ───────────────
log_msg("  Plot A: marker heatmap...")

# Z-score each marker across IC samples
z_hm <- marker_scores |>
  group_by(gene) |>
  mutate(z = { mu <- mean(H_Score); sg <- sd(H_Score)
               if (sg > 1e-9) (H_Score - mu) / sg else 0 }) |>
  ungroup()

gene_order <- c(
  "KRT7", "KRT20",
  "NKX2-1", "SYP",                          # Lung
  "TG", "PAX8",                              # Thyroid (+shared)
  "GATA3", "SOX10", "ESR1", "PGR", "SCGB1D2", # Breast (+shared)
  "CDX2",                                    # GI/Colorectal
  "WT1",                                     # Ovary/Endometrial
  "PAX2", "AMACR", "MME",                   # Renal
  "SATB2",                                   # Colorectal
  "ARG1", "HMGCS2",                          # Hepatocellular
  "FOLH1", "NKX3-1",                        # Prostate
  "TP63"                                     # Bladder
)
gene_order <- intersect(gene_order, unique(z_hm$gene))

z_wide <- z_hm |>
  filter(gene %in% gene_order) |>
  dplyr::select(sample_name, gene, z) |>
  pivot_wider(names_from = gene, values_from = z, values_fill = 0)

mat_hm <- as.matrix(z_wide[, gene_order[gene_order %in% colnames(z_wide)], drop = FALSE])
rownames(mat_hm) <- z_wide$sample_name

# Positivity overlay (asterisk in cell text)
pos_wide <- marker_scores |>
  filter(gene %in% colnames(mat_hm)) |>
  dplyr::select(sample_name, gene, is_positive) |>
  pivot_wider(names_from = gene, values_from = is_positive, values_fill = FALSE) |>
  arrange(match(sample_name, rownames(mat_hm)))
pos_mat <- as.matrix(pos_wide[, colnames(mat_hm), drop = FALSE])
rownames(pos_mat) <- pos_wide$sample_name
display_mat <- matrix(
  ifelse(pos_mat[rownames(mat_hm), colnames(mat_hm)], "*", ""),
  nrow = nrow(mat_hm), dimnames = dimnames(mat_hm)
)

# Annotation: panel source and marker category
ann_col <- data.frame(
  Panel = ifelse(colnames(mat_hm) %in% MARKERS_MT, "MultiTissue", "5K"),
  Category = dplyr::case_when(
    colnames(mat_hm) %in% c("KRT7","KRT20")           ~ "Primary",
    colnames(mat_hm) %in% c("NKX2-1","SYP")           ~ "Lung / SCLC",
    colnames(mat_hm) %in% c("TG")                     ~ "Thyroid",
    colnames(mat_hm) %in% c("PAX8","PAX2")            ~ "Thyroid / Renal / Gyn",
    colnames(mat_hm) %in% c("GATA3","SOX10","ESR1","PGR","SCGB1D2") ~ "Breast",
    colnames(mat_hm) %in% c("WT1")                    ~ "Gyn / Ovary",
    colnames(mat_hm) %in% c("CDX2","SATB2")           ~ "GI / Colorectal",
    colnames(mat_hm) %in% c("AMACR","MME")            ~ "Renal",
    colnames(mat_hm) %in% c("ARG1","HMGCS2")          ~ "Hepatocellular",
    colnames(mat_hm) %in% c("FOLH1","NKX3-1")        ~ "Prostate",
    colnames(mat_hm) %in% c("TP63")                   ~ "Bladder",
    TRUE ~ "Other"
  ),
  row.names = colnames(mat_hm), stringsAsFactors = FALSE
)

ann_row <- data.frame(
  CK_pattern  = classification$CK_pattern[match(rownames(mat_hm), classification$sample_name)],
  Prediction  = classification$prediction_1[match(rownames(mat_hm), classification$sample_name)],
  row.names   = rownames(mat_hm), stringsAsFactors = FALSE
)

CAT_COLS <- c("Primary" = "#2C2C2C", "Lung / SCLC" = "#E41A1C",
              "Thyroid" = "#984EA3", "Thyroid / Renal / Gyn" = "#A29BD4",
              "Breast" = "#F781BF", "Gyn / Ovary" = "#FFDB58",
              "GI / Colorectal" = "#8C564B", "Renal" = "#377EB8",
              "Hepatocellular" = "#BCBD22", "Prostate" = "#1F77B4",
              "Bladder" = "#17BECF", "Other" = "#AAAAAA")
PANEL_COLS <- c("5K" = "#4CAF50", "MultiTissue" = "#FF9800")
CK_COLS    <- c("CK7+/CK20-" = "#D73027","CK7+/CK20+" = "#FC8D59",
                "CK7-/CK20+" = "#4575B4","CK7-/CK20-" = "#A6A6A6")
PRED_LABELS <- unique(ann_row$Prediction)
PRED_COLS   <- CANCER_COLS[PRED_LABELS]
PRED_COLS[is.na(PRED_COLS)] <- "#AAAAAA"

png(file.path(PLOT_DIR, "marker_heatmap_IC.png"),
    width = 16, height = 6, units = "in", res = 200, bg = "white")
pheatmap(
  mat_hm, cluster_rows = FALSE, cluster_cols = FALSE,
  display_numbers = display_mat, number_color = "black", fontsize_number = 11,
  annotation_col = ann_col, annotation_row = ann_row,
  annotation_colors = list(
    Panel      = PANEL_COLS,
    Category   = CAT_COLS,
    CK_pattern = CK_COLS,
    Prediction = PRED_COLS
  ),
  color      = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks     = seq(-3, 3, length.out = 101),
  border_color = "grey80", fontsize = 9, fontsize_col = 8, fontsize_row = 9,
  main   = "Tumor-of-Origin Markers — IC Cohort (z-scored H-Score)\n* = positive call (pct_positive > 5%)",
  silent = TRUE
)
dev.off()
log_msg("  Plot A saved.")


# ── Plot B: Prediction confidence bars ────────────────────────────────────────
log_msg("  Plot B: prediction bars...")

# Build long table of all scores for top-3 per sample
all_scores <- lapply(IC_SAMPLES, function(sname) {
  scores_ct <- sapply(CANCER_PROFILES, function(prof) score_sample(sname, prof))
  names(scores_ct) <- sapply(CANCER_PROFILES, `[[`, "label")
  data.frame(sample_name = sname, cancer_type = names(scores_ct),
             score = round(scores_ct, 3), stringsAsFactors = FALSE)
}) |> dplyr::bind_rows()

# Keep all types that have score > 0 in at least one sample
types_with_signal <- all_scores |>
  group_by(cancer_type) |>
  summarise(max_score = max(score)) |>
  filter(max_score > 0) |>
  pull(cancer_type)

plot_scores <- all_scores |>
  filter(cancer_type %in% types_with_signal) |>
  mutate(
    sample_name = factor(sample_name, levels = IC_SAMPLES),
    cancer_type = factor(cancer_type, levels = rev(unique(cancer_type)))
  )

p_bars <- ggplot(plot_scores, aes(x = score, y = cancer_type, fill = cancer_type)) +
  geom_col(show.legend = FALSE) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = CANCER_COLS, na.value = "#AAAAAA") +
  scale_x_continuous(limits = c(0, 1.4), breaks = c(0, 0.5, 1.0)) +
  facet_wrap(~ sample_name, nrow = 1) +
  labs(title    = "Tumor-of-Origin Prediction Scores — IC Cohort",
       subtitle = "Score = positive markers / total expected − 0.5×penalised negatives + 0.3×supporting\nDashed line = 0.5 confidence threshold",
       x = "Score", y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    plot.title     = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle  = element_text(hjust = 0.5, colour = "grey40", size = 8),
    strip.text     = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95"),
    axis.text.y    = element_text(size = 8)
  )

ggsave(file.path(PLOT_DIR, "prediction_bars_IC.png"),
       p_bars, width = 14, height = 6, dpi = 200, bg = "white")
log_msg("  Plot B saved.")


# ── Plot C: Classification summary table ──────────────────────────────────────
log_msg("  Plot C: classification table...")

# Positive markers per sample
pos_per_sample <- marker_scores |>
  filter(is_positive) |>
  group_by(sample_name) |>
  summarise(positive = paste(sort(gene), collapse = "\n"), .groups = "drop")

table_df <- classification |>
  left_join(pos_per_sample, by = "sample_name") |>
  dplyr::select(sample_name, CK_pattern, prediction_1, score_1, prediction_2, score_2, positive) |>
  mutate(score_1 = round(score_1, 2), score_2 = round(score_2, 2))

colnames(table_df) <- c("Sample", "CK7/CK20", "Prediction 1", "Score 1",
                         "Prediction 2", "Score 2", "Positive markers")

png(file.path(PLOT_DIR, "classification_table_IC.png"),
    width = 16, height = 4, units = "in", res = 200, bg = "white")
grid.newpage()
gt <- gridExtra::tableGrob(table_df, rows = NULL,
        theme = gridExtra::ttheme_default(
          base_size = 9,
          core    = list(fg_params = list(fontsize = 8)),
          colhead = list(fg_params = list(fontsize = 9, fontface = "bold"))
        ))
grid.draw(gt)
dev.off()
log_msg("  Plot C saved.")


# ── Plot D: Spatial expression of key markers ─────────────────────────────────
log_msg(strrep("─", 60))
log_msg("  Plot D: spatial marker plots...")

SPATIAL_DIR  <- file.path(PLOT_DIR, "spatial")
HE_THUMB_DIR <- file.path(PLOT_DIR, "he_thumbs")
dir.create(SPATIAL_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(HE_THUMB_DIR, recursive = TRUE, showWarnings = FALSE)

HE_DIR <- file.path(MT_ROOT, "HE_POST_XENIUM_20X")  # folder of 20x H&E TIF images
HE_SOURCE_MAP <- c(
  IC_002T1 = file.path(ANALYSIS_ROOT, "New_Plots", "IC-002T1_HE_thumb.png"),
  IC_004T1 = file.path(ANALYSIS_ROOT, "New_Plots", "IC-004T1_HE_thumb.png"),
  IC_005T1 = file.path(HE_DIR, "STR1339.tif"),
  IC_017T2 = file.path(HE_DIR, "STR1340.tif"),
  IC_029T2 = file.path(HE_DIR, "STR1341.tif")
)

get_he_thumb <- function(sname, thumb_dir) {
  src <- HE_SOURCE_MAP[sname]
  if (is.na(src)) return(NULL)
  if (grepl("\\.png$", src, ignore.case = TRUE)) {
    return(if (file.exists(src)) src else NULL)
  }
  thumb_png <- file.path(thumb_dir, paste0(sname, "_he_thumb.png"))
  if (file.exists(thumb_png)) return(thumb_png)
  if (!file.exists(src)) return(NULL)
  log_msg(sprintf("    Generating H&E thumbnail: %s ...", sname))
  cmd <- sprintf('convert "%s[0]" -thumbnail 1200x -normalize -colorspace sRGB "%s"',
                 src, thumb_png)
  ret <- system(cmd, intern = FALSE, ignore.stderr = TRUE)
  if (ret == 0 && file.exists(thumb_png)) thumb_png else NULL
}

dark_theme <- theme_void(base_size = 10) +
  theme(
    panel.background  = element_rect(fill = "black", colour = NA),
    plot.background   = element_rect(fill = "black", colour = NA),
    legend.background = element_rect(fill = "black", colour = NA),
    legend.key        = element_rect(fill = "black", colour = NA),
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 9, colour = "white"),
    legend.position   = "bottom",
    legend.text       = element_text(colour = "white", size = 7),
    legend.title      = element_text(colour = "white", size = 7)
  )

# Fixed set of spatial marker panels per sample:
# H&E | CK7 (MT) | CK20 | Top-3 positive markers (from 5K)
# For CK7 we draw from MT coordinates — show as a separate 5K-only panel
SPATIAL_GENES_5K <- c("KRT20", "NKX2-1", "GATA3", "SOX10", "ESR1", "PGR",
                       "CDX2", "PAX8", "PAX2", "SATB2", "FOLH1", "NKX3-1",
                       "ARG1", "TP63", "WT1", "SYP", "TG", "AMACR", "MME",
                       "SCGB1D2")

spatial_meta <- cohort_meta |> dplyr::select(cell_id, sample_name, x, y, Compartment_binary)

build_marker_panels <- function(sname) {
  bpc_dir <- file.path(ANALYSIS_ROOT, "per_sample", sname, "bpcells", "counts")
  if (!dir.exists(bpc_dir)) return(NULL)

  bpc         <- open_matrix_dir(bpc_dir)
  genes_avail <- intersect(SPATIAL_GENES_5K, rownames(bpc))
  meta_s      <- spatial_meta |> filter(sample_name == sname)
  cells_use   <- intersect(meta_s$cell_id, colnames(bpc))
  if (length(cells_use) == 0) return(NULL)

  expr_mat <- as.matrix(bpc[genes_avail, cells_use, drop = FALSE])
  plot_df  <- meta_s |> filter(cell_id %in% cells_use) |>
    left_join(
      as.data.frame(t(log1p(expr_mat))) |> tibble::rownames_to_column("cell_id"),
      by = "cell_id"
    ) |> filter(!is.na(x), !is.na(y))
  if (nrow(plot_df) == 0) return(NULL)

  pt <- if (nrow(plot_df) > 100000) 0.08 else if (nrow(plot_df) > 50000) 0.15 else 0.3

  # Identify which markers are positive for this sample (prioritise for display)
  pos_genes_sample <- marker_scores |>
    filter(sample_name == sname, is_positive, gene %in% genes_avail) |>
    arrange(desc(pct_pos)) |>
    pull(gene)

  # Always show: KRT20; then top-positive up to 5 additional slots
  show_genes <- unique(c("KRT20", pos_genes_sample,
                          genes_avail[!genes_avail %in% c("KRT20", pos_genes_sample)]))[1:6]
  show_genes <- show_genes[show_genes %in% colnames(plot_df)]

  # H&E panel
  thumb <- get_he_thumb(sname, HE_THUMB_DIR)
  if (!is.null(thumb) && file.exists(thumb)) {
    img <- png::readPNG(thumb)
    if (length(dim(img)) == 3) img <- img[nrow(img):1, , , drop = FALSE]
    else                        img <- img[nrow(img):1,  ,  drop = FALSE]
    if (length(dim(img)) == 3) {
      wmask <- img[,,1] > 0.92 & img[,,2] > 0.92 & img[,,3] > 0.92
      img[,,1][wmask] <- 0; img[,,2][wmask] <- 0; img[,,3][wmask] <- 0
    }
    p_he <- ggplot() +
      annotation_raster(as.raster(img), xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
      xlim(0,1) + ylim(0,1) + coord_fixed(nrow(img)/ncol(img)) +
      labs(title = "H&E") + dark_theme
  } else {
    p_he <- ggplot() +
      annotate("text", x=0.5, y=0.5, label="H&E\n(unavailable)", colour="grey60", size=3.5) +
      xlim(0,1)+ylim(0,1)+coord_fixed()+labs(title="H&E")+dark_theme
  }

  # Marker spatial panels
  marker_panels <- lapply(show_genes, function(g) {
    q95 <- max(quantile(plot_df[[g]], 0.95, na.rm = TRUE), 0.01)
    is_pos_sample <- marker_scores |>
      filter(sample_name == sname, gene == g) |> pull(is_positive)
    is_pos_sample <- if (length(is_pos_sample)) is_pos_sample[1] else FALSE
    title_txt <- if (is_pos_sample) paste0(g, " ✓") else g

    ggplot(plot_df |> arrange(.data[[g]]),
           aes(x = x, y = y, colour = pmin(.data[[g]], q95))) +
      geom_point(size = pt, stroke = 0, alpha = 0.85) +
      scale_colour_viridis_c(option = "inferno", begin = 0.05, name = g,
        limits = c(0, q95),
        guide = guide_colourbar(barwidth = 4, barheight = 0.5,
                                title.position = "top", title.hjust = 0.5)) +
      coord_equal() + labs(title = title_txt) + dark_theme
  })

  c(list(he = p_he), setNames(marker_panels, show_genes))
}

for (sname in IC_SAMPLES) {
  log_msg(sprintf("    %s ...", sname))
  panels <- build_marker_panels(sname)
  if (is.null(panels)) next

  all_panels <- panels
  p1_labelled <- all_panels[[1]] +
    labs(y = gsub("_", "-", sname)) +
    theme(axis.title.y = element_text(colour = "white", size = 11, face = "bold",
                                      angle = 90, vjust = 0.5, margin = margin(r = 6)))
  all_panels[[1]] <- p1_labelled

  fig <- patchwork::wrap_plots(all_panels, nrow = 1) +
    patchwork::plot_annotation(
      theme = theme(plot.background = element_rect(fill = "black", colour = NA))
    ) & theme(plot.background = element_rect(fill = "black", colour = NA))

  ggsave(file.path(SPATIAL_DIR, paste0(sname, "_spatial_origin.png")),
         fig, width = 7 * length(all_panels), height = 6, dpi = 200, bg = "black")
  log_msg(sprintf("    Saved: %s", paste0(sname, "_spatial_origin.png")))
}

# Consolidated spatial: H&E + KRT20 + top-positive marker per sample
log_msg("  Building consolidated spatial overview...")
ic_rows_consol <- lapply(IC_SAMPLES, function(sname) {
  panels <- build_marker_panels(sname)
  if (is.null(panels)) return(NULL)
  # take H&E + first 4 marker panels
  sel <- panels[1:min(5, length(panels))]
  p1_lab <- sel[[1]] +
    labs(y = gsub("_", "-", sname)) +
    theme(axis.title.y = element_text(colour = "white", size = 11, face = "bold",
                                      angle = 90, vjust = 0.5, margin = margin(r = 6)))
  sel[[1]] <- p1_lab
  patchwork::wrap_plots(sel, nrow = 1) +
    patchwork::plot_annotation(
      theme = theme(plot.background = element_rect(fill = "black", colour = NA))
    ) & theme(plot.background = element_rect(fill = "black", colour = NA))
})
ic_rows_consol <- Filter(Negate(is.null), ic_rows_consol)

if (length(ic_rows_consol) > 0) {
  ic_consol <- patchwork::wrap_plots(ic_rows_consol, ncol = 1) +
    patchwork::plot_annotation(
      title   = "IC Cohort — Tumor-of-Origin Marker Expression",
      caption = "H&E | CK20 | positive markers shown first | * = positive call",
      theme   = theme(
        plot.title   = element_text(face = "bold", hjust = 0.5, size = 16, colour = "white"),
        plot.caption = element_text(hjust = 0.5, size = 9, colour = "grey60"),
        plot.background = element_rect(fill = "black", colour = NA)
      )
    ) & theme(plot.background = element_rect(fill = "black", colour = NA))

  ggsave(file.path(PLOT_DIR, "IC_combined_origin_spatial.png"), ic_consol,
         width = 35, height = length(ic_rows_consol) * 6, dpi = 200, bg = "black")
  log_msg("  Consolidated spatial saved.")
}

log_msg("  Plot D done.")


# ── Final summary ─────────────────────────────────────────────────────────────

log_msg(strrep("=", 60))
log_msg("SCRIPT 12c COMPLETE — Tumor-of-Origin Prediction")
log_msg(sprintf("Markers used: %d 5K + %d MT = %d total",
                length(MARKERS_5K), length(MARKERS_MT),
                length(MARKERS_5K) + length(MARKERS_MT)))
log_msg("Unavailable markers: KRT19, NAPSA, SCGB2A2 (absent from both panels)")
log_cat("\nFinal predictions:\n")
log_cat(capture.output(print(
  classification |> dplyr::select(sample_name, CK_pattern, prediction_1, score_1)
)), "\n")
log_msg(strrep("=", 60))

close(.log_con)
