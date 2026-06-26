# =============================================================================
# 12b_emt_stemness_IC_only.R
# -----------------------------------------------------------------------------
# PURPOSE : IC-specific EMT & Stemness analysis incorporating genes from the
#           Xenium MultiTissue panel that were absent from the 5K panel.
#
#   Additional genes recovered from MultiTissue panel (541 genes):
#     ALDH1A3  → EMT_T1_POS (T1 progenitor/stemness)
#     LY6D     → EMT_T2_STROMA (T2 inflammatory)
#
#   All other absent genes (VIM, SPARC, S100A4, EGR1/JUN/FOS, etc.) are not
#   present in either panel.
#
#   Key differences vs script 12:
#     • IC samples only (n=5); z-scoring across these 5 samples
#     • MT-panel H-Scores use ALL MT cells (no tumour/stroma filter — MT cells
#       are from a separate Xenium run and lack CancerFinder3 labels)
#     • 5K-panel H-Scores use tumour compartment as before
#
# INPUT  : cohort_slim_v2.rds
#          per_sample/<IC>/bpcells/counts  (5K panel)
#          IC_MultiTissue/<STR>/cell_feature_matrix  (MT panel)
#
# OUTPUT : EMT_Stemness/
#   ├── emt_gene_hscores_IC.csv
#   ├── emt_stemness_classification_IC.csv
#   └── PLOTs/12b_emt_stemness_IC/
#       ├── emt_gene_heatmap_IC.png
#       ├── emt_scatter_IC.png
#       ├── emt_score_bars_IC.png
#       ├── stemness_score_bars_IC.png
#       └── spatial/<sample>_spatial_emt_stemness_IC.png
# =============================================================================

# ── 0. Configuration ──────────────────────────────────────────────────────────

ANALYSIS_ROOT <- "."
COHORT_V2_RDS <- file.path(ANALYSIS_ROOT, "cohort_slim_v2.rds")
OUT_DIR       <- file.path(ANALYSIS_ROOT, "EMT_Stemness")
PLOT_DIR      <- file.path(ANALYSIS_ROOT, "PLOTs", "12b_emt_stemness_IC")
LOG_PATH      <- file.path(ANALYSIS_ROOT, "logs", "12b_emt_stemness_IC.log")
MT_ROOT       <- "path/to/IC_MultiTissue"  # directory containing per-sample Xenium output for IC samples

dir.create(OUT_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ANALYSIS_ROOT, "logs"), recursive = TRUE, showWarnings = FALSE)

# IC sample → MultiTissue directory
IC_MT_MAP <- list(
  IC_002T1 = file.path(MT_ROOT, "output-XETG00289__0104567__STR1337__20260514__101648_IC_002-T1"),
  IC_004T1 = file.path(MT_ROOT, "output-XETG00289__0104567__STR1338__20260514__101648_IC_004-T1"),
  IC_005T1 = file.path(MT_ROOT, "output-XETG00289__0104567__STR1339__20260514__101648"),
  IC_017T2 = file.path(MT_ROOT, "output-XETG00289__0104567__STR1340__20260514__101648"),
  IC_029T2 = file.path(MT_ROOT, "output-XETG00289__0104567__STR1341__20260514__101648")
)

IC_SAMPLES <- names(IC_MT_MAP)

# ── Positivity threshold (inherited from Script 03 / Script 12) ───────────────
NEGCTRL_RDS <- file.path(ANALYSIS_ROOT, "Classification", "negcontrol_expression.rds")
if (file.exists(NEGCTRL_RDS)) {
  negctrl_df          <- readRDS(NEGCTRL_RDS)
  negctrl_df$neg_per_probe <- negctrl_df$neg_sum_raw / negctrl_df$n_probes
  negctrl_p99         <- quantile(negctrl_df$neg_per_probe, 0.99, na.rm = TRUE)
  if (is.na(negctrl_p99) || negctrl_p99 == 0)
    negctrl_p99 <- mean(negctrl_df$neg_per_probe, na.rm = TRUE) +
                   2 * sd(negctrl_df$neg_per_probe, na.rm = TRUE)
  H_THRESHOLD_RAW <- negctrl_p99
  rm(negctrl_df)
} else {
  H_THRESHOLD_RAW <- 1
}

# ── Gene sets ─────────────────────────────────────────────────────────────────

# ALDH1A3 added (T1 progenitor) — from MT panel
EMT_T1_POS <- c(
  "SNAI1", "SNAI2", "TWIST1", "ZEB1", "PRRX1",
  "CDH2", "TNC", "PALLD", "MYLK", "PDPN", "NRP2",
  "MMP13", "MMP14", "LOXL1",
  "NDRG1", "IGF1", "CDC20",
  "ALDH1A3"
)
EMT_T1_EPI    <- c("EPCAM", "CLDN7")
EMT_T2_TUMOR  <- c("SNAI1", "KLF4", "MAFB", "POSTN", "TIMP2", "CCND1", "RB1", "NF1")
# LY6D added (T2 inflammatory stroma) — from MT panel
EMT_T2_STROMA <- c("CCRL2", "NOTCH2", "IRF7", "CXCL16", "LY6D")
STEMNESS_GENES <- c("KLF5", "MYBL2", "NFE2L3", "ILF3", "HMGA1", "HMGB3")

MT_GENES        <- c("ALDH1A3", "LY6D")
ALL_EMT_GENES   <- unique(c(EMT_T1_POS, EMT_T1_EPI, EMT_T2_TUMOR, EMT_T2_STROMA, STEMNESS_GENES))
ALL_5K_GENES    <- setdiff(ALL_EMT_GENES, MT_GENES)


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
  library(ggrepel)
  library(png)
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

log_msg("Script 12b started — IC-only EMT & Stemness (+ MultiTissue panel genes)")
log_cat(sprintf("H-Score threshold (raw): %.5f\n", H_THRESHOLD_RAW))
log_cat(sprintf("5K-panel genes: %d | MT-panel additions: %s\n",
                length(ALL_5K_GENES), paste(MT_GENES, collapse = ", ")))


# ── 3. Load cohort metadata ───────────────────────────────────────────────────

log_msg("Loading cohort metadata...")
cohort_meta <- readRDS(COHORT_V2_RDS) |>
  dplyr::select(cell_id, sample_name, Compartment_binary, x, y) |>
  filter(sample_name %in% IC_SAMPLES)

log_cat(sprintf("IC cells: %d across %d samples\n",
                nrow(cohort_meta), length(IC_SAMPLES)))


# ── 4a. H-Scores from 5K BPCells (tumour compartment) ────────────────────────

log_msg(strrep("─", 60))
log_msg("Extracting 5K-panel H-Scores (tumour compartment)...")

compute_hscores_5k <- function(sname) {
  bpc_dir <- file.path(ANALYSIS_ROOT, "per_sample", sname, "bpcells", "counts")
  if (!dir.exists(bpc_dir)) { log_msg(sprintf("  [SKIP] %s", sname)); return(NULL) }

  bpc         <- open_matrix_dir(bpc_dir)
  genes_avail <- intersect(ALL_5K_GENES, rownames(bpc))
  meta_s      <- cohort_meta |> filter(sample_name == sname, Compartment_binary == "Tumor")
  cells_use   <- intersect(meta_s$cell_id, colnames(bpc))

  if (length(cells_use) == 0 || length(genes_avail) == 0) return(NULL)

  expr <- t(as.matrix(bpc[genes_avail, cells_use, drop = FALSE]))

  lapply(genes_avail, function(g) {
    vals       <- expr[, g]
    n_total    <- length(vals)
    n_positive <- sum(vals > H_THRESHOLD_RAW)
    pct_pos    <- 100 * n_positive / n_total
    mean_int   <- if (n_positive > 0) mean(vals[vals > H_THRESHOLD_RAW]) else 0
    data.frame(sample_name = sname, gene = g, panel = "5K",
               compartment = "Tumor", n_cells = n_total,
               n_positive = n_positive, pct_pos = pct_pos,
               mean_int = mean_int, H_Score = pct_pos * mean_int,
               stringsAsFactors = FALSE)
  }) |> dplyr::bind_rows()
}

hscores_5k <- lapply(IC_SAMPLES, function(s) {
  log_msg(sprintf("  %s ...", s)); compute_hscores_5k(s)
}) |> dplyr::bind_rows()


# ── 4b. H-Scores from MultiTissue MTX (all cells — no compartment labels) ────

log_msg(strrep("─", 60))
log_msg("Extracting MultiTissue-panel H-Scores (all cells, no compartment filter)...")

compute_hscores_mt <- function(sname) {
  mt_dir  <- IC_MT_MAP[[sname]]
  mtx_dir <- file.path(mt_dir, "cell_feature_matrix")
  if (!dir.exists(mtx_dir)) { log_msg(sprintf("  [SKIP] MT dir missing: %s", sname)); return(NULL) }

  barcodes <- readLines(gzfile(file.path(mtx_dir, "barcodes.tsv.gz")))
  features <- read.table(gzfile(file.path(mtx_dir, "features.tsv.gz")),
                         sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  mat      <- as(Matrix::readMM(gzfile(file.path(mtx_dir, "matrix.mtx.gz"))), "dgCMatrix")

  lapply(MT_GENES, function(g) {
    gi <- which(features$V2 == g)
    if (length(gi) == 0) return(NULL)
    vals       <- as.numeric(mat[gi[1], ])
    n_total    <- length(vals)
    n_positive <- sum(vals > H_THRESHOLD_RAW)
    pct_pos    <- 100 * n_positive / n_total
    mean_int   <- if (n_positive > 0) mean(vals[vals > H_THRESHOLD_RAW]) else 0
    data.frame(sample_name = sname, gene = g, panel = "MultiTissue",
               compartment = "All", n_cells = n_total,
               n_positive = n_positive, pct_pos = pct_pos,
               mean_int = mean_int, H_Score = pct_pos * mean_int,
               stringsAsFactors = FALSE)
  }) |> dplyr::bind_rows()
}

hscores_mt <- lapply(IC_SAMPLES, function(s) {
  log_msg(sprintf("  %s ...", s)); compute_hscores_mt(s)
}) |> dplyr::bind_rows()


# ── 4c. Combine and save H-Scores ─────────────────────────────────────────────

emt_hscores <- dplyr::bind_rows(hscores_5k, hscores_mt) |>
  mutate(H_Score = round(H_Score, 4))

write.csv(emt_hscores, file.path(OUT_DIR, "emt_gene_hscores_IC.csv"), row.names = FALSE)
log_msg("IC gene H-Scores saved.")
log_cat(sprintf("\nGenes scored (5K): %s\n",
                paste(sort(unique(hscores_5k$gene)), collapse = ", ")))
log_cat(sprintf("Genes scored (MT): %s\n\n",
                paste(sort(unique(hscores_mt$gene)), collapse = ", ")))


# ── 5. Z-scored EMT and Stemness scores (across IC samples only) ──────────────

log_msg(strrep("─", 60))
log_msg("Computing EMT-T1, EMT-T2, and Stemness scores (z-scored across IC n=5)...")

# For MT genes the compartment label is "All"; treat them alongside tumour scores
pivot_wide <- function(panel_label, comp_label) {
  emt_hscores |>
    filter(panel == panel_label, compartment == comp_label) |>
    dplyr::select(sample_name, gene, H_Score) |>
    pivot_wider(names_from = gene, values_from = H_Score, values_fill = 0)
}

# 5K tumour-compartment wide matrix + MT genes wide matrix (all cells)
wide_tumor <- pivot_wide("5K", "Tumor")
wide_mt    <- pivot_wide("MultiTissue", "All")

# Join MT columns into tumour-wide frame (same sample rows, new gene columns)
wide_combined <- wide_tumor |>
  left_join(wide_mt |> dplyr::select(sample_name, any_of(MT_GENES)),
            by = "sample_name")

zscore_cols <- function(df) {
  gene_cols <- setdiff(colnames(df), "sample_name")
  for (g in gene_cols) {
    mu <- mean(df[[g]], na.rm = TRUE); sg <- sd(df[[g]], na.rm = TRUE)
    df[[g]] <- if (sg > 1e-9) (df[[g]] - mu) / sg else 0
  }
  df
}

z_combined <- zscore_cols(wide_combined)

mean_z <- function(zdf, genes) {
  avail <- intersect(genes, colnames(zdf))
  if (length(avail) == 0) return(rep(0, nrow(zdf)))
  rowMeans(as.matrix(zdf[, avail, drop = FALSE]), na.rm = TRUE)
}

scores <- data.frame(sample_name = z_combined$sample_name)

scores$T1_positive_score   <- mean_z(z_combined, intersect(EMT_T1_POS, colnames(z_combined)))
scores$T1_epithelial_score <- mean_z(z_combined, intersect(EMT_T1_EPI, colnames(z_combined)))
scores$EMT_T1_score        <- scores$T1_positive_score - scores$T1_epithelial_score

scores$EMT_T2_tumor_score  <- mean_z(z_combined, intersect(EMT_T2_TUMOR,  colnames(z_combined)))
scores$EMT_T2_stroma_score <- mean_z(z_combined, intersect(EMT_T2_STROMA, colnames(z_combined)))
scores$EMT_T2_score        <- (scores$EMT_T2_tumor_score + scores$EMT_T2_stroma_score) / 2

scores$Stemness_score      <- mean_z(z_combined, intersect(STEMNESS_GENES, colnames(z_combined)))
scores$EMT_delta           <- scores$EMT_T1_score - scores$EMT_T2_score

MIXED_THRESHOLD <- 0.5
scores$EMT_type <- dplyr::case_when(
  scores$EMT_T1_score <= 0 & scores$EMT_T2_score <= 0                   ~ "No_EMT",
  abs(scores$EMT_delta) <= MIXED_THRESHOLD &
    scores$EMT_T1_score > 0 & scores$EMT_T2_score > 0                   ~ "Mixed",
  scores$EMT_T1_score > scores$EMT_T2_score & scores$EMT_T1_score > 0   ~ "EMT_T1",
  scores$EMT_T2_score > scores$EMT_T1_score & scores$EMT_T2_score > 0   ~ "EMT_T2",
  TRUE                                                                   ~ "No_EMT"
)
scores$Stemness_class <- ifelse(scores$Stemness_score > 0, "Stemness_High", "Stemness_Low")
scores$Combined_class <- paste0(scores$EMT_type, " / ", scores$Stemness_class)
scores <- scores |> mutate(across(where(is.numeric), ~ round(.x, 4)))

write.csv(scores, file.path(OUT_DIR, "emt_stemness_classification_IC.csv"), row.names = FALSE)

log_cat("\nIC EMT & Stemness classification (with MT panel genes):\n")
log_cat(capture.output(print(
  scores |> dplyr::select(sample_name, EMT_T1_score, EMT_T2_score, EMT_delta,
                           Stemness_score, EMT_type, Stemness_class)
)), "\n")
log_cat("\nEMT-type counts:\n")
log_cat(capture.output(print(table(scores$EMT_type))), "\n")
log_cat("\nStemness counts:\n")
log_cat(capture.output(print(table(scores$Stemness_class))), "\n\n")


# ── 6. Plots ──────────────────────────────────────────────────────────────────

log_msg(strrep("─", 60))
log_msg("Generating plots...")

EMT_TYPE_COLS <- c(EMT_T1 = "#D73027", EMT_T2 = "#4575B4",
                   Mixed  = "#984EA3", No_EMT = "#A6A6A6")
STEM_COLS     <- c(Stemness_High = "#1B7837", Stemness_Low = "#A6D96A")
PANEL_COLS    <- c("T1 Positive" = "#E41A1C", "T1 Epithelial" = "#4DAF4A",
                   "T2 Tumor"   = "#377EB8", "T2 Stroma"    = "#A6CEE3",
                   "Stemness"   = "#984EA3")

# ── Plot A: Gene H-Score heatmap ──────────────────────────────────────────────
log_msg("  Plot A: gene H-Score heatmap...")

gene_order_5k <- c(
  intersect(EMT_T1_POS,   unique(hscores_5k$gene)),
  intersect(EMT_T1_EPI,   unique(hscores_5k$gene)),
  intersect(setdiff(EMT_T2_TUMOR, EMT_T1_POS), unique(hscores_5k$gene)),
  intersect(setdiff(EMT_T2_STROMA, MT_GENES),  unique(hscores_5k$gene)),
  intersect(STEMNESS_GENES, unique(hscores_5k$gene))
)
gene_order_mt <- intersect(MT_GENES, unique(hscores_mt$gene))
gene_order    <- c(gene_order_5k, gene_order_mt)

cat_map <- c(
  setNames(rep("T1 Positive",   length(intersect(EMT_T1_POS, gene_order))),
           intersect(EMT_T1_POS, gene_order)),
  setNames(rep("T1 Epithelial", length(intersect(EMT_T1_EPI, gene_order))),
           intersect(EMT_T1_EPI, gene_order)),
  setNames(rep("T2 Tumor",      length(intersect(setdiff(EMT_T2_TUMOR, EMT_T1_POS), gene_order))),
           intersect(setdiff(EMT_T2_TUMOR, EMT_T1_POS), gene_order)),
  setNames(rep("T2 Stroma",     length(intersect(EMT_T2_STROMA, gene_order))),
           intersect(EMT_T2_STROMA, gene_order)),
  setNames(rep("Stemness",      length(intersect(STEMNESS_GENES, gene_order))),
           intersect(STEMNESS_GENES, gene_order))
)

# Build z-score heatmap matrix (5K tumour genes + MT genes side by side)
z_hm_5k <- hscores_5k |>
  filter(gene %in% gene_order_5k) |>
  dplyr::select(sample_name, gene, H_Score) |>
  group_by(gene) |>
  mutate(z = { mu <- mean(H_Score); sg <- sd(H_Score)
               if (sg > 1e-9) (H_Score - mu) / sg else 0 }) |>
  ungroup()

z_hm_mt <- hscores_mt |>
  filter(gene %in% gene_order_mt) |>
  dplyr::select(sample_name, gene, H_Score) |>
  group_by(gene) |>
  mutate(z = { mu <- mean(H_Score); sg <- sd(H_Score)
               if (sg > 1e-9) (H_Score - mu) / sg else 0 }) |>
  ungroup()

z_hm_all <- dplyr::bind_rows(z_hm_5k, z_hm_mt)

z_wide_hm <- z_hm_all |>
  filter(gene %in% gene_order) |>
  dplyr::select(sample_name, gene, z) |>
  pivot_wider(names_from = gene, values_from = z, values_fill = 0)

gene_order <- intersect(gene_order, colnames(z_wide_hm))
mat_hm     <- as.matrix(z_wide_hm[, gene_order])
rownames(mat_hm) <- z_wide_hm$sample_name
mat_hm <- mat_hm[intersect(
  scores |> arrange(desc(EMT_T1_score)) |> pull(sample_name),
  rownames(mat_hm)), , drop = FALSE]

ann_col <- data.frame(Category = factor(cat_map[gene_order], levels = unique(cat_map)),
                      row.names = gene_order)
ann_row <- data.frame(
  EMT_type = scores$EMT_type[match(rownames(mat_hm), scores$sample_name)],
  Stemness = scores$Stemness_class[match(rownames(mat_hm), scores$sample_name)],
  row.names = rownames(mat_hm)
)

# Mark MT-panel genes with asterisk in column labels
col_labels           <- gene_order
col_labels[gene_order %in% MT_GENES] <- paste0(gene_order[gene_order %in% MT_GENES], "*")

png(file.path(PLOT_DIR, "emt_gene_heatmap_IC.png"),
    width = 18, height = 5, units = "in", res = 200, bg = "white")
pheatmap(
  mat_hm, cluster_rows = FALSE, cluster_cols = FALSE,
  labels_col    = col_labels,
  annotation_col = ann_col, annotation_row = ann_row,
  annotation_colors = list(Category = PANEL_COLS, EMT_type = EMT_TYPE_COLS,
                           Stemness = STEM_COLS),
  color      = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks     = seq(-3, 3, length.out = 101),
  border_color = "grey80", fontsize = 9, fontsize_col = 7, fontsize_row = 9,
  main   = "EMT & Stemness Gene Expression — IC cases (z-scored H-Score)\n* = MultiTissue panel gene",
  silent = TRUE
)
dev.off()
log_msg("  Plot A saved.")


# ── Plot B: EMT scatter ────────────────────────────────────────────────────────
log_msg("  Plot B: EMT scatter...")

p_scatter <- ggplot(scores,
  aes(x = EMT_T1_score, y = EMT_T2_score,
      colour = EMT_type, size = pmax(Stemness_score + 2, 0.2),
      label = sample_name)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_point(alpha = 0.85, stroke = 0.5) +
  ggrepel::geom_text_repel(size = 3.5, max.overlaps = 20, show.legend = FALSE,
                            segment.size = 0.3) +
  scale_colour_manual(values = EMT_TYPE_COLS, name = "EMT type") +
  scale_size_continuous(name = "Stemness score\n(shifted +2)", range = c(2, 8)) +
  labs(title    = "IC Cohort — EMT Subtype Classification",
       subtitle = "Scores z-scaled across IC samples only | * ALDH1A3 & LY6D from MultiTissue panel",
       x = "EMT-T1 score  (Invasive / Embryonic)",
       y = "EMT-T2 score  (Inflammatory / Adult)") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = 8),
        legend.position = "right")

ggsave(file.path(PLOT_DIR, "emt_scatter_IC.png"),
       p_scatter, width = 8, height = 6, dpi = 200, bg = "white")
log_msg("  Plot B saved.")


# ── Plot C: Score bar chart ────────────────────────────────────────────────────
log_msg("  Plot C: score bar chart...")

scores_long <- scores |>
  dplyr::select(sample_name, EMT_T1_score, EMT_T2_score, Stemness_score, EMT_type) |>
  pivot_longer(cols = c(EMT_T1_score, EMT_T2_score, Stemness_score),
               names_to = "score_type", values_to = "score") |>
  mutate(
    score_type = recode(score_type,
      EMT_T1_score   = "EMT-T1 (Invasive)",
      EMT_T2_score   = "EMT-T2 (Inflammatory)",
      Stemness_score = "Stemness"),
    score_type  = factor(score_type, levels = c("EMT-T1 (Invasive)",
                                                 "EMT-T2 (Inflammatory)", "Stemness")),
    sample_name = factor(sample_name,
      levels = scores |> arrange(desc(EMT_T1_score)) |> pull(sample_name))
  )

p_bars <- ggplot(scores_long, aes(x = sample_name, y = score, fill = score_type)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
  scale_fill_manual(values = c("EMT-T1 (Invasive)"     = "#D73027",
                                "EMT-T2 (Inflammatory)" = "#4575B4",
                                "Stemness"              = "#984EA3"), name = NULL) +
  labs(title    = "EMT & Stemness Scores — IC Cohort",
       subtitle = "Z-scored H-Score | ALDH1A3 & LY6D added from MultiTissue panel",
       x = NULL, y = "Score (z-scored)") +
  theme_classic(base_size = 11) +
  theme(plot.title  = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "top")

ggsave(file.path(PLOT_DIR, "emt_score_bars_IC.png"),
       p_bars, width = 8, height = 5, dpi = 200, bg = "white")
log_msg("  Plot C saved.")


# ── Plot D: Stemness bar chart ────────────────────────────────────────────────
log_msg("  Plot D: stemness bar chart...")

stem_avail <- intersect(STEMNESS_GENES, colnames(z_combined))
p_stem <- scores |>
  mutate(sample_name = factor(sample_name,
    levels = arrange(scores, desc(Stemness_score))$sample_name)) |>
  ggplot(aes(x = sample_name, y = Stemness_score, fill = Stemness_class)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey30") +
  scale_fill_manual(values = STEM_COLS, name = NULL) +
  labs(title    = "Stemness Score — IC Cohort",
       subtitle = paste0("Mean z-scored H-Score: ", paste(stem_avail, collapse = ", ")),
       x = NULL, y = "Stemness score (z-scored)") +
  theme_classic(base_size = 11) +
  theme(plot.title  = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "top")

ggsave(file.path(PLOT_DIR, "stemness_score_bars_IC.png"),
       p_stem, width = 7, height = 5, dpi = 200, bg = "white")
log_msg("  Plot D saved.")


# ── Plot E: Per-sample spatial panels (5 facets) ─────────────────────────────
# H&E | Tumour/Stroma | EMT-T1 | EMT-T2 | Stemness
# Per-cell scores from 5K BPCells (ALDH1A3/LY6D from MT are in a different
# coordinate space so are not overlaid on the per-cell scatter)
log_msg(strrep("─", 60))
log_msg("  Plot E: 5-facet spatial panels (IC only)...")

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

get_he_thumb <- function(sname, source_path, thumb_dir) {
  if (grepl("\\.png$", source_path, ignore.case = TRUE))
    return(if (file.exists(source_path)) source_path else NULL)
  thumb_png <- file.path(thumb_dir, paste0(sname, "_he_thumb.png"))
  if (file.exists(thumb_png)) return(thumb_png)
  if (!file.exists(source_path)) return(NULL)
  log_msg(sprintf("    Generating H&E thumbnail for %s ...", sname))
  cmd <- sprintf('convert "%s[0]" -thumbnail 1200x -normalize -colorspace sRGB "%s"',
                 source_path, thumb_png)
  ret <- system(cmd, intern = FALSE, ignore.stderr = TRUE)
  if (ret == 0 && file.exists(thumb_png)) thumb_png else NULL
}

ALL_EMT_SPATIAL_5K <- setdiff(
  unique(c(EMT_T1_POS, EMT_T1_EPI, EMT_T2_TUMOR, EMT_T2_STROMA, STEMNESS_GENES)),
  MT_GENES
)

dark_theme <- theme_void(base_size = 10) +
  theme(
    panel.background  = element_rect(fill = "black", colour = NA),
    plot.background   = element_rect(fill = "black", colour = NA),
    legend.background = element_rect(fill = "black", colour = NA),
    legend.key        = element_rect(fill = "black", colour = NA),
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 10,
                                     colour = "white"),
    legend.position   = "bottom",
    legend.text       = element_text(colour = "white", size = 8),
    legend.title      = element_text(colour = "white", size = 8)
  )

spatial_meta <- cohort_meta |> dplyr::select(cell_id, sample_name, x, y, Compartment_binary)

build_spatial_panels <- function(sname) {
  bpc_dir <- file.path(ANALYSIS_ROOT, "per_sample", sname, "bpcells", "counts")
  if (!dir.exists(bpc_dir)) return(NULL)

  bpc         <- open_matrix_dir(bpc_dir)
  genes_avail <- intersect(ALL_EMT_SPATIAL_5K, rownames(bpc))
  meta_s      <- spatial_meta |> filter(sample_name == sname)
  cells_use   <- intersect(meta_s$cell_id, colnames(bpc))
  if (length(cells_use) == 0 || length(genes_avail) == 0) return(NULL)

  expr_mat <- as.matrix(bpc[genes_avail, cells_use, drop = FALSE])

  mean_log1p <- function(genes) {
    g <- intersect(genes, rownames(expr_mat))
    if (length(g) == 0) return(rep(0, length(cells_use)))
    as.numeric(colMeans(log1p(expr_mat[g, , drop = FALSE])))
  }

  scores_df <- data.frame(
    cell_id  = cells_use,
    emt_t1   = mean_log1p(setdiff(EMT_T1_POS, MT_GENES)) -
               mean_log1p(EMT_T1_EPI),
    emt_t2   = mean_log1p(unique(c(setdiff(EMT_T2_TUMOR,  MT_GENES),
                                   setdiff(EMT_T2_STROMA, MT_GENES)))),
    stemness = mean_log1p(STEMNESS_GENES),
    stringsAsFactors = FALSE
  )

  plot_df <- meta_s |>
    filter(cell_id %in% cells_use) |>
    left_join(scores_df, by = "cell_id") |>
    filter(!is.na(x), !is.na(y))
  if (nrow(plot_df) == 0) return(NULL)

  pt <- if (nrow(plot_df) > 100000) 0.08 else if (nrow(plot_df) > 50000) 0.15 else 0.3

  # Panel 1 — H&E
  src   <- HE_SOURCE_MAP[sname]
  thumb <- if (!is.na(src)) get_he_thumb(sname, src, HE_THUMB_DIR) else NULL
  if (!is.null(thumb) && file.exists(thumb)) {
    img <- png::readPNG(thumb)
    if (length(dim(img)) == 3) img <- img[nrow(img):1, , , drop = FALSE]
    else                        img <- img[nrow(img):1,  ,  drop = FALSE]
    if (length(dim(img)) == 3) {
      wmask <- img[,,1] > 0.92 & img[,,2] > 0.92 & img[,,3] > 0.92
      img[,,1][wmask] <- 0; img[,,2][wmask] <- 0; img[,,3][wmask] <- 0
    }
    p1 <- ggplot() +
      annotation_raster(as.raster(img), xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
      xlim(0, 1) + ylim(0, 1) +
      coord_fixed(ratio = nrow(img) / ncol(img)) +
      labs(title = "H&E") + dark_theme
  } else {
    p1 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "H&E\n(not available)",
               colour = "grey60", size = 4, hjust = 0.5) +
      xlim(0, 1) + ylim(0, 1) + coord_fixed() +
      labs(title = "H&E") + dark_theme
  }

  # Panel 2 — Tumour / Stroma
  p2 <- ggplot(plot_df, aes(x = x, y = y, colour = Compartment_binary)) +
    geom_point(size = pt, stroke = 0, alpha = 0.8) +
    scale_colour_manual(values = c(Tumor = "#FF6B6B", Stroma = "#4CC9F0"), name = NULL) +
    coord_equal() + labs(title = "Tumour / Stroma") + dark_theme +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

  # Panel 3 — EMT-T1 (diverging)
  t1_abs <- max(abs(quantile(plot_df$emt_t1, c(0.02, 0.98), na.rm = TRUE)), 0.01)
  p3 <- ggplot(plot_df |> arrange(abs(emt_t1)), aes(x = x, y = y, colour = emt_t1)) +
    geom_point(size = pt, stroke = 0, alpha = 0.85) +
    scale_colour_gradientn(
      colours = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      limits = c(-t1_abs, t1_abs), oob = scales::squish,
      name = "EMT-T1\n(log1p)",
      guide = guide_colourbar(barwidth = 6, barheight = 0.5,
                              title.position = "top", title.hjust = 0.5)
    ) + coord_equal() + labs(title = "EMT-T1 (Invasive)") + dark_theme

  # Panel 4 — EMT-T2 (plasma)
  t2_q95 <- max(quantile(plot_df$emt_t2, 0.95, na.rm = TRUE), 0.01)
  p4 <- ggplot(plot_df |> arrange(emt_t2),
               aes(x = x, y = y, colour = pmin(emt_t2, t2_q95))) +
    geom_point(size = pt, stroke = 0, alpha = 0.85) +
    scale_colour_viridis_c(option = "plasma", name = "EMT-T2\n(log1p)",
      limits = c(0, t2_q95), begin = 0.05,
      guide = guide_colourbar(barwidth = 6, barheight = 0.5,
                              title.position = "top", title.hjust = 0.5)
    ) + coord_equal() + labs(title = "EMT-T2 (Inflammatory)") + dark_theme

  # Panel 5 — Stemness (magma)
  stem_q95 <- max(quantile(plot_df$stemness, 0.95, na.rm = TRUE), 0.01)
  p5 <- ggplot(
    plot_df |> arrange(stemness) |> mutate(stemness = pmin(stemness, stem_q95)),
    aes(x = x, y = y, colour = stemness)
  ) +
    geom_point(size = pt, stroke = 0, alpha = 0.85) +
    scale_colour_viridis_c(option = "magma", name = "Stemness\n(log1p)",
      limits = c(0, stem_q95), begin = 0.05,
      guide = guide_colourbar(barwidth = 6, barheight = 0.5,
                              title.position = "top", title.hjust = 0.5)
    ) + coord_equal() + labs(title = "Stemness") + dark_theme

  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5)
}

# Per-sample individual plots
for (sname in IC_SAMPLES) {
  log_msg(sprintf("    %s ...", sname))
  panels <- build_spatial_panels(sname)
  if (is.null(panels)) next

  p1_labelled <- panels$p1 +
    labs(y = gsub("_", "-", sname)) +
    theme(axis.title.y = element_text(colour = "white", size = 11, face = "bold",
                                      angle = 90, vjust = 0.5, margin = margin(r = 6)))

  fig <- (p1_labelled + panels$p2 + panels$p3 + panels$p4 + panels$p5) +
    patchwork::plot_layout(nrow = 1) +
    patchwork::plot_annotation(
      theme = theme(plot.background = element_rect(fill = "black", colour = NA))
    ) & theme(plot.background = element_rect(fill = "black", colour = NA))

  ggsave(file.path(SPATIAL_DIR, paste0(sname, "_spatial_emt_stemness_IC.png")),
         fig, width = 25, height = 6, dpi = 200, bg = "black")
  log_msg(sprintf("    Saved: %s", paste0(sname, "_spatial_emt_stemness_IC.png")))
}

# Consolidated IC plot (no H&E, sample name on y-axis)
log_msg("  Building consolidated IC spatial plot (no H&E)...")
ic_rows <- lapply(IC_SAMPLES, function(sname) {
  panels <- build_spatial_panels(sname)
  if (is.null(panels)) return(NULL)
  p2_labelled <- panels$p2 +
    labs(y = gsub("_", "-", sname)) +
    theme(axis.title.y = element_text(colour = "white", size = 11, face = "bold",
                                      angle = 90, vjust = 0.5, margin = margin(r = 6)))
  (p2_labelled + panels$p3 + panels$p4 + panels$p5) +
    patchwork::plot_layout(nrow = 1) +
    patchwork::plot_annotation(
      theme = theme(plot.background = element_rect(fill = "black", colour = NA))
    ) & theme(plot.background = element_rect(fill = "black", colour = NA))
})
ic_rows <- Filter(Negate(is.null), ic_rows)

if (length(ic_rows) > 0) {
  ic_combined <- patchwork::wrap_plots(ic_rows, ncol = 1) +
    patchwork::plot_annotation(
      title   = "IC Cohort  —  EMT & Stemness Spatial Overview",
      caption = "Tumour/Stroma (CF3)  |  EMT-T1 (RdBu)  |  EMT-T2 (plasma)  |  Stemness (magma)  |  Per-cell scores from 5K panel",
      theme   = theme(
        plot.title   = element_text(face = "bold", hjust = 0.5, size = 16, colour = "white"),
        plot.caption = element_text(hjust = 0.5, size = 9, colour = "grey60"),
        plot.background = element_rect(fill = "black", colour = NA)
      )
    ) & theme(plot.background = element_rect(fill = "black", colour = NA))

  ggsave(file.path(PLOT_DIR, "IC_combined_spatial.png"), ic_combined,
         width = 20, height = length(ic_rows) * 6, dpi = 200, bg = "black")
  log_msg("  Consolidated IC spatial plot saved.")
}

log_msg("  Plot E done.")


# ── Final summary ─────────────────────────────────────────────────────────────

log_msg(strrep("=", 60))
log_msg("SCRIPT 12b COMPLETE")
log_msg(sprintf("5K-panel genes scored: %d", length(unique(hscores_5k$gene))))
log_msg(sprintf("MT-panel genes added:  %s", paste(MT_GENES, collapse = ", ")))
log_msg(sprintf("Outputs: %s", OUT_DIR))
log_msg("  emt_gene_hscores_IC.csv")
log_msg("  emt_stemness_classification_IC.csv")
log_msg(sprintf("  PLOTs/12b_emt_stemness_IC/"))
log_msg(strrep("=", 60))

close(.log_con)
