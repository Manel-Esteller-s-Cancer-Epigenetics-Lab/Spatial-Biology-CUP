# =============================================================================
# 03_threshold_classification.R
# -----------------------------------------------------------------------------
# PURPOSE : Define gene-level positivity thresholds and classify samples.
#
# PHASE 1 — Cell-level threshold (t_cell) per gene, three methods:
#   A. Negative Control Probes : 99th percentile of NegControlProbe
#      expression across all cells × all samples. Primary method —
#      directly controls false-discovery at the transcript level.
#   B. Stromal Floor           : 95th percentile of stromal cell expression.
#      Applied to epithelial/tumor markers as biological validation.
#   C. GMM (mclust)            : Two-component Gaussian mixture on tumor
#      cell expression. Applied to ICI markers where stroma is not a
#      valid negative control.
#
# PHASE 2 — Spatial H-Score per sample per gene:
#   pct_pos   = % tumor cells with expression > t_cell  (penetrance)
#   mean_int  = mean expression of positive tumor cells  (intensity)
#   H_Score   = pct_pos × mean_int                       (composite)
#
# PHASE 3 — Sample classification per gene:
#   1D K-means (k=2,3) on H_Score across samples.
#   Optimal k chosen by gap statistic.
#   Output: Negative / Low / High classification per sample per gene.
#
# INPUT  : cohort_slim.rds  (from 01b_assemble_cohort.R)
#          Raw Xenium cell_feature_matrix.h5 per sample (for NegControl)
#
# OUTPUT : <ANALYSIS_ROOT>/Classification/
#   ├── negcontrol_expression.rds   per-cell NegControl expression
#   ├── threshold_table.csv         t_cell per gene × method
#   ├── hscores_per_sample.csv      H-Score + pct_pos + mean_int per sample × gene
#   ├── sample_classification.csv   Negative/Low/High per sample × gene
#   ├── plots/
#   │   ├── negcontrol_distribution.png
#   │   ├── gmm_fits_ICI.png
#   │   ├── stromal_floor_epithelial.png
#   │   ├── hscore_heatmap.png
#   │   └── sample_classification_heatmap.png
#   └── 03_threshold_classification.log
#
# USAGE  : Rscript 03_threshold_classification.R
# =============================================================================


# ── 0. Configuration ──────────────────────────────────────────────────────────

ANALYSIS_ROOT <- "."
COHORT_RDS    <- file.path(ANALYSIS_ROOT, "cohort_slim_v2.rds")
OUT_DIR       <- file.path(ANALYSIS_ROOT, "Classification")
PLOT_DIR      <- file.path(ANALYSIS_ROOT, "PLOTs", "03_classification")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Sample sheet — needed to locate raw h5 files for NegControl extraction
SAMPLE_SHEET  <- file.path(ANALYSIS_ROOT, "sample_sheet.csv")

# Marker gene groups
GENES_EPITHELIAL <- c("KRT20", "MET", "ERBB2", "KIT", "FGFR2", "CCNE1", "MYC")
GENES_ICI_CPS    <- c("CD274")                              # PD-L1: Combined Positive Score (tumor + stroma)
GENES_ICI_STROMA <- c("LAG3", "CTLA4", "TIGIT", "PDCD1")   # T-cell exhaustion markers: stroma only
GENES_ICI        <- c(GENES_ICI_CPS, GENES_ICI_STROMA)
ALL_GENES        <- c(GENES_EPITHELIAL, GENES_ICI)

# Samples to exclude from analysis (not processed in cohort)

# NegControl probe name pattern (as seen in features.tsv.gz)
NEGCTRL_PATTERN  <- "^NegControlProbe"

# Thresholds
NEGCTRL_PCT      <- 0.99   # 99th percentile of NegControl distribution
STROMAL_PCT      <- 0.95   # 95th percentile of stromal expression
GMM_COMPONENTS   <- 2      # two-component mixture for ICI genes
GMM_MAX_CELLS    <- 100000 # subsample for GMM fitting — mclust is slow on >500k

# K-means classification
KMAX             <- 3      # maximum k to test (gap statistic chooses 2 or 3)


# ── 1. Libraries ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(mclust)    # GMM fitting    — install.packages("mclust")
  library(cluster)   # gap statistic  — install.packages("cluster")
  library(knitr)
})

options(Seurat.object.assay.version = "v5")

# Gene display name helper (mirrors Script 02)
GENE_DISPLAY <- c(
  CD274 = "PD-L1 (CD274)",
  PDCD1 = "PD-1 (PDCD1)",
  KRT20 = "CK20 (KRT20)",
  ERBB2 = "ERBB2 (HER2)"
)
label_gene <- function(g) ifelse(g %in% names(GENE_DISPLAY), GENE_DISPLAY[g], g)


# ── 2. Logging ────────────────────────────────────────────────────────────────

log_path <- file.path(OUT_DIR, "03_threshold_classification.log")
.log_con <- file(log_path, open = "wt")

log_msg <- function(...) {
  txt <- paste0(format(Sys.time(), "[%H:%M:%S] "), paste(..., sep = ""))
  message(txt)
  tryCatch(writeLines(txt, con = .log_con), error = function(e) NULL)
}
log_cat <- function(...) {
  txt <- paste(..., sep = "")
  cat(txt)
  tryCatch(writeLines(trimws(txt,"right"), con = .log_con),
           error = function(e) NULL)
}

log_msg("Script 03 started")
log_msg(sprintf("Output directory: %s", OUT_DIR))


# ── 3. Load cohort slim ───────────────────────────────────────────────────────

log_msg("Loading cohort_slim.rds...")
cohort <- readRDS(COHORT_RDS)

cohort_tumor  <- cohort |> filter(Compartment_binary == "Tumor")
cohort_stroma <- cohort |> filter(Compartment_binary == "Stroma")

log_cat(sprintf(
  "Cohort: %d samples | %d total cells | %d tumor | %d stroma\n",
  length(unique(cohort$sample_name)),
  nrow(cohort),
  nrow(cohort_tumor),
  nrow(cohort_stroma)
))

samples <- read.csv(SAMPLE_SHEET, stringsAsFactors = FALSE) 



# =============================================================================
# PHASE 1A — Negative Control Threshold
# Extract NegControlProbe expression from raw Xenium h5 files,
# compute 99th percentile across all cells × all samples.
# =============================================================================

log_msg(strrep("─", 60))
log_msg("PHASE 1A: Negative Control Probe extraction")
log_msg(strrep("─", 60))

negctrl_rds <- file.path(OUT_DIR, "negcontrol_expression.rds")

if (file.exists(negctrl_rds)) {
  log_msg("Checkpoint found — loading negcontrol_expression.rds")
  negctrl_all <- readRDS(negctrl_rds)
} else {

  negctrl_list <- vector("list", nrow(samples))

  for (i in seq_len(nrow(samples))) {
    sname    <- samples$sample_name[i]
    data_dir <- samples$data_dir[i]
    h5_path  <- file.path(data_dir, "cell_feature_matrix.h5")

    log_msg(sprintf("  [%d/%d] %s", i, nrow(samples), sname))
    log_cat(sprintf("    h5 path: %s\n", h5_path))

    if (!file.exists(h5_path)) {
      log_msg(sprintf("    h5 NOT FOUND — skipping (check drive is mounted)"))
      next
    }

    log_msg(sprintf("    h5 found — reading..."))

    tryCatch({
      # Read10X_h5 returns a named list when multiple feature types present
      mat_list <- Seurat::Read10X_h5(
        filename        = h5_path,
        use.names       = TRUE,
        unique.features = TRUE
      )

      log_cat(sprintf("    mat_list class: %s\n", class(mat_list)))

      if (is.list(mat_list)) {
        log_cat(sprintf("    Feature types in h5: %s\n",
                        paste(names(mat_list), collapse = " | ")))
        neg_key <- grep("Negative Control", names(mat_list),
                        value = TRUE, ignore.case = TRUE)
        if (length(neg_key) == 0) {
          log_msg(sprintf("    No 'Negative Control' key — skipping"))
          next
        }
        neg_mat <- mat_list[[neg_key[1]]]
        rm(mat_list)
      } else {
        # Single matrix — filter rows matching NegControlProbe pattern
        log_cat(sprintf("    Single matrix: %d features\n", nrow(mat_list)))
        neg_rows <- grepl(NEGCTRL_PATTERN, rownames(mat_list))
        log_cat(sprintf("    NegControl rows found: %d\n", sum(neg_rows)))
        if (sum(neg_rows) == 0) {
          log_msg("    No NegControl rows — skipping")
          next
        }
        neg_mat <- mat_list[neg_rows, , drop = FALSE]
        rm(mat_list)
      }

      gc(verbose = FALSE)
      log_cat(sprintf("    NegControl probes: %d | Cells: %d\n",
                      nrow(neg_mat), ncol(neg_mat)))

      neg_sum <- Matrix::colSums(neg_mat)

      negctrl_list[[i]] <- data.frame(
        sample_name  = sname,
        cell_id      = names(neg_sum),
        neg_sum_raw  = as.numeric(neg_sum),
        n_probes     = nrow(neg_mat),
        stringsAsFactors = FALSE
      )

      rm(neg_mat, neg_sum)
      gc(verbose = FALSE)

    }, error = function(e) {
      log_msg(sprintf("    ERROR: %s", conditionMessage(e)))
    })
  }

  negctrl_all <- dplyr::bind_rows(Filter(Negate(is.null), negctrl_list))

  if (nrow(negctrl_all) > 0) {
    saveRDS(negctrl_all, negctrl_rds)
    log_msg(sprintf("NegControl data saved: %d cells × %d samples",
                    nrow(negctrl_all),
                    length(unique(negctrl_all$sample_name))))
  } else {
    log_msg("NegControl extraction returned 0 cells — checkpoint NOT saved.")
    log_msg("Ensure h5 data files are accessible.")
  }
}

# Compute 99th percentile of per-probe-normalised background
# neg_sum / n_probes gives mean counts per NegControl probe per cell
# This is the per-probe background level — threshold applied per gene
# Guard against empty extraction (e.g. drive not mounted)
if (nrow(negctrl_all) == 0) {
  log_msg("WARNING: No NegControl data extracted — check that raw Xenium h5 files")
  log_msg("  are accessible. Falling back to mean+2SD threshold = log1p(1).")
  negctrl_threshold_raw <- 1
  negctrl_threshold     <- log1p(1)
} else {

negctrl_all <- negctrl_all |>
  mutate(neg_per_probe = neg_sum_raw / n_probes)

# 99th percentile of raw per-probe counts across all cells × all samples
negctrl_threshold_raw <- quantile(negctrl_all$neg_per_probe,
                                   NEGCTRL_PCT, na.rm = TRUE)

# Convert to log-normalised scale: log1p(count / n_probes * scale_factor)
# We use scale_factor = 1 (counts are already per-probe) and match the
# log1p(x * 10000 / lib_size) normalisation used by NormalizeData.
# Since NegControl probes are not in the library size calculation,
# we use log1p(negctrl_threshold_raw) as a conservative direct threshold.
# This is appropriate because NegControl counts are already very low
# and log1p(0) = 0 means any expression > 0 clears this bar — so we
# use the raw count directly as the per-probe noise floor.
negctrl_threshold <- log1p(negctrl_threshold_raw)

# If 99th pctile is 0 (most cells have zero NegControl counts — expected for
# Xenium), use mean + 2SD of per-probe counts as the threshold instead.
# This is still derived from NegControl biology, not an arbitrary number.
if (is.nan(negctrl_threshold) || negctrl_threshold == 0) {
  neg_mean <- mean(negctrl_all$neg_per_probe, na.rm = TRUE)
  neg_sd   <- sd(negctrl_all$neg_per_probe,   na.rm = TRUE)
  negctrl_threshold_raw <- neg_mean + 2 * neg_sd
  negctrl_threshold     <- log1p(negctrl_threshold_raw)
  log_msg(sprintf(
    "99th pctile = 0 (expected for sparse data) — using mean+2SD: raw=%.4f log1p=%.4f",
    negctrl_threshold_raw, negctrl_threshold
  ))
}

# Final safety floor — if still 0 after mean+2SD (all NegControl = 0),
# fall back to log1p(1) = 0.693 which corresponds to ~1 raw count
if (is.nan(negctrl_threshold) || negctrl_threshold == 0) {
  negctrl_threshold <- log1p(1)
  log_msg("All NegControl counts are 0 — using log1p(1) = 0.693 as floor")
}

log_cat(sprintf(
  "NegControl threshold: raw=%.4f  log-norm=%.4f\n",
  negctrl_threshold_raw, negctrl_threshold
))

# Plot negative control distribution
p_negctrl <- ggplot(
  negctrl_all |> filter(neg_per_probe <= quantile(neg_per_probe, 0.999)),
  aes(x = neg_per_probe)
) +
  geom_histogram(bins = 80, fill = "#4393C3", alpha = 0.8, color = "white") +
  geom_vline(xintercept = negctrl_threshold_raw,
             color = "#B22222", linewidth = 1, linetype = "dashed") +
  annotate("text",
           x = negctrl_threshold_raw * 1.5,
           y = Inf, vjust = 1.5,
           label = sprintf("99th pctile = %.3f", negctrl_threshold_raw),
           color = "#B22222", size = 3.5, fontface = "bold") +
  facet_wrap(~ sample_name, scales = "free_y") +
  labs(
    title    = "Negative Control Probe Expression — Per-Probe Counts",
    subtitle = "Red dashed line = 99th percentile threshold (t_cell anchor)",
    x        = "Mean counts per NegControl probe per cell",
    y        = "Cell count"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    strip.text  = element_text(face = "bold", size = 7)
  )

  ggsave(file.path(PLOT_DIR, "negcontrol_distribution.png"),
         p_negctrl, width = 16, height = 10, dpi = 300, bg = "white")
  log_msg("NegControl distribution plot saved.")

} # end else (negctrl_all not empty)


# =============================================================================
# PHASE 1B — Stromal Floor Threshold (epithelial markers)
# =============================================================================

log_msg(strrep("─", 60))
log_msg("PHASE 1B: Stromal floor thresholds (epithelial markers)")
log_msg(strrep("─", 60))

stromal_thresholds <- cohort_stroma |>
  dplyr::select(all_of(GENES_EPITHELIAL)) |>
  summarise(across(everything(),
                   list(
                     p95   = ~ quantile(.x, STROMAL_PCT, na.rm = TRUE),
                     mean  = ~ mean(.x, na.rm = TRUE),
                     sd    = ~ sd(.x, na.rm = TRUE)
                   ))) |>
  pivot_longer(everything(),
               names_to  = c("gene", ".value"),
               names_pattern = "(.+)_(p95|mean|sd)") |>
  rename(stromal_floor = p95)

log_cat("Stromal floor thresholds (95th pctile of stroma):\n")
log_cat(capture.output(print(stromal_thresholds)), "\n")

# Plot stromal floor per gene
stroma_long <- cohort_stroma |>
  dplyr::select(all_of(GENES_EPITHELIAL)) |>
  pivot_longer(everything(), names_to = "gene", values_to = "expression") |>
  filter(!is.na(expression))

p_stromal <- ggplot(stroma_long,
    aes(x = expression)) +
  geom_density(fill = "#4393C3", alpha = 0.6, color = NA) +
  geom_vline(
    data    = stromal_thresholds,
    mapping = aes(xintercept = stromal_floor),
    color   = "#B22222", linewidth = 0.8, linetype = "dashed"
  ) +
  geom_vline(xintercept = negctrl_threshold,
             color = "darkgreen", linewidth = 0.8, linetype = "dotted") +
  facet_wrap(~ gene, scales = "free", ncol = 4) +
  labs(
    title    = "Stromal Expression Distribution — Epithelial Markers",
    subtitle = "Red dashed = stromal floor (95th pctile) | Green dotted = NegControl threshold",
    x        = "Log-normalised expression",
    y        = "Density"
  ) +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold"))

ggsave(file.path(PLOT_DIR, "stromal_floor_epithelial.png"),
       p_stromal, width = 14, height = 6, dpi = 300, bg = "white")
log_msg("Stromal floor plot saved.")


# =============================================================================
# PHASE 1C — GMM Threshold (ICI markers)
# Fit two-component Gaussian mixture to tumor cell expression.
# Threshold = intersection of the two components.
# =============================================================================

log_msg(strrep("─", 60))
log_msg("PHASE 1C: GMM thresholds (ICI markers)")
log_msg(strrep("─", 60))

gmm_thresholds <- list()
gmm_plot_list  <- list()

for (g in GENES_ICI) {

  # Source compartment: CD274 pools tumor + stroma (single GMM for CPS);
  # LAG3/CTLA4/TIGIT/PDCD1 use stroma only (T-cell exhaustion markers).
  if (g %in% GENES_ICI_CPS) {
    expr_source     <- c(cohort_tumor[[g]], cohort_stroma[[g]])
    gmm_compartment <- "tumor + stroma (CPS)"
  } else {
    expr_source     <- cohort_stroma[[g]]
    gmm_compartment <- "stroma"
  }

  # Filter to expressing cells (> 0) before GMM — fitting on all cells
  # (90%+ zeros) gives a degenerate model.
  expr_vals_all <- expr_source[!is.na(expr_source) &
                                 is.finite(expr_source) &
                                 expr_source > 0]

  if (length(expr_vals_all) < 100) {
    log_msg(sprintf("  %s: too few expressing cells (%d) — using NegControl threshold",
                    g, length(expr_vals_all)))
    gmm_thresholds[[g]] <- negctrl_threshold
    next
  }

  # Subsample for GMM fitting — mclust is slow on >100k cells
  set.seed(42)
  if (length(expr_vals_all) > GMM_MAX_CELLS) {
    expr_vals <- sample(expr_vals_all, GMM_MAX_CELLS)
    log_msg(sprintf("  Fitting GMM for %s (%d expressing, subsampled to %d)...",
                    g, length(expr_vals_all), GMM_MAX_CELLS))
  } else {
    expr_vals <- expr_vals_all
    log_msg(sprintf("  Fitting GMM for %s (%d expressing cells, %s)...",
                    g, length(expr_vals), gmm_compartment))
  }

  # Fit GMM — mclust automatically selects best model via BIC
  # modelNames = "V" allows unequal variances (more realistic)
  set.seed(42)
  gmm_fit <- tryCatch(
    mclust::Mclust(expr_vals, G = GMM_COMPONENTS,
                   modelNames = "V", verbose = FALSE),
    error = function(e) {
      log_msg(sprintf("    GMM failed for %s: %s — using NegControl threshold",
                      g, conditionMessage(e)))
      NULL
    }
  )

  if (is.null(gmm_fit)) {
    gmm_thresholds[[g]] <- negctrl_threshold
    next
  }

  # Extract component parameters
  means <- gmm_fit$parameters$mean
  vars  <- gmm_fit$parameters$variance$sigmasq
  props <- gmm_fit$parameters$pro

  # Sort components by mean (component 1 = low/background, 2 = expressing)
  ord    <- order(means)
  means  <- means[ord]
  vars   <- vars[ord]
  props  <- props[ord]

  # Find intersection of the two Gaussian components analytically
  # Solve: props[1]*N(x|mu1,sigma1) = props[2]*N(x|mu2,sigma2)
  # Numerically find the crossover point between the two components
  x_seq    <- seq(min(expr_vals), max(expr_vals), length.out = 5000)
  dens1    <- props[1] * dnorm(x_seq, means[1], sqrt(vars[1]))
  dens2    <- props[2] * dnorm(x_seq, means[2], sqrt(vars[2]))
  cross_idx <- which(diff(sign(dens2 - dens1)) > 0)

  if (length(cross_idx) == 0) {
    log_msg(sprintf("    No intersection found for %s — using NegControl threshold", g))
    gmm_thresh <- negctrl_threshold
  } else {
    gmm_thresh <- x_seq[cross_idx[1]]
    log_cat(sprintf("    %s: GMM threshold = %.4f (means: %.3f, %.3f)\n",
                    g, gmm_thresh, means[1], means[2]))
  }

  # Use the more conservative threshold between GMM and NegControl
  final_thresh <- max(gmm_thresh, negctrl_threshold)
  log_cat(sprintf("    %s: final threshold = %.4f (GMM=%.4f, NegCtrl=%.4f)\n",
                  g, final_thresh, gmm_thresh, negctrl_threshold))
  gmm_thresholds[[g]] <- final_thresh

  # Build plot for this gene
  df_dens <- data.frame(x = x_seq, dens1 = dens1, dens2 = dens2,
                         total = dens1 + dens2)
  p_gmm <- ggplot() +
    geom_histogram(
      data    = data.frame(expr = expr_vals),
      mapping = aes(x = expr, y = after_stat(density)),
      bins    = 80, fill = "grey80", color = "white", alpha = 0.7
    ) +
    geom_line(data = df_dens, aes(x = x, y = dens1),
              color = "#4393C3", linewidth = 0.8) +
    geom_line(data = df_dens, aes(x = x, y = dens2),
              color = "#B22222", linewidth = 0.8) +
    geom_line(data = df_dens, aes(x = x, y = total),
              color = "black", linewidth = 0.6, linetype = "dashed") +
    geom_vline(xintercept = final_thresh,
               color = "darkgreen", linewidth = 1) +
    geom_vline(xintercept = negctrl_threshold,
               color = "orange", linewidth = 0.7, linetype = "dotted") +
    annotate("text", x = final_thresh, y = Inf, vjust = 1.5, hjust = -0.1,
             label = sprintf("t=%.3f", final_thresh),
             color = "darkgreen", size = 3, fontface = "bold") +
    labs(title    = label_gene(g),
         subtitle = "Blue=background | Red=expressing | Green=final threshold",
         x        = sprintf("Log-norm expression (%s)", gmm_compartment),
         y        = "Density") +
    theme_classic(base_size = 9) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 7))

  gmm_plot_list[[g]] <- p_gmm
}

if (length(gmm_plot_list) > 0) {
  p_gmm_all <- wrap_plots(gmm_plot_list, ncol = min(3, length(gmm_plot_list))) +
    plot_annotation(title = "GMM Fits — ICI Markers (CD274: tumor + stroma | LAG3/CTLA4/TIGIT/PDCD1: stroma)")
  ggsave(file.path(PLOT_DIR, "gmm_fits_ICI.png"),
         p_gmm_all, width = 14, height = 8, dpi = 300, bg = "white")
  log_msg("GMM fit plots saved.")
} else {
  log_msg("WARNING: No GMM plots produced — all genes fell back to NegControl threshold.")
}


# =============================================================================
# PHASE 1 — Compile threshold table
# Primary: NegControl for all genes
# Validation: Stromal floor (epithelial) | GMM intersection (ICI)
# Final t_cell: max(NegControl, Stromal/GMM) — most conservative wins
# =============================================================================

log_msg(strrep("─", 60))
log_msg("Compiling threshold table...")
log_msg(strrep("─", 60))

threshold_table <- bind_rows(
  # Epithelial markers
  stromal_thresholds |>
    dplyr::select(gene, stromal_floor) |>
    mutate(
      negctrl_threshold = negctrl_threshold,
      gmm_threshold     = NA_real_,
      method_primary    = "NegControl",
      method_secondary  = "StromalFloor",
      t_cell = pmax(negctrl_threshold, stromal_floor)
    ),
  # ICI markers
  data.frame(
    gene             = GENES_ICI,
    stromal_floor    = NA_real_,
    negctrl_threshold = negctrl_threshold,
    gmm_threshold    = unlist(gmm_thresholds),
    method_primary   = "NegControl",
    method_secondary = "GMM",
    t_cell           = pmax(negctrl_threshold, unlist(gmm_thresholds)),
    stringsAsFactors = FALSE
  )
) |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))

write.csv(threshold_table,
          file.path(OUT_DIR, "threshold_table.csv"),
          row.names = FALSE)

log_cat("Final thresholds (t_cell):\n")
log_cat(capture.output(
  print(threshold_table |> dplyr::select(gene, t_cell, method_secondary))
), "\n")


# =============================================================================
# PHASE 2 — Spatial H-Score per sample per gene
# Computed on TUMOR cells only.
# pct_pos  = % tumor cells with expression > t_cell
# mean_int = mean expression of positive tumor cells
# H_Score  = pct_pos × mean_int
# =============================================================================

log_msg(strrep("─", 60))
log_msg("PHASE 2: Computing Spatial H-Scores...")
log_msg(strrep("─", 60))

# Build a named vector of thresholds for easy lookup
t_cell_vec <- setNames(threshold_table$t_cell, threshold_table$gene)

hscore_rows <- list()

for (g in ALL_GENES) {
  t <- t_cell_vec[g]

  if (g %in% GENES_ICI_CPS) {
    # CD274: Combined Positive Score — positives from tumor + stroma,
    # denominator = tumor cells (matches clinical CPS definition).
    tumor_expr  <- cohort_tumor  |> dplyr::select(sample_name, expr = all_of(g)) |> filter(!is.na(expr))
    stroma_expr <- cohort_stroma |> dplyr::select(sample_name, expr = all_of(g)) |> filter(!is.na(expr))

    denom <- tumor_expr |>
      group_by(sample_name) |>
      summarise(n_cells_denom = n(), .groups = "drop")

    pos_both <- dplyr::bind_rows(tumor_expr, stroma_expr) |>
      filter(expr > t) |>
      group_by(sample_name) |>
      summarise(n_positive = n(), mean_int = mean(expr), .groups = "drop")

    gene_data <- denom |>
      left_join(pos_both, by = "sample_name") |>
      mutate(
        n_positive  = coalesce(n_positive, 0L),
        mean_int    = coalesce(mean_int, 0),
        pct_pos     = 100 * n_positive / n_cells_denom,
        H_Score     = pct_pos * mean_int,
        t_cell_used = t,
        compartment = "CPS",
        gene        = g
      )

  } else if (g %in% GENES_ICI_STROMA) {
    # LAG3, CTLA4, TIGIT, PDCD1: stroma cells only.
    gene_data <- cohort_stroma |>
      dplyr::select(sample_name, expr = all_of(g)) |>
      filter(!is.na(expr)) |>
      group_by(sample_name) |>
      summarise(
        n_cells_denom = n(),
        n_positive    = sum(expr > t),
        pct_pos       = 100 * n_positive / n_cells_denom,
        mean_int      = ifelse(n_positive > 0, mean(expr[expr > t]), 0),
        H_Score       = pct_pos * mean_int,
        t_cell_used   = t,
        compartment   = "Stroma",
        .groups       = "drop"
      ) |>
      mutate(gene = g)

  } else {
    # Epithelial markers: tumor cells only.
    gene_data <- cohort_tumor |>
      dplyr::select(sample_name, expr = all_of(g)) |>
      filter(!is.na(expr)) |>
      group_by(sample_name) |>
      summarise(
        n_cells_denom = n(),
        n_positive    = sum(expr > t),
        pct_pos       = 100 * n_positive / n_cells_denom,
        mean_int      = ifelse(n_positive > 0, mean(expr[expr > t]), 0),
        H_Score       = pct_pos * mean_int,
        t_cell_used   = t,
        compartment   = "Tumor",
        .groups       = "drop"
      ) |>
      mutate(gene = g)
  }

  hscore_rows[[g]] <- gene_data
}

hscores <- dplyr::bind_rows(hscore_rows) |>
  dplyr::select(sample_name, gene, compartment, n_cells_denom, n_positive,
                pct_pos, mean_int, H_Score, t_cell_used) |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))

write.csv(hscores,
          file.path(OUT_DIR, "hscores_per_sample.csv"),
          row.names = FALSE)

log_msg(sprintf("H-Scores computed: %d rows", nrow(hscores)))

# ── H-Score heatmap ───────────────────────────────────────────────────────────
hscore_wide <- hscores |>
  dplyr::select(sample_name, gene, H_Score) |>
  pivot_wider(names_from = gene, values_from = H_Score)

p_hscore_heat <- hscores |>
  mutate(gene = factor(gene, levels = ALL_GENES)) |>
  ggplot(aes(x = gene, y = sample_name, fill = H_Score)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = round(H_Score, 1)),
            size = 2.5, color = "grey10", fontface = "bold") +
  scale_fill_gradientn(
    name   = "H-Score\n(pct × mean)",
    colors = c("grey95", "#FDBB84", "#E34A33", "#7F0000"),
    limits = c(0, NA)
  ) +
  scale_x_discrete(labels = function(x) ifelse(x == "CD274", "PD-L1",
                                        ifelse(x == "PDCD1", "PD-1",
                                        ifelse(x == "ERBB2", "HER2", x)))) +
  labs(
    title    = "Spatial H-Score per Sample × Marker",
    subtitle = "H-Score = % cells > t_cell × mean expression of positive cells (compartment is gene-specific)",
    x        = NULL, y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8),
    axis.text.x   = element_text(angle = 35, hjust = 1, face = "bold"),
    axis.text.y   = element_text(size = 9)
  )

ggsave(file.path(PLOT_DIR, "hscore_heatmap.png"),
       p_hscore_heat, width = 12, height = 7, dpi = 300, bg = "white")
log_msg("H-Score heatmap saved.")


# =============================================================================
# PHASE 3 — Sample classification via K-means on H-Score
# For each gene: test k=2 and k=3, choose via gap statistic.
# Labels: k=2 → Negative / Positive
#         k=3 → Negative / Low / High
# =============================================================================

log_msg(strrep("─", 60))
log_msg("PHASE 3: K-means classification...")
log_msg(strrep("─", 60))

classify_gene <- function(gene_name, hscore_df, kmax = KMAX) {

  dat <- hscore_df |>
    filter(gene == gene_name) |>
    arrange(sample_name)

  x      <- matrix(dat$H_Score, ncol = 1)
  n_samp <- nrow(dat)

  if (n_samp < 4) {
    log_msg(sprintf("  %s: too few samples (%d) for gap statistic — using k=2",
                    gene_name, n_samp))
    best_k <- 2
  } else {
    # Gap statistic to choose optimal k
    set.seed(42)
    gap <- cluster::clusGap(x, FUN = kmeans, K.max = min(kmax, n_samp - 1),
                            B = 100, verbose = FALSE)
    best_k <- cluster::maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"],
                              method = "Tibs2001SEmax")
    best_k <- max(2, min(best_k, kmax))  # enforce 2 ≤ k ≤ kmax
    log_cat(sprintf("  %s: optimal k = %d (gap statistic)\n",
                    gene_name, best_k))
  }

  set.seed(42)
  km    <- kmeans(x, centers = best_k, nstart = 50, iter.max = 100)
  cents <- sort(km$centers[, 1])

  # Label clusters by H-Score rank (lowest = Negative)
  if (best_k == 2) {
    labels <- c("Negative", "Positive")
  } else {
    labels <- c("Negative", "Low", "High")
  }

  # Map cluster IDs to ordered labels
  cluster_order <- order(km$centers[, 1])
  label_map     <- setNames(labels[seq_along(cluster_order)],
                             cluster_order)
  dat$classification <- label_map[as.character(km$cluster)]
  dat$k_used         <- best_k
  dat$classification <- factor(dat$classification,
                                levels = c("Negative","Low","Positive","High"))
  dat
}

classification_list <- lapply(ALL_GENES, classify_gene,
                               hscore_df = hscores)
classification      <- dplyr::bind_rows(classification_list)

write.csv(classification,
          file.path(OUT_DIR, "sample_classification.csv"),
          row.names = FALSE)

log_msg("Classification complete.")
# Print classification summary — one row per sample, one col per gene
class_wide <- classification |>
  dplyr::select(sample_name, gene, classification) |>
  pivot_wider(names_from = gene, values_from = classification)

log_cat(capture.output(print(class_wide, n = 20)), "\n")

# Also print H-Score summary separately
hscore_wide_print <- hscores |>
  dplyr::select(sample_name, gene, H_Score) |>
  pivot_wider(names_from = gene, values_from = H_Score) |>
  mutate(across(where(is.numeric), ~ round(.x, 2)))

log_cat("\nH-Scores per sample:\n")
log_cat(capture.output(print(hscore_wide_print, n = 20)), "\n")

# ── Classification heatmap ────────────────────────────────────────────────────
class_cols <- c(
  "Negative" = "#2166AC",
  "Low"      = "#FEE08B",
  "Positive" = "#D73027",
  "High"     = "#7F0000"
)

p_class_heat <- classification |>
  mutate(gene = factor(gene, levels = ALL_GENES)) |>
  ggplot(aes(x = gene, y = sample_name, fill = classification)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(H_Score, 1)),
            size = 2.5, color = "grey10") +
  scale_fill_manual(values = class_cols, name = "Classification",
                    drop = FALSE) +
  scale_x_discrete(labels = function(x) ifelse(x == "CD274", "PD-L1",
                                        ifelse(x == "PDCD1", "PD-1",
                                        ifelse(x == "ERBB2", "HER2", x)))) +
  labs(
    title    = "Sample Classification per Marker",
    subtitle = "Based on K-means clustering of Spatial H-Score (k chosen by gap statistic)\nH-Score value shown in each cell",
    x        = NULL, y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8),
    axis.text.x   = element_text(angle = 35, hjust = 1, face = "bold"),
    axis.text.y   = element_text(size = 9),
    legend.position = "right"
  )

ggsave(file.path(PLOT_DIR, "sample_classification_heatmap.png"),
       p_class_heat, width = 12, height = 7, dpi = 300, bg = "white")
log_msg("Classification heatmap saved.")


# ── Summary printout ──────────────────────────────────────────────────────────

log_msg(strrep("=", 60))
log_msg("SCRIPT 03 COMPLETE")
log_msg(sprintf("Outputs written to: %s", OUT_DIR))
log_msg(sprintf("  threshold_table.csv        : %d genes", nrow(threshold_table)))
log_msg(sprintf("  hscores_per_sample.csv     : %d rows", nrow(hscores)))
log_msg(sprintf("  sample_classification.csv  : %d rows", nrow(classification)))
log_msg(strrep("=", 60))

close(.log_con)
