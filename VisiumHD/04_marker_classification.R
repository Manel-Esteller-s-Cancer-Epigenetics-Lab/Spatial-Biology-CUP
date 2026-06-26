# Visium HD 16 µm — Marker Classification

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(mclust)

# Paths
project_dir <- "/home/irsp_carles/phd/projects/CUP/Visium/scripts"
assay_name  <- "Spatial.016um"
image_name  <- "slice1.016um"
out_dir     <- file.path(project_dir, "0.visiumHD/4.marker_classification/results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sample_info <- tibble(
  sample  = c("L02", "L03", "L04_1O"),
  rds     = file.path(project_dir,
              "0.visiumHD/5.UCD/UCD_annotated_objects_16um",
              paste0(c("L02", "L03", "L04_1O"), "_16um_seu_ucd_annotated.rds")),
  ucd_tumor_only_dir = file.path(project_dir,
              paste0("0.visiumHD/5.UCD/UCD_output_", c("L02", "L03", "L04_1O"), "_16um_TumorOnly"))
)

sample_labels <- c(L02 = "CUP_062", L03 = "CUP_138", L04_1O = "CUP_204_1O")
sample_order  <- sample_info$sample

relabel_samples <- function(x) {
  x_chr <- as.character(x)
  out <- unname(sample_labels[x_chr])
  out[is.na(out)] <- x_chr[is.na(out)]
  out
}

# Load annotated Seurat objects
seu_list <- setNames(lapply(sample_info$rds, readRDS), sample_info$sample)
for (nm in names(seu_list)) {
  DefaultAssay(seu_list[[nm]]) <- assay_name
  seu_list[[nm]]$tumor_manual  <- seu_list[[nm]]$ucd_compartment
}

# Marker catalog
# candidates: gene name(s) to try in order (first found is used)
# threshold_method: "Stromal P95" or "GMM"

marker_catalog <- tibble(
  marker = c(
    "KRT7", "KRT20",
    "MET", "ERBB2", "KIT", "FGFR2", "CCNE1", "MYC",
    "CD274", "LAG3", "CTLA4", "TIGIT", "PDCD1",
    "CCL4", "CCL5", "CXCL9", "CXCL10", "CXCL11", "CXCR3",
    "TNFSF10", "TNF",
    "CDX2", "PAX8", "TTF1", "NAPSA", "SYP", "Thyroglobulin",
    "GATA3", "ESR1", "PGR", "CK19", "WT1",
    "PAX2", "AMACR", "CD10", "TP63", "SATB2",
    "ARG1", "HepPar1", "PSMA", "NKX3-1"
  ),
  candidates = list(
    c("KRT7"), c("KRT20"),
    c("MET"), c("ERBB2"), c("KIT"), c("FGFR2"), c("CCNE1"), c("MYC"),
    c("CD274"), c("LAG3"), c("CTLA4"), c("TIGIT"), c("PDCD1"),
    c("CCL4"), c("CCL5"), c("CXCL9"), c("CXCL10"), c("CXCL11"), c("CXCR3"),
    c("TNFSF10"), c("TNF"),
    c("CDX2"), c("PAX8"), c("TTF1", "NKX2-1"), c("NAPSA"), c("SYP"), c("TG"),
    c("GATA3"), c("ESR1"), c("PGR"), c("KRT19"), c("WT1"),
    c("PAX2"), c("AMACR"), c("MME"), c("TP63"), c("SATB2"),
    c("ARG1"), c("CPS1"), c("FOLH1"), c("NKX3-1")
  ),
  marker_group = c(
    rep("Actionable target", 8),
    rep("Immune checkpoint", 5),
    rep("Effector cell traffic", 6),
    rep("Antitumour cytokine", 2),
    rep("Histological marker", 20)
  ),
  threshold_method = c(
    rep("Stromal P95", 8),
    "GMM", "GMM", "GMM", "GMM", "Stromal P95",  # immune checkpoints
    rep("GMM", 6),                               # effector cell traffic
    rep("GMM", 2),                               # antitumour cytokines
    rep("Stromal P95", 20)                       # histological markers
  )
)

# Immune microenvironment markers scored in Tumor + Stroma combined
immune_context_markers  <- c("CD274", "LAG3", "CTLA4", "TIGIT",
                              "CCL4", "CCL5", "CXCL9", "CXCL10", "CXCL11", "CXCR3",
                              "TNFSF10", "TNF")
# PDCD1 evaluated in cancer/tumor compartment only
cancer_only_immune_markers <- c("PDCD1")

# Classification thresholds
min_positive_bins          <- 10
min_total_positive_expr    <- 10
strong_positive_bins       <- 50
strong_total_positive_expr <- 100
binary_evidence_cutoff     <- 0.75
hscore_positive_cutoff     <- 1

# Scoring functions
resolve_marker <- function(counts, candidates) {
  present <- candidates[candidates %in% rownames(counts)]
  if (length(present) == 0) return(NA_character_)
  present[1]
}

get_stroma_threshold <- function(vals_stroma) {
  vals_stroma <- vals_stroma[is.finite(vals_stroma)]
  if (length(vals_stroma) == 0) return(NA_real_)
  thr <- as.numeric(quantile(vals_stroma, 0.95, na.rm = TRUE))
  if (!is.finite(thr)) NA_real_ else thr
}

get_gmm_threshold <- function(values) {
  values <- values[is.finite(values) & values > 0]
  if (length(values) < 50 || length(unique(values)) < 2) return(0)
  log_vals <- log1p(values)
  model <- tryCatch(Mclust(log_vals, G = 2, verbose = FALSE), error = function(e) NULL)
  if (is.null(model) || is.null(model$parameters$mean)) return(0)
  mu     <- model$parameters$mean
  sigmasq <- model$parameters$variance$sigmasq
  sigma  <- if (length(sigmasq) == 1) rep(sqrt(sigmasq), 2) else sqrt(sigmasq)
  if (length(mu) < 2 || length(sigma) < 2) return(0)
  ord    <- order(mu)
  mu     <- mu[ord]; sigma <- sigma[ord]
  x_seq  <- seq(min(log_vals), max(log_vals), length.out = 1000)
  d1     <- dnorm(x_seq, mu[1], sigma[1])
  d2     <- dnorm(x_seq, mu[2], sigma[2])
  idx    <- which(diff(sign(d1 - d2)) != 0)
  threshold <- expm1(if (length(idx) == 0) mean(mu) else x_seq[idx[1]])
  if (!is.finite(threshold)) 0 else threshold
}

score_values_from_threshold <- function(vals, threshold) {
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) {
    return(tibble(n_total_bins = 0L, n_positive_bins = 0L, pct_positive_bins = 0,
                  mean_positive_expr = 0, total_positive_expr = 0, hscore = 0))
  }
  if (is.na(threshold)) threshold <- 0
  positive           <- vals > threshold
  n_positive_bins    <- sum(positive)
  pct_positive_bins  <- mean(positive) * 100
  mean_positive_expr <- if (n_positive_bins > 0) mean(vals[positive], na.rm = TRUE) else 0
  total_positive_expr <- if (n_positive_bins > 0) sum(vals[positive], na.rm = TRUE)  else 0
  tibble(
    n_total_bins        = length(vals),
    n_positive_bins     = n_positive_bins,
    pct_positive_bins   = round(pct_positive_bins,    3),
    mean_positive_expr  = round(mean_positive_expr,   3),
    total_positive_expr = round(total_positive_expr,  3),
    hscore              = round(pct_positive_bins * mean_positive_expr, 6)
  )
}

score_marker_visium <- function(seu, sample_label, marker, candidates,
                                marker_group, threshold_method, assay_name) {
  counts        <- GetAssayData(seu, assay = assay_name, layer = "counts")
  gene_resolved <- resolve_marker(counts, candidates)
  tumor_cells   <- colnames(seu)[seu$tumor_manual == "Tumor"]
  stroma_cells  <- colnames(seu)[seu$tumor_manual == "Stroma"]

  scoring_cells      <- tumor_cells
  threshold_cells    <- tumor_cells
  compartment_label  <- "Tumor"

  if (marker %in% immune_context_markers) {
    scoring_cells   <- c(tumor_cells, stroma_cells)
    threshold_cells <- scoring_cells
    compartment_label <- "Tumor + stroma"
  }

  if (is.na(gene_resolved)) {
    return(tibble(
      sample = sample_label, marker = marker, gene_resolved = NA_character_,
      marker_group = marker_group, compartment = compartment_label,
      threshold_method = threshold_method, threshold_source = NA_character_,
      threshold = NA_real_, n_total_bins = length(scoring_cells),
      n_positive_bins = 0L, pct_positive_bins = 0,
      mean_positive_expr = 0, total_positive_expr = 0, hscore = 0,
      n_tumor_bins = length(tumor_cells), n_stroma_bins = length(stroma_cells),
      detected = FALSE))
  }

  vals_scoring   <- as.numeric(counts[gene_resolved, scoring_cells])
  vals_threshold <- as.numeric(counts[gene_resolved, threshold_cells])
  vals_stroma    <- as.numeric(counts[gene_resolved, stroma_cells])

  if (threshold_method == "GMM") {
    threshold        <- get_gmm_threshold(vals_threshold)
    threshold_source <- if (marker %in% immune_context_markers)
      "Tumour/stroma-bin GMM" else "Tumour-bin GMM"
  } else {
    threshold        <- get_stroma_threshold(vals_stroma)
    threshold_source <- "Sample UCD-derived stromal P95"
  }

  if (marker %in% cancer_only_immune_markers)
    threshold_source <- paste0(threshold_source, " / cancer compartment only")

  if (length(threshold) == 0 || is.na(threshold) || !is.finite(threshold)) {
    threshold        <- 0
    threshold_source <- paste0(threshold_source, " / fallback 0")
  }

  score <- score_values_from_threshold(vals_scoring, threshold)

  tibble(
    sample = sample_label, marker = marker, gene_resolved = gene_resolved,
    marker_group = marker_group, compartment = compartment_label,
    threshold_method = threshold_method, threshold_source = threshold_source,
    threshold         = round(threshold, 6),
    n_total_bins      = score$n_total_bins,
    n_positive_bins   = score$n_positive_bins,
    pct_positive_bins = score$pct_positive_bins,
    mean_positive_expr = score$mean_positive_expr,
    total_positive_expr = score$total_positive_expr,
    hscore            = score$hscore,
    n_tumor_bins      = length(tumor_cells),
    n_stroma_bins     = length(stroma_cells),
    detected          = TRUE
  )
}

# Run marker scoring
message("Scoring markers across all samples...")

marker_scores_df <- bind_rows(lapply(names(seu_list), function(s) {
  message("  ", s)
  seu <- seu_list[[s]]
  bind_rows(lapply(seq_len(nrow(marker_catalog)), function(i) {
    score_marker_visium(
      seu              = seu,
      sample_label     = s,
      marker           = marker_catalog$marker[i],
      candidates       = marker_catalog$candidates[[i]],
      marker_group     = marker_catalog$marker_group[i],
      threshold_method = marker_catalog$threshold_method[i],
      assay_name       = assay_name)
  }))
}))

# Evidence score and binary classification
marker_scores_df <- marker_scores_df %>%
  mutate(
    marker_group_plot = case_when(
      marker %in% c("KRT7", "KRT20") ~ "Epithelial",
      marker_group == "Actionable target" ~ "Actionable target",
      TRUE ~ marker_group
    ),
    log_pct_positive_bins   = log1p(pct_positive_bins),
    log_total_positive_expr = log1p(total_positive_expr),
    log_hscore              = log1p(hscore)
  ) %>%
  group_by(marker_group_plot) %>%
  mutate(
    z_pct_positive_bins   = as.numeric(scale(log_pct_positive_bins)),
    z_total_positive_expr = as.numeric(scale(log_total_positive_expr)),
    z_hscore              = as.numeric(scale(log_hscore)),
    evidence_score        = rowMeans(
      cbind(z_pct_positive_bins, z_total_positive_expr, z_hscore), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    evidence_score  = ifelse(is.finite(evidence_score), evidence_score, NA_real_),
    minimum_signal  = detected & n_positive_bins >= min_positive_bins & total_positive_expr >= min_total_positive_expr,
    strong_signal   = detected & (n_positive_bins >= strong_positive_bins | total_positive_expr >= strong_total_positive_expr),
    evidence_positive = !is.na(evidence_score) & evidence_score >= binary_evidence_cutoff,
    hscore_positive = detected & hscore >= hscore_positive_cutoff,
    status = case_when(
      hscore_positive                          ~ "Positive",
      minimum_signal & (strong_signal | evidence_positive) ~ "Positive",
      TRUE                                     ~ "Negative"
    ),
    sample_label = relabel_samples(sample)
  )

# KRT7/KRT20 immunophenotypic interpretation
krt_status <- marker_scores_df %>%
  filter(marker %in% c("KRT7", "KRT20")) %>%
  select(sample, sample_label, marker, status) %>%
  pivot_wider(names_from = marker, values_from = status)

is_pos <- function(x) x == "Positive"
is_neg <- function(x) x == "Negative"

krt_interpretation <- krt_status %>%
  mutate(
    ck_pattern = case_when(
      is_pos(KRT7) & is_neg(KRT20)  ~ "CK7+/CK20-",
      is_pos(KRT7) & is_pos(KRT20)  ~ "CK7+/CK20+",
      is_neg(KRT7) & is_pos(KRT20)  ~ "CK7-/CK20+",
      is_neg(KRT7) & is_neg(KRT20)  ~ "CK7-/CK20-",
      TRUE                          ~ "Indeterminate"
    ),
    likely_origin = case_when(
      ck_pattern == "CK7+/CK20-" ~ "Lung / Thyroid / Breast / Upper GI-pancreatobiliary / Gynaecological / Renal papillary / Bladder",
      ck_pattern == "CK7+/CK20+" ~ "Bladder / Upper GI-pancreatobiliary / Rectum",
      ck_pattern == "CK7-/CK20+" ~ "Colorectal / Upper GI",
      ck_pattern == "CK7-/CK20-" ~ "Renal / Hepatocellular / Prostate / Gastric / SCLC",
      TRUE                       ~ "Indeterminate"
    )
  )

# Composite immune activity score
score_marker_compartment <- function(seu, sample_label, marker, candidates,
                                     marker_group, assay_name) {
  counts        <- GetAssayData(seu, assay = assay_name, layer = "counts")
  gene_resolved <- resolve_marker(counts, candidates)
  tumor_cells   <- colnames(seu)[seu$tumor_manual == "Tumor"]
  stroma_cells  <- colnames(seu)[seu$tumor_manual == "Stroma"]
  if (is.na(gene_resolved)) {
    return(tibble(sample = sample_label, marker = marker, marker_group = marker_group,
                  compartment = c("Tumor", "Stroma"), hscore = NA_real_, detected = FALSE))
  }
  vals_tumor  <- as.numeric(counts[gene_resolved, tumor_cells])
  vals_stroma <- as.numeric(counts[gene_resolved, stroma_cells])
  threshold   <- get_gmm_threshold(c(vals_tumor, vals_stroma))
  tibble(
    sample      = sample_label,
    marker      = marker,
    marker_group = marker_group,
    compartment = c("Tumor", "Stroma"),
    hscore      = c(score_values_from_threshold(vals_tumor,  threshold)$hscore,
                    score_values_from_threshold(vals_stroma, threshold)$hscore),
    detected    = TRUE
  )
}

immune_compartment_scores <- bind_rows(lapply(names(seu_list), function(s) {
  seu <- seu_list[[s]]
  marker_catalog %>%
    filter(marker_group %in% c("Effector cell traffic", "Antitumour cytokine")) %>%
    split(seq_len(nrow(.))) %>%
    lapply(function(row) {
      score_marker_compartment(seu, s, row$marker, row$candidates[[1]], row$marker_group, assay_name)
    }) %>% bind_rows()
}))

immune_signature_scores <- immune_compartment_scores %>%
  group_by(sample, marker_group, compartment) %>%
  summarise(signature_hscore = mean(hscore, na.rm = TRUE),
            n_markers_detected = sum(detected), .groups = "drop")

traffic_scores <- immune_signature_scores %>%
  filter(marker_group == "Effector cell traffic", compartment == "Tumor") %>%
  select(sample, traffic_score = signature_hscore, traffic_n = n_markers_detected)

cytokine_scores <- immune_signature_scores %>%
  filter(marker_group == "Antitumour cytokine", compartment == "Tumor") %>%
  select(sample, cytokine_score = signature_hscore, cytokine_n = n_markers_detected)

ici_scores <- marker_scores_df %>%
  filter(marker_group == "Immune checkpoint") %>%
  group_by(sample) %>%
  summarise(ici_score = mean(hscore, na.rm = TRUE), ici_n = sum(detected), .groups = "drop")

composite_immune_score <- traffic_scores %>%
  full_join(cytokine_scores, by = "sample") %>%
  full_join(ici_scores, by = "sample") %>%
  mutate(
    composite_score = rowMeans(
      as.matrix(select(., traffic_score, cytokine_score, ici_score)), na.rm = TRUE),
    sample_label = relabel_samples(sample)
  ) %>%
  arrange(desc(composite_score))

# Weighted tumor-origin scores
origin_marker_weights <- tribble(
  ~candidate_primary,          ~marker,       ~weight,
  "Lung",                      "TTF1",         1,
  "Lung",                      "NAPSA",        1,
  "Lung/SCLC",                 "SYP",          1,
  "Thyroid",                   "Thyroglobulin", 1,
  "Thyroid",                   "TTF1",         0.5,
  "Thyroid",                   "PAX8",         0.5,
  "Breast",                    "GATA3",        1,
  "Breast",                    "ESR1",         1,
  "Breast",                    "PGR",          1,
  "Upper GI / Pancreatobiliary","CDX2",        0.5,
  "Upper GI / Pancreatobiliary","CK19",        1,
  "Gynaecological",            "PAX8",         0.5,
  "Gynaecological",            "ESR1",         0.5,
  "Gynaecological",            "PGR",          0.5,
  "Gynaecological",            "WT1",          1,
  "Renal",                     "PAX8",         0.5,
  "Renal",                     "PAX2",         1,
  "Renal",                     "AMACR",        1,
  "Renal",                     "CD10",         1,
  "Bladder",                   "GATA3",        1,
  "Bladder",                   "TP63",         1,
  "Colorectal / Rectum",       "CDX2",         1,
  "Colorectal / Rectum",       "SATB2",        1,
  "Hepatocellular",            "ARG1",         1,
  "Hepatocellular",            "HepPar1",      1,
  "Prostate",                  "PSMA",         1,
  "Prostate",                  "NKX3-1",       1,
  "Gastric",                   "CDX2",         0.5,
  "SCLC",                      "TTF1",         0.5,
  "SCLC",                      "SYP",          1
)

origin_scores <- marker_scores_df %>%
  select(sample, sample_label, marker, hscore) %>%
  inner_join(origin_marker_weights, by = "marker") %>%
  group_by(sample, sample_label, candidate_primary) %>%
  summarise(
    weighted_origin_score = sum(weight * hscore, na.rm = TRUE) / sum(weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(sample, desc(weighted_origin_score))

# UCD tumor-of-origin prediction
tcga_21 <- c("luad", "coad", "brca", "cesc", "chol", "escc", "gbm", "hnsc",
             "lihc", "blca", "lusc", "meso", "ovca", "pdac", "prad", "rcc",
             "skcm", "stad", "thca", "ucec", "uvm")

get_ucd_tumor_means <- function(seu, ucd_cancer_file, sample_label, tcga_21) {
  ucd_cancer      <- read.delim(ucd_cancer_file, row.names = 1, check.names = FALSE)
  tumor_barcodes  <- colnames(seu)[seu$tumor_manual == "Tumor"]
  common_barcodes <- intersect(rownames(ucd_cancer), tumor_barcodes)
  ucd_tumor       <- ucd_cancer[common_barcodes, , drop = FALSE]
  ucd_means       <- colMeans(ucd_tumor, na.rm = TRUE)

  tibble(cancer_type = tcga_21) %>%
    left_join(
      tibble(cancer_type = names(ucd_means), mean_score = as.numeric(ucd_means)) %>%
        filter(cancer_type %in% tcga_21),
      by = "cancer_type"
    ) %>%
    mutate(mean_score = replace_na(mean_score, 0),
           sample = sample_label)
}

ucd_scores_long <- bind_rows(lapply(seq_len(nrow(sample_info)), function(i) {
  s         <- sample_info$sample[i]
  ucd_file  <- file.path(sample_info$ucd_tumor_only_dir[i], "ucdbase_cancer.csv")
  if (!file.exists(ucd_file)) {
    warning("UCD cancer file not found for ", s, ": ", ucd_file)
    return(NULL)
  }
  get_ucd_tumor_means(seu_list[[s]], ucd_file, s, tcga_21)
}))

ucd_top_candidates <- ucd_scores_long %>%
  group_by(sample) %>%
  slice_max(mean_score, n = 5, with_ties = FALSE) %>%
  mutate(rank = row_number(), sample_label = relabel_samples(sample)) %>%
  ungroup()

# Integrated tumor-origin summary
krt_summary <- krt_interpretation %>%
  select(sample, ck_pattern, krt_likely_origin = likely_origin)

weighted_origin_top3 <- origin_scores %>%
  group_by(sample) %>%
  slice_max(weighted_origin_score, n = 3, with_ties = FALSE) %>%
  mutate(rank  = paste0("weighted_top", row_number()),
         label = paste0(candidate_primary, " (", sprintf("%.3f", weighted_origin_score), ")")) %>%
  select(sample, rank, label) %>%
  pivot_wider(names_from = rank, values_from = label)

ucd_top3_summary <- ucd_scores_long %>%
  group_by(sample) %>%
  slice_max(mean_score, n = 3, with_ties = FALSE) %>%
  mutate(rank  = paste0("ucd_top", row_number()),
         label = paste0(toupper(cancer_type), " (", sprintf("%.4f", mean_score), ")")) %>%
  select(sample, rank, label) %>%
  pivot_wider(names_from = rank, values_from = label)

integrated_origin_summary <- krt_summary %>%
  full_join(weighted_origin_top3, by = "sample") %>%
  full_join(ucd_top3_summary,     by = "sample") %>%
  arrange(factor(sample, levels = sample_order)) %>%
  mutate(sample_label = relabel_samples(sample)) %>%
  relocate(sample_label, .after = sample)

# Save outputs
write.csv(marker_scores_df,         file.path(out_dir, "marker_scores.csv"),         row.names = FALSE)
write.csv(krt_interpretation,       file.path(out_dir, "krt_interpretation.csv"),    row.names = FALSE)
write.csv(composite_immune_score,   file.path(out_dir, "composite_immune_score.csv"),row.names = FALSE)
write.csv(origin_scores,            file.path(out_dir, "origin_scores.csv"),          row.names = FALSE)
write.csv(ucd_scores_long,          file.path(out_dir, "ucd_origin_scores.csv"),      row.names = FALSE)
write.csv(integrated_origin_summary,file.path(out_dir, "integrated_origin_summary.csv"), row.names = FALSE)

message("\nMarker classification complete. Results saved to: ", out_dir)
