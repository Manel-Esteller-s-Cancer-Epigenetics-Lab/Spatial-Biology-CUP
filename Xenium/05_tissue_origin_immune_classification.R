# =============================================================================
# 05_tissue_origin_immune_classification.R
# -----------------------------------------------------------------------------
# PURPOSE : Two complementary analyses on the expanded marker set:
#
# PART A: Tissue of Origin Scoring
#   For each CUP sample, compute a quantitative "origin score" for each
#   candidate primary site based on the Spatial H-Score of site-specific
#   marker genes. Samples are ranked by their top candidate primary site.
#
#   Candidate primary sites and their marker genes:
#     Colorectal         : CDX2, SATB2, KRT20
#     Liver (HCC)        : CPS1, ARG1
#     Lung adenocarcinoma: NKX2.1
#     Thyroid            : TG, NKX2.1
#     Neuroendocrine     : SYP
#     Breast             : GATA3, ESR1, PGR
#     Gynaecological     : PAX8, WT1, PAX2
#     Melanoma           : SOX10
#     Squamous/Basal     : TP63
#     Prostate           : FOLH1, NKX3.1, AMACR
#     Renal              : PAX8, PAX2, MME
#     Urothelial         : GATA3, TP63
#
# PART B: Immune Signature Scoring
#   Effector Cell Traffic  : CXCL9, CXCL10, CXCL11, CXCR3
#   Antitumor Cytokines    : TNFSF10, TNF
#   Immune Checkpoint      : CD274, LAG3, CTLA4, TIGIT, PDCD1 (from Script 03)
#
#   These are analysed in tumour AND stroma separately.
#   A composite "Immune Activity Score" combines all three signatures.
#
# INPUT  : cohort_slim_v2.rds
#          hscores_per_sample.csv (from Script 03, for ICI genes)
#
# OUTPUT : Classification/
#   ├── origin_scores.csv
#   ├── origin_classification.csv
#   ├── immune_signature_scores.csv
#   ├── plots/
#   │   ├── origin_score_heatmap.png
#   │   ├── origin_top_candidate.png
#   │   ├── immune_traffic_heatmap.png
#   │   ├── immune_cytokine_heatmap.png
#   │   └── immune_composite_score.png
#   └── 05_classification.log
# =============================================================================

# ── 0. Configuration ──────────────────────────────────────────────────────────

ANALYSIS_ROOT <- "."
COHORT_V2_RDS <- file.path(ANALYSIS_ROOT, "cohort_slim_v2.rds")
HSCORES_CSV <- file.path(
  ANALYSIS_ROOT,
  "Classification",
  "hscores_per_sample.csv"
)
OUT_DIR <- file.path(ANALYSIS_ROOT, "Classification")
PLOT_DIR <- file.path(ANALYSIS_ROOT, "PLOTs", "05_origin_immune")
LOG_PATH <- file.path(ANALYSIS_ROOT, "logs", "05_classification.log")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)


# ── Gene display names ─────────────────────────────────────────────────────────
GENE_DISPLAY <- c(
  NKX2.1 = "TTF1 (NKX2.1)",
  NKX3.1 = "NKX3.1",
  CPS1 = "HepPar1 (CPS1)",
  FOLH1 = "PSMA (FOLH1)",
  MME = "CD10 (MME)",
  KRT20 = "CK20 (KRT20)",
  ERBB2 = "HER2 (ERBB2)",
  CD274 = "PD-L1 (CD274)",
  PDCD1 = "PD-1 (PDCD1)",
  TNFSF10 = "TRAIL (TNFSF10)",
  SYP = "Synaptophysin (SYP)",
  TG = "Thyroglobulin (TG)"
)
label_gene <- function(g) {
  ifelse(g %in% names(GENE_DISPLAY), GENE_DISPLAY[g], g)
}

# ── Tissue of origin marker map ────────────────────────────────────────────────
# Each site: list of marker genes with weights
# Weight = 1 for primary markers, 0.5 for secondary/shared markers
ORIGIN_MAP <- list(
  Colorectal = list(
    genes = c("CDX2", "SATB2", "KRT20"),
    weights = c(1.0, 1.0, 0.5)
  ),
  Liver_HCC = list(genes = c("CPS1", "ARG1"), weights = c(1.0, 1.0)),
  Lung_Adeno = list(genes = c("NKX2.1"), weights = c(1.0)),
  Thyroid = list(genes = c("TG", "NKX2.1"), weights = c(1.0, 0.5)),
  Neuroendocrine = list(genes = c("SYP"), weights = c(1.0)),
  Breast = list(genes = c("GATA3", "ESR1", "PGR"), weights = c(1.0, 1.0, 1.0)),
  Gynaecological = list(
    genes = c("PAX8", "WT1", "PAX2"),
    weights = c(1.0, 1.0, 0.5)
  ),
  Melanoma = list(genes = c("SOX10"), weights = c(1.0)),
  Squamous = list(genes = c("TP63"), weights = c(1.0)),
  Prostate = list(
    genes = c("FOLH1", "NKX3.1", "AMACR"),
    weights = c(1.0, 1.0, 1.0)
  ),
  Renal = list(genes = c("PAX8", "PAX2", "MME"), weights = c(1.0, 0.5, 0.5)),
  Urothelial = list(genes = c("GATA3", "TP63"), weights = c(1.0, 0.5))
)

# ── Immune signature gene groups ───────────────────────────────────────────────
GENES_TRAFFIC <- c("CXCL9", "CXCL10", "CXCL11", "CXCR3")
GENES_CYTOKINES <- c("TNFSF10", "TNF")
GENES_ICI <- c("CD274", "LAG3", "CTLA4", "TIGIT", "PDCD1")

# ── Threshold (consistent with Script 03) ─────────────────────────────────────
# Immune/secreted genes: use log1p(1) = 0.693 as conservative t_cell
# (same NegControl-derived fallback used in Script 03)
T_CELL_IMMUNE <- log1p(1)


# ── 1. Libraries ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(knitr)
  library(mclust)
})


# ── 2. Logging ────────────────────────────────────────────────────────────────

.log_con <- file(LOG_PATH, open = "wt")
log_msg <- function(...) {
  txt <- paste0(format(Sys.time(), "[%H:%M:%S] "), paste(..., sep = ""))
  message(txt)
  tryCatch(writeLines(txt, con = .log_con), error = function(e) NULL)
}
log_cat <- function(...) {
  txt <- paste(..., sep = "")
  cat(txt)
  tryCatch(
    writeLines(trimws(txt, "right"), con = .log_con),
    error = function(e) NULL
  )
}

log_msg("Script 05 started")


# ── 3. Load data ──────────────────────────────────────────────────────────────

log_msg("Loading cohort_slim_v2.rds...")
cohort <- readRDS(COHORT_V2_RDS) 

cohort_tumor <- cohort |> filter(Compartment_binary == "Tumor")
cohort_stroma <- cohort |> filter(Compartment_binary == "Stroma")

log_cat(sprintf(
  "Cohort: %d samples | %d cells | %d tumor | %d stroma\n",
  length(unique(cohort$sample_name)),
  nrow(cohort),
  nrow(cohort_tumor),
  nrow(cohort_stroma)
))

# Load existing ICI H-Scores from Script 03
hscores_ici <- read.csv(HSCORES_CSV, stringsAsFactors = FALSE) |>
  filter(gene %in% GENES_ICI)


# =============================================================================
# PART A: Tissue of Origin H-Scores and Classification
# =============================================================================

log_msg(strrep("─", 60))
log_msg("PART A: Tissue of Origin Scoring")
log_msg(strrep("─", 60))

# ── A1. Compute H-Score per sample per origin marker gene ─────────────────────
# Use stromal floor threshold for all histological markers (= 0.693)
T_CELL_ORIGIN <- log1p(1)

# All unique origin marker genes
all_origin_genes <- unique(unlist(lapply(ORIGIN_MAP, `[[`, "genes")))
all_origin_genes <- all_origin_genes[all_origin_genes %in% colnames(cohort)]

log_cat(sprintf(
  "Origin marker genes found in data: %d / %d\n",
  length(all_origin_genes),
  length(unique(unlist(lapply(ORIGIN_MAP, `[[`, "genes"))))
))

# H-Score per sample per gene (tumor cells only)
origin_hscore_list <- lapply(all_origin_genes, function(g) {
  cohort_tumor |>
    dplyr::select(sample_name, expr = all_of(g)) |>
    filter(!is.na(expr)) |>
    group_by(sample_name) |>
    summarise(
      n_tumor = n(),
      n_positive = sum(expr > T_CELL_ORIGIN),
      pct_pos = 100 * n_positive / n_tumor,
      mean_int = ifelse(n_positive > 0, mean(expr[expr > T_CELL_ORIGIN]), 0),
      H_Score = pct_pos * mean_int,
      .groups = "drop"
    ) |>
    mutate(gene = g)
})

origin_hscores <- dplyr::bind_rows(origin_hscore_list)

write.csv(
  origin_hscores,
  file.path(OUT_DIR, "origin_gene_hscores.csv"),
  row.names = FALSE
)
log_msg("Origin gene H-Scores computed.")

# ── A2. Compute weighted origin score per sample per site ─────────────────────
sample_names <- sort(unique(cohort$sample_name))
site_names <- names(ORIGIN_MAP)

origin_scores <- expand.grid(
  sample_name = sample_names,
  site = site_names,
  stringsAsFactors = FALSE
)

origin_scores$origin_score <- mapply(
  function(sname, site) {
    genes <- ORIGIN_MAP[[site]]$genes
    weights <- ORIGIN_MAP[[site]]$weights

    # Get H-Scores for genes present in data
    gene_scores <- sapply(seq_along(genes), function(i) {
      g <- genes[i]
      w <- weights[i]
      if (!g %in% all_origin_genes) {
        return(0)
      }
      hs <- origin_hscores |>
        filter(sample_name == sname, gene == g) |>
        pull(H_Score)
      if (length(hs) == 0) {
        return(0)
      }
      hs[1] * w
    })

    # Weighted sum normalised by sum of weights for present genes
    present_weights <- weights[genes %in% all_origin_genes]
    if (sum(present_weights) == 0) {
      return(0)
    }
    sum(gene_scores) / sum(present_weights)
  },
  origin_scores$sample_name,
  origin_scores$site
)

origin_scores <- origin_scores |>
  mutate(origin_score = round(origin_score, 3)) |>
  arrange(sample_name, desc(origin_score))

write.csv(
  origin_scores,
  file.path(OUT_DIR, "origin_scores.csv"),
  row.names = FALSE
)

log_cat("\nWeighted origin scores (top 3 per sample):\n")
top3 <- origin_scores |>
  group_by(sample_name) |>
  slice_max(order_by = origin_score, n = 3) |>
  ungroup()
log_cat(capture.output(print(top3, n = 50)), "\n")

# ── A3. Top candidate primary site per sample ──────────────────────────────────
origin_top <- origin_scores |>
  group_by(sample_name) |>
  slice_max(order_by = origin_score, n = 1) |>
  ungroup() |>
  rename(top_site = site, top_score = origin_score)

# Add second and third candidates
origin_ranked <- origin_scores |>
  group_by(sample_name) |>
  mutate(rank = rank(-origin_score, ties.method = "first")) |>
  filter(rank <= 3) |>
  mutate(label = paste0(site, " (", round(origin_score, 1), ")")) |>
  summarise(
    top1 = label[rank == 1],
    top2 = ifelse(any(rank == 2), label[rank == 2], NA),
    top3 = ifelse(any(rank == 3), label[rank == 3], NA),
    .groups = "drop"
  )

write.csv(
  origin_ranked,
  file.path(OUT_DIR, "origin_classification.csv"),
  row.names = FALSE
)

log_cat("\nTop candidate primary sites:\n")
log_cat(capture.output(print(origin_ranked)), "\n")

# ── A4. Origin score heatmap ──────────────────────────────────────────────────
log_msg("Plotting origin score heatmap...")

# Normalise scores within each site for colour scaling
# (raw scores have different ranges across sites)
origin_scores_norm <- origin_scores |>
  group_by(site) |>
  mutate(score_norm = scales::rescale(origin_score, to = c(0, 1))) |>
  ungroup() |>
  mutate(site = factor(site, levels = site_names))

p_origin_heat <- ggplot(
  origin_scores_norm,
  aes(x = site, y = sample_name, fill = score_norm)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = ifelse(origin_score > 0, round(origin_score, 1), "")),
    size = 2.8,
    color = "grey10",
    fontface = "bold"
  ) +
  scale_fill_gradientn(
    name = "Normalised\norigin score",
    colors = c("grey97", "#FEE08B", "#FC8D59", "#D73027", "#7F0000"),
    limits = c(0, 1)
  ) +
  scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
  labs(
    title = "Tissue of Origin Score",
    subtitle = paste0(
      "Weighted Spatial H-Score per candidate primary site\n",
      "Raw score shown in each cell | Colour = normalised within site"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9, face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

ggsave(
  file.path(PLOT_DIR, "origin_score_heatmap.png"),
  p_origin_heat,
  width = 14,
  height = 7,
  dpi = 300,
  bg = "white"
)
log_msg("Origin score heatmap saved.")

# ── A5. Top candidate bar chart ───────────────────────────────────────────────
origin_top_plot <- origin_scores |>
  group_by(sample_name) |>
  mutate(rank = rank(-origin_score, ties.method = "first")) |>
  filter(rank <= 3) |>
  ungroup() |>
  mutate(
    site = factor(site, levels = site_names),
    rank_label = paste0("#", rank)
  )

# Colour palette for sites
site_colours <- setNames(
  c(
    "#E41A1C",
    "#377EB8",
    "#4DAF4A",
    "#984EA3",
    "#FF7F00",
    "#A65628",
    "#F781BF",
    "#FFFF33",
    "#999999",
    "#66C2A5",
    "#FC8D62",
    "#8DA0CB"
  ),
  site_names
)

p_top_site <- ggplot(
  origin_top_plot,
  aes(
    x = reorder(sample_name, -origin_score * (rank == 1)),
    y = origin_score,
    fill = site,
    alpha = factor(rank)
  )
) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(
    aes(label = gsub("_", " ", site)),
    position = position_dodge(0.8),
    vjust = -0.3,
    size = 2.5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = site_colours, name = "Candidate site") +
  scale_alpha_manual(
    values = c("1" = 1, "2" = 0.6, "3" = 0.35),
    name = "Rank"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Top 3 Candidate Primary Sites per CUP Sample",
    subtitle = "Bar height = weighted origin score | Opacity = rank",
    x = NULL,
    y = "Weighted Origin Score"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(PLOT_DIR, "origin_top_candidate.png"),
  p_top_site,
  width = 14,
  height = 6,
  dpi = 300,
  bg = "white"
)
log_msg("Top candidate plot saved.")


# =============================================================================
# PART B: Immune Signature Scoring
# =============================================================================

log_msg(strrep("─", 60))
log_msg("PART B: Immune Signature Scoring")
log_msg(strrep("─", 60))

# ── B1. H-Score per sample per immune gene (Tumor AND Stroma) ─────────────────
compute_immune_hscore <- function(
  gene_vec,
  compartment_df,
  t_cell = T_CELL_IMMUNE,
  compartment_label
) {
  lapply(gene_vec, function(g) {
    if (!g %in% colnames(compartment_df)) {
      return(data.frame(
        sample_name = character(),
        gene = character(),
        compartment = character(),
        pct_pos = numeric(),
        mean_int = numeric(),
        H_Score = numeric()
      ))
    }
    compartment_df |>
      dplyr::select(sample_name, expr = all_of(g)) |>
      filter(!is.na(expr)) |>
      group_by(sample_name) |>
      summarise(
        n_cells = n(),
        n_positive = sum(expr > t_cell),
        pct_pos = 100 * n_positive / n_cells,
        mean_int = ifelse(n_positive > 0, mean(expr[expr > t_cell]), 0),
        H_Score = pct_pos * mean_int,
        .groups = "drop"
      ) |>
      mutate(gene = g, compartment = compartment_label)
  }) |>
    dplyr::bind_rows()
}

all_immune_genes <- c(GENES_TRAFFIC, GENES_CYTOKINES)

immune_tumor <- compute_immune_hscore(
  all_immune_genes,
  cohort_tumor,
  T_CELL_IMMUNE,
  "Tumor"
)
immune_stroma <- compute_immune_hscore(
  all_immune_genes,
  cohort_stroma,
  T_CELL_IMMUNE,
  "Stroma"
)
immune_hscores <- dplyr::bind_rows(immune_tumor, immune_stroma)

write.csv(
  immune_hscores,
  file.path(OUT_DIR, "immune_signature_scores.csv"),
  row.names = FALSE
)
log_msg("Immune H-Scores computed.")

# ── B2. Effector Cell Traffic heatmap ─────────────────────────────────────────
plot_immune_heatmap <- function(genes, title_str, filename) {
  dat <- immune_hscores |>
    filter(gene %in% genes) |>
    mutate(
      gene = factor(
        sapply(gene, label_gene),
        levels = sapply(genes, label_gene)
      ),
      compartment = factor(compartment, levels = c("Tumor", "Stroma"))
    )

  p <- ggplot(dat, aes(x = gene, y = sample_name, fill = H_Score)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(
      aes(label = round(H_Score, 1)),
      size = 2.5,
      color = "grey10",
      fontface = "bold"
    ) +
    scale_fill_gradientn(
      name = "H-Score",
      colors = c("grey97", "#C6DBEF", "#6BAED6", "#2171B5", "#084594"),
      limits = c(0, NA)
    ) +
    facet_wrap(~compartment, ncol = 2) +
    labs(
      title = title_str,
      subtitle = "H-Score = % cells > t_cell × mean expression of positive cells",
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8),
      axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold", size = 11)
    )

  ggsave(
    file.path(PLOT_DIR, filename),
    p,
    width = 12,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  log_msg(sprintf("Saved: %s", filename))
  p
}

p_traffic <- plot_immune_heatmap(
  GENES_TRAFFIC,
  "Effector Cell Traffic Signature (CXCL9, CXCL10, CXCL11, CXCR3)",
  "immune_traffic_heatmap.png"
)

p_cytokine <- plot_immune_heatmap(
  GENES_CYTOKINES,
  "Antitumor Cytokine Signature (TRAIL/TNFSF10, TNF)",
  "immune_cytokine_heatmap.png"
)

# ── B3. Composite Immune Activity Score ───────────────────────────────────────
# Three components, each normalised 0-1 within the cohort:
#   1. Effector Cell Traffic  (mean H-Score of CXCL9/10/11/CXCR3 in stroma)
#   2. Antitumor Cytokines    (mean H-Score of TNFSF10/TNF in tumor)
#   3. Immune Checkpoint      (mean H-Score of ICI genes in tumor, from Script 03)
log_msg("Computing composite Immune Activity Score...")

# Component 1: Traffic score (stroma - these chemokines attract T cells from stroma)
traffic_score <- immune_stroma |>
  filter(gene %in% GENES_TRAFFIC) |>
  group_by(sample_name) |>
  summarise(traffic_score = mean(H_Score, na.rm = TRUE), .groups = "drop")

# Component 2: Cytokine score (tumor - antitumor cytokines from tumor cells)
cytokine_score <- immune_tumor |>
  filter(gene %in% GENES_CYTOKINES) |>
  group_by(sample_name) |>
  summarise(cytokine_score = mean(H_Score, na.rm = TRUE), .groups = "drop")

# Component 3: ICI score (CD274: CPS tumor+stroma; LAG3/CTLA4/TIGIT/PDCD1: stroma — from Script 03)
ici_score <- hscores_ici |>
  group_by(sample_name) |>
  summarise(ici_score = mean(H_Score, na.rm = TRUE), .groups = "drop")

# Combine and compute composite score
composite <- traffic_score |>
  left_join(cytokine_score, by = "sample_name") |>
  left_join(ici_score, by = "sample_name") |>
  mutate(
    # Normalise each component 0-1
    traffic_norm = scales::rescale(traffic_score, to = c(0, 1)),
    cytokine_norm = scales::rescale(cytokine_score, to = c(0, 1)),
    ici_norm = scales::rescale(ici_score, to = c(0, 1)),
    # Equal-weight composite
    immune_activity_score = (traffic_norm + cytokine_norm + ici_norm) / 3
  ) |>
  arrange(desc(immune_activity_score))

write.csv(
  composite,
  file.path(OUT_DIR, "immune_composite_score.csv"),
  row.names = FALSE
)

log_cat("\nComposite Immune Activity Scores:\n")
log_cat(
  capture.output(
    print(
      composite |>
        dplyr::select(
          sample_name,
          traffic_score,
          cytokine_score,
          ici_score,
          immune_activity_score
        ) |>
        mutate(across(where(is.numeric), ~ round(.x, 3)))
    )
  ),
  "\n"
)

# ── B4. Composite score plot ───────────────────────────────────────────────────
composite_long <- composite |>
  dplyr::select(sample_name, traffic_norm, cytokine_norm, ici_norm) |>
  pivot_longer(
    cols = c(traffic_norm, cytokine_norm, ici_norm),
    names_to = "component",
    values_to = "score"
  ) |>
  mutate(
    component = recode(
      component,
      "traffic_norm" = "Effector Cell Traffic",
      "cytokine_norm" = "Antitumor Cytokines",
      "ici_norm" = "Immune Checkpoints"
    ),
    sample_name = factor(sample_name, levels = composite$sample_name)
  )

p_composite <- ggplot(
  composite_long,
  aes(x = sample_name, y = score, fill = component)
) +
  geom_col(width = 0.75, color = "white", linewidth = 0.25) +
  geom_point(
    data = composite |>
      mutate(sample_name = factor(sample_name, levels = composite$sample_name)),
    aes(x = sample_name, y = immune_activity_score),
    inherit.aes = FALSE,
    shape = 23,
    size = 3.5,
    fill = "grey10",
    color = "white",
    stroke = 0.5
  ) +
  scale_fill_manual(
    values = c(
      "Effector Cell Traffic" = "#2166AC",
      "Antitumor Cytokines" = "#D73027",
      "Immune Checkpoints" = "#4DAF4A"
    ),
    name = "Signature"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    labels = function(x) round(x, 2)
  ) +
  labs(
    title = "Composite Immune Activity Score",
    subtitle = paste0(
      "Stacked bar = normalised component scores (0-1 each)\n",
      "Diamond = composite score (mean of 3 components)"
    ),
    x = NULL,
    y = "Normalised Score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(PLOT_DIR, "immune_composite_score.png"),
  p_composite,
  width = 12,
  height = 6,
  dpi = 300,
  bg = "white"
)
log_msg("Composite immune score plot saved.")


# =============================================================================
# PART C: Combined Classification Heatmap (all 37 markers)
# K-means Negative/Low/Positive/High for new markers, merged with Script 03
# =============================================================================

log_msg(strrep("─", 60))
log_msg("PART C: Combined classification heatmap (all 37 markers)")
log_msg(strrep("─", 60))

# ── C1. Compute H-Scores for all new marker genes (tumor cells only) ──────────
# Use same threshold as Script 03: log1p(1) = 0.693
T_CELL_NEW <- log1p(1)

all_new_genes <- c(
  GENES_TRAFFIC,
  GENES_CYTOKINES, # immune signatures
  all_origin_genes # histological markers
)
# deduplicate
all_new_genes <- unique(all_new_genes)

log_msg(sprintf(
  "Computing H-Scores for %d new genes...",
  length(all_new_genes)
))

new_hscore_list <- lapply(all_new_genes, function(g) {
  if (!g %in% colnames(cohort_tumor)) {
    return(NULL)
  }
  cohort_tumor |>
    dplyr::select(sample_name, expr = all_of(g)) |>
    filter(!is.na(expr)) |>
    group_by(sample_name) |>
    summarise(
      n_tumor = n(),
      n_positive = sum(expr > T_CELL_NEW),
      pct_pos = 100 * n_positive / n_tumor,
      mean_int = ifelse(n_positive > 0, mean(expr[expr > T_CELL_NEW]), 0),
      H_Score = pct_pos * mean_int,
      t_cell = T_CELL_NEW,
      .groups = "drop"
    ) |>
    mutate(gene = g)
})

new_hscores <- dplyr::bind_rows(Filter(Negate(is.null), new_hscore_list))

write.csv(
  new_hscores,
  file.path(OUT_DIR, "new_markers_hscores.csv"),
  row.names = FALSE
)

# ── C2. K-means classification for new markers ────────────────────────────────
classify_1d <- function(gene_name, hscore_df, kmax = 3) {
  dat <- hscore_df |>
    filter(gene == gene_name) |>
    arrange(sample_name)

  if (nrow(dat) < 3) {
    dat$classification <- "Negative"
    dat$k_used <- 1
    return(dat)
  }

  x <- matrix(dat$H_Score, ncol = 1)
  n_samp <- nrow(dat)

  # Gap statistic to choose k
  best_k <- 2L # default fallback
  set.seed(42)
  tryCatch(
    {
      gap <- cluster::clusGap(
        x,
        FUN = kmeans,
        K.max = min(kmax, n_samp - 1),
        B = 50,
        verbose = FALSE
      )
      k_opt <- cluster::maxSE(
        gap$Tab[, "gap"],
        gap$Tab[, "SE.sim"],
        method = "Tibs2001SEmax"
      )
      best_k <- max(2L, min(as.integer(k_opt), kmax))
    },
    error = function(e) {
      # best_k stays at default 2
    }
  )

  set.seed(42)
  km <- kmeans(x, centers = best_k, nstart = 50, iter.max = 100)

  if (best_k == 2) {
    labels <- c("Negative", "Positive")
  } else {
    labels <- c("Negative", "Low", "High")
  }

  cluster_order <- order(km$centers[, 1])
  label_map <- setNames(labels[seq_along(cluster_order)], cluster_order)
  dat$classification <- label_map[as.character(km$cluster)]
  dat$k_used <- best_k
  dat$classification <- factor(
    dat$classification,
    levels = c("Negative", "Low", "Positive", "High")
  )
  dat
}

log_msg("Running K-means classification for new markers...")
new_class_list <- lapply(all_new_genes, classify_1d, hscore_df = new_hscores)
new_classification <- dplyr::bind_rows(new_class_list)

write.csv(
  new_classification,
  file.path(OUT_DIR, "new_markers_classification.csv"),
  row.names = FALSE
)

# ── C3. Load existing Script 03 classification and combine ────────────────────
existing_class <- read.csv(
  file.path(OUT_DIR, "sample_classification.csv"),
  stringsAsFactors = FALSE
) |>
  dplyr::select(sample_name, gene, H_Score, classification, k_used)

# Combine - new markers + existing markers
combined_class <- dplyr::bind_rows(
  existing_class,
  new_classification |>
    dplyr::select(sample_name, gene, H_Score, classification, k_used)
)

# ── C4. Define gene display order for the heatmap ─────────────────────────────
# Group by biological category for interpretability
GENE_ORDER <- c(
  # Actionable targets (from Script 03)
  "MET",
  "ERBB2",
  "KIT",
  "FGFR2",
  "CCNE1",
  "MYC",
  # ICI (from Script 03)
  "CD274",
  "LAG3",
  "CTLA4",
  "TIGIT",
  "PDCD1",
  # Effector Cell Traffic
  "CXCL9",
  "CXCL10",
  "CXCL11",
  "CXCR3",
  # Antitumor Cytokines
  "TNFSF10",
  "TNF",
  # Histological origin markers - grouped by site
  "KRT20", # CRC/general epithelial
  "CDX2",
  "SATB2", # Colorectal
  "CPS1",
  "ARG1", # Liver
  "NKX2.1",
  "TG",
  "SYP", # Lung/Thyroid/NE
  "GATA3",
  "ESR1",
  "PGR", # Breast
  "PAX8",
  "WT1",
  "PAX2", # Gynaecological/Renal
  "SOX10", # Melanoma
  "TP63", # Squamous
  "AMACR",
  "FOLH1",
  "NKX3.1", # Prostate
  "MME" # Renal/shared
)

# Keep only genes present in combined_class
GENE_ORDER <- GENE_ORDER[GENE_ORDER %in% combined_class$gene]

# Display labels
display_labels <- sapply(GENE_ORDER, label_gene)

# Group annotation for faceting
GENE_GROUPS <- data.frame(
  gene = GENE_ORDER,
  group = c(
    rep("Actionable Targets", 6),
    rep("Immune Checkpoints", 5),
    rep("Effector Cell Traffic", 4),
    rep("Antitumor Cytokines", 2),
    rep("Histological Markers", length(GENE_ORDER) - 17)
  ),
  stringsAsFactors = FALSE
)
GENE_GROUPS <- GENE_GROUPS[GENE_GROUPS$gene %in% combined_class$gene, ]

# ── C5. Build combined classification heatmap ─────────────────────────────────
CLASS_COLS <- c(
  "Negative" = "#2166AC",
  "Low" = "#FEE08B",
  "Positive" = "#D73027",
  "High" = "#7F0000"
)

# heat_dat built from combined_dedup below after deduplication
log_msg("Combined 37-marker classification heatmap saved.")


# ── C6. Deduplicate and build combined classification table ───────────────────
# KRT20 and possibly other genes appear in both Script 03 output and
# new_hscores. Deduplicate keeping Script 03 values (they came first in bind_rows)
combined_dedup <- combined_class |>
  filter(gene %in% GENE_ORDER) |>
  dplyr::select(sample_name, gene, H_Score, classification) |>
  group_by(sample_name, gene) |>
  dplyr::slice(1) |>
  ungroup()

# Rebuild heat_dat from deduplicated data
heat_dat <- combined_dedup |>
  mutate(
    gene = factor(gene, levels = GENE_ORDER),
    classification = factor(
      as.character(classification),
      levels = c("Negative", "Low", "Positive", "High")
    ),
    H_Score_label = round(H_Score, 1)
  ) |>
  left_join(GENE_GROUPS, by = "gene")

# Pivot to wide classification table
combined_wide_class <- combined_dedup |>
  dplyr::select(sample_name, gene, classification) |>
  mutate(
    gene = factor(gene, levels = GENE_ORDER),
    classification = as.character(classification)
  ) |>
  pivot_wider(
    names_from = gene,
    values_from = classification,
    values_fn = dplyr::first
  )

# Pivot to wide H-Score table
combined_wide_hscore <- combined_dedup |>
  dplyr::select(sample_name, gene, H_Score) |>
  mutate(gene = factor(gene, levels = GENE_ORDER)) |>
  pivot_wider(
    names_from = gene,
    values_from = H_Score,
    values_fn = function(x) round(x[1], 2)
  )

log_cat("\nCombined classification (all 37 markers):\n")
log_cat(capture.output(print(combined_wide_class, width = 200)), "\n")

write.csv(
  combined_wide_class,
  file.path(OUT_DIR, "combined_classification_37markers.csv"),
  row.names = FALSE
)
write.csv(
  combined_wide_hscore,
  file.path(OUT_DIR, "combined_hscores_37markers.csv"),
  row.names = FALSE
)
log_msg("Combined classification and H-Score tables saved.")

# ── C7. Build and save two separate classification heatmaps ──────────────────
# Heatmap 1: Molecular markers (Actionable + ICI + Immune signatures) - 18 genes
# Heatmap 2: Histological origin markers - 19 genes

# Shared heatmap theme
heatmap_theme <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        face = "bold",
        size = 8
      ),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold", size = 9, color = "white"),
      strip.background = element_rect(fill = "grey30", color = NA),
      legend.position = "right",
      panel.spacing = unit(0.4, "lines")
    )
}

SUBTITLE_NOTE <- paste0(
  "K-means clustering of Spatial H-Score (gap statistic) | ",
  "H-Score in each cell | ",
  "Blue = Negative | Yellow = Low | Red = Positive | Dark red = High"
)

# ── Heatmap 1: Molecular markers ──────────────────────────────────────────────
MOLECULAR_GROUPS <- c(
  "Actionable Targets",
  "Immune Checkpoints",
  "Effector Cell Traffic",
  "Antitumor Cytokines"
)

heat_mol <- heat_dat |>
  filter(group %in% MOLECULAR_GROUPS) |>
  mutate(group = factor(group, levels = MOLECULAR_GROUPS))

display_labels_mol <- setNames(
  sapply(
    levels(heat_mol$gene)[levels(heat_mol$gene) %in% heat_mol$gene],
    label_gene
  ),
  levels(heat_mol$gene)[levels(heat_mol$gene) %in% heat_mol$gene]
)

p_mol <- ggplot(
  heat_mol,
  aes(x = gene, y = sample_name, fill = classification)
) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(
    aes(label = H_Score_label),
    size = 2.5,
    color = "grey10",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = CLASS_COLS,
    name = "Classification",
    drop = FALSE
  ) +
  scale_x_discrete(labels = display_labels) +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  labs(
    title = "Sample Classification: Molecular Markers (18 genes)",
    subtitle = SUBTITLE_NOTE,
    x = NULL,
    y = NULL
  ) +
  heatmap_theme()

ggsave(
  file.path(PLOT_DIR, "classification_heatmap_molecular_markers.png"),
  p_mol,
  width = 14,
  height = 7,
  dpi = 300,
  bg = "white"
)
log_msg("Molecular markers classification heatmap saved.")

# ── Heatmap 2: Histological origin markers ─────────────────────────────────────
heat_hist <- heat_dat |>
  filter(group == "Histological Markers") |>
  mutate(group = factor(group))

p_hist <- ggplot(
  heat_hist,
  aes(x = gene, y = sample_name, fill = classification)
) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(
    aes(label = H_Score_label),
    size = 2.5,
    color = "grey10",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = CLASS_COLS,
    name = "Classification",
    drop = FALSE
  ) +
  scale_x_discrete(labels = display_labels) +
  labs(
    title = "Sample Classification: Histological Markers (19 genes)",
    subtitle = paste0(
      SUBTITLE_NOTE,
      "
",
      "Columns ordered by candidate primary site: ",
      "CRC | Liver | Lung/Thyroid/NE | Breast | Gynaecological | ",
      "Melanoma | Squamous | Prostate | Renal"
    ),
    x = NULL,
    y = NULL
  ) +
  heatmap_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8.5)
  )

ggsave(
  file.path(PLOT_DIR, "classification_heatmap_histological_markers.png"),
  p_hist,
  width = 14,
  height = 7,
  dpi = 300,
  bg = "white"
)
log_msg("Histological markers classification heatmap saved.")


# ── Final summary ──────────────────────────────────────────────────────────────

log_msg(strrep("=", 60))
log_msg("SCRIPT 05 COMPLETE")
log_msg(sprintf("Outputs: %s", OUT_DIR))
log_msg("  origin_gene_hscores.csv")
log_msg("  origin_scores.csv")
log_msg("  origin_classification.csv")
log_msg("  immune_signature_scores.csv")
log_msg("  immune_composite_score.csv")
log_msg("  plots/origin_score_heatmap.png")
log_msg("  plots/origin_top_candidate.png")
log_msg("  plots/immune_traffic_heatmap.png")
log_msg("  plots/immune_cytokine_heatmap.png")
log_msg("  plots/immune_composite_score.png")
log_msg("  plots/classification_heatmap_molecular_markers.png")
log_msg("  plots/classification_heatmap_histological_markers.png")
log_msg("  new_markers_hscores.csv")
log_msg("  new_markers_classification.csv")
log_msg("  combined_classification_37markers.csv")
log_msg(strrep("=", 60))

close(.log_con)
