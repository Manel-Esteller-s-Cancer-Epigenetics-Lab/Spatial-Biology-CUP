# =============================================================================
# 09_immune_infiltration_markers.R
# -----------------------------------------------------------------------------
# PURPOSE : Characterise immune cell infiltration in the tumour compartment
#           of CUP/IC samples using marker-based scoring from the Xenium 5K
#           spatial transcriptomics data.
#
# APPROACH:
#   1. Extract confirmed immune marker genes from BPCells on-disk matrices
#      (same fast approach as Script 01c - no Xenium reload)
#   2. Compute per-cell expression in tumour cells only (Compartment = Tumor)
#   3. Score each immune population per sample using Spatial H-Score
#      (penetrance × intensity), consistent with Scripts 03/05
#   4. Classify samples as enriched / not enriched per population
#   5. Visualise: heatmap, bar charts, immune landscape summary
#
# IMMUNE POPULATIONS (based on panel availability):
#   ✅ B cells            : CD19, MS4A1 (CD20)
#   ✅ cDC1               : CLEC9A, XCR1
#   ⚠️  General DCs        : ITGAX (CD11c) [HLA genes absent from panel]
#   ✅ TAM (general)      : CD68
#   ✅ TAM M2-like        : CD163, MRC1
#   ✅ TAM M1-like        : IL1B, TNF [TNF already in cohort_slim_v2]
#   ⚠️  NK cells           : NCR1, PRF1 [NKG7/GNLY absent]
#   ⚠️  MDSCs              : ARG1 [S100A8/S100A9 absent]
#   ⚠️  Neutrophils (TANs) : CXCR2 [S100A8/S100A9 absent]
#   ✗  Plasma cells       : NO markers in panel (SDC1/Ig genes all absent)
#
# GENES TO EXTRACT (13 new + TNF already in cohort_slim_v2):
#   CD19, MS4A1, CLEC9A, XCR1, ITGAX, CD68, CD163, MRC1,
#   IL1B, NCR1, PRF1, ARG1, CXCR2
#
# INPUT  : cohort_slim_v2.rds (for Compartment_binary, TNF already present)
#          per_sample/SAMPLE/bpcells/counts/ (for new gene extraction)
#
# OUTPUT : Analysis/ImmuneInfiltration/
#   ├── immune_marker_expression.rds    (per-cell expression, tumor only)
#   ├── immune_hscores.csv              (H-Score per gene per sample)
#   ├── immune_population_scores.csv    (aggregate score per population)
#   ├── immune_enrichment_class.csv     (enriched/not enriched per population)
#   └── plots/
#       ├── immune_heatmap_hscores.png
#       ├── immune_heatmap_enrichment.png
#       ├── immune_population_barplot.png
#       └── immune_landscape_summary.png
# =============================================================================

# ── 0. Configuration ──────────────────────────────────────────────────────────

ANALYSIS_ROOT <- "."
COHORT_V2_RDS <- file.path(ANALYSIS_ROOT, "cohort_slim_v2.rds")
SAMPLE_SHEET <- file.path(ANALYSIS_ROOT, "sample_sheet.csv")
OUT_DIR <- file.path(ANALYSIS_ROOT, "ImmuneInfiltration")
PLOT_DIR <- file.path(ANALYSIS_ROOT, "PLOTs", "09_immune_infiltration")
LOG_PATH <- file.path(ANALYSIS_ROOT, "logs", "09_immune_infiltration.log")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLES_EXCLUDE <- c("IC_002T1")

# ── Genes to extract from BPCells ─────────────────────────────────────────────
# TNF is already in cohort_slim_v2 - skip it here
NEW_IMMUNE_GENES <- c(
  # B cells
  "CD19",
  "MS4A1",
  # cDC1
  "CLEC9A",
  "XCR1",
  # General DC
  "ITGAX",
  # TAM general
  "CD68",
  # TAM M2
  "CD163",
  "MRC1",
  # TAM M1 (TNF already present)
  "IL1B",
  # NK cells
  "NCR1",
  "PRF1",
  # MDSC proxy
  "ARG1",
  # Neutrophil proxy
  "CXCR2"
)

# ARG1 is already in cohort_slim_v2 - remove duplicates
# (will be checked at runtime)

# ── Immune population definitions ──────────────────────────────────────────────
# Each population: list of marker genes + weights
# Weight = 1.0 for defining markers, 0.5 for supporting markers
IMMUNE_POPULATIONS <- list(
  "B cells" = list(
    genes = c("CD19", "MS4A1"),
    weights = c(1.0, 1.0),
    color = "#377EB8",
    note = "Both primary B cell markers"
  ),

  "cDC1" = list(
    genes = c("CLEC9A", "XCR1"),
    weights = c(1.0, 1.0),
    color = "#FF7F00",
    note = "Both cDC1-defining markers"
  ),

  "General DCs" = list(
    genes = c("ITGAX"),
    weights = c(1.0),
    color = "#FDBF6F",
    note = "CD11c only; HLA genes absent from panel"
  ),

  "TAM (general)" = list(
    genes = c("CD68"),
    weights = c(1.0),
    color = "#6A3D9A",
    note = "Pan-macrophage marker"
  ),

  "TAM M2-like" = list(
    genes = c("CD163", "MRC1"),
    weights = c(1.0, 1.0),
    color = "#CAB2D6",
    note = "Immunosuppressive macrophage phenotype"
  ),

  "TAM M1-like" = list(
    genes = c("IL1B", "TNF"),
    weights = c(1.0, 1.0),
    color = "#E31A1C",
    note = "Pro-inflammatory macrophage phenotype"
  ),

  "NK cells" = list(
    genes = c("NCR1", "PRF1"),
    weights = c(1.0, 1.0),
    color = "#33A02C",
    note = "NKG7/GNLY absent from panel"
  ),

  "MDSC proxy" = list(
    genes = c("ARG1"),
    weights = c(1.0),
    color = "#B15928",
    note = "ARG1 only; S100A8/S100A9 absent from panel - interpret with caution"
  ),

  "Neutrophils (TAN)" = list(
    genes = c("CXCR2"),
    weights = c(1.0),
    color = "#A65628",
    note = "CXCR2 only; S100A8/S100A9 absent - interpret with caution"
  )
)

# Positivity threshold - consistent with Scripts 03/05
T_CELL <- log1p(1) # = 0.693; NegControl-derived baseline


# ── 1. Libraries ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(BPCells)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

options(Seurat.object.assay.version = "v5")


# ── 2. Logging ────────────────────────────────────────────────────────────────

dir.create(
  file.path(ANALYSIS_ROOT, "logs"),
  recursive = TRUE,
  showWarnings = FALSE
)
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

log_msg("Script 09 - Immune Infiltration Marker Analysis - started")


# ── 3. Load cohort slim v2 ────────────────────────────────────────────────────

log_msg("Loading cohort_slim_v2.rds...")
cohort <- readRDS(COHORT_V2_RDS) |>
  filter(!sample_name %in% SAMPLES_EXCLUDE)

# Identify which immune genes are already in cohort_slim_v2
already_present <- NEW_IMMUNE_GENES[NEW_IMMUNE_GENES %in% colnames(cohort)]
genes_to_extract <- NEW_IMMUNE_GENES[!NEW_IMMUNE_GENES %in% colnames(cohort)]

log_cat(sprintf(
  "Already in cohort_slim_v2: %s\n",
  paste(already_present, collapse = ", ")
))
log_cat(sprintf(
  "Need to extract from BPCells: %s\n",
  paste(genes_to_extract, collapse = ", ")
))

# Also check for TNF and ARG1 (may already be present)
all_immune_genes_needed <- unique(c(NEW_IMMUNE_GENES, "TNF", "ARG1"))

# Tumour cells only - our analysis is restricted to tumour compartment
cohort_tumor <- cohort |>
  filter(Compartment_binary == "Tumor")

sample_names <- sort(unique(cohort$sample_name))
samples <- read.csv(SAMPLE_SHEET, stringsAsFactors = FALSE) |>
  filter(!sample_name %in% SAMPLES_EXCLUDE)

log_cat(sprintf(
  "Tumour cells: %d across %d samples\n",
  nrow(cohort_tumor),
  length(sample_names)
))


# ── 4. Extract new immune genes from BPCells ──────────────────────────────────

if (length(genes_to_extract) > 0) {
  log_msg(strrep("─", 60))
  log_msg(sprintf(
    "Extracting %d new immune genes from BPCells...",
    length(genes_to_extract)
  ))

  new_expr_list <- vector("list", nrow(samples))

  for (i in seq_len(nrow(samples))) {
    sname <- samples$sample_name[i]
    bp_path <- file.path(
      ANALYSIS_ROOT,
      "per_sample",
      sname,
      "bpcells",
      "counts"
    )

    log_msg(sprintf("  [%d/%d] %s", i, nrow(samples), sname))

    if (!dir.exists(bp_path)) {
      log_msg("    BPCells dir not found - SKIPPING")
      next
    }

    tryCatch(
      {
        bp_mat <- BPCells::open_matrix_dir(bp_path)

        # Only extract genes that are in this sample's matrix
        genes_avail <- intersect(genes_to_extract, rownames(bp_mat))
        genes_absent <- setdiff(genes_to_extract, rownames(bp_mat))

        if (length(genes_absent) > 0) {
          log_cat(sprintf(
            "    Absent from matrix: %s\n",
            paste(genes_absent, collapse = ", ")
          ))
        }

        # Compute library sizes from full matrix
        lib_sizes <- BPCells::colSums(bp_mat)
        lib_sizes[lib_sizes == 0] <- 1

        # Extract and normalise
        raw_slice <- bp_mat[genes_avail, , drop = FALSE]
        raw_sparse <- as(raw_slice, "dgCMatrix")

        # Intersect cells with tumor cells for this sample
        tumor_barcodes <- cohort_tumor |>
          filter(sample_name == sname) |>
          pull(cell_id)

        shared <- intersect(colnames(raw_sparse), tumor_barcodes)

        if (length(shared) == 0) {
          log_msg("    No shared tumor cells - SKIPPING")
          next
        }

        # Subset to tumor cells
        raw_tumor <- raw_sparse[, shared, drop = FALSE]
        ls_tumor <- lib_sizes[shared]

        # log1p CPM normalise
        norm_mat <- sweep(raw_tumor, 2, ls_tumor / 1e6, FUN = "/")
        norm_mat@x <- log1p(norm_mat@x)

        # Build data frame with check.names=FALSE for gene names
        expr_df <- as.data.frame(
          setNames(
            lapply(seq_len(nrow(norm_mat)), function(j) {
              as.numeric(norm_mat[j, ])
            }),
            rownames(norm_mat)
          ),
          check.names = FALSE
        )
        expr_df$cell_id <- shared
        expr_df$sample_name <- sname

        # Add NA for absent genes
        for (g in genes_absent) {
          expr_df[[g]] <- NA_real_
        }

        new_expr_list[[i]] <- expr_df

        rm(bp_mat, raw_slice, raw_sparse, raw_tumor, norm_mat, expr_df)
        gc(verbose = FALSE)

        log_cat(sprintf(
          "    Extracted %d tumor cells × %d genes\n",
          length(shared),
          length(genes_avail)
        ))
      },
      error = function(e) {
        log_msg(sprintf("    ERROR: %s", conditionMessage(e)))
      }
    )
  }

  new_expr_all <- dplyr::bind_rows(Filter(Negate(is.null), new_expr_list)) |>
    filter(!sample_name %in% SAMPLES_EXCLUDE)

  log_cat(sprintf(
    "New immune expression table: %d cells × %d cols\n",
    nrow(new_expr_all),
    ncol(new_expr_all)
  ))

  # Save checkpoint
  saveRDS(new_expr_all, file.path(OUT_DIR, "immune_marker_expression.rds"))
  log_msg("Immune marker expression saved.")
} else {
  log_msg(
    "All immune genes already in cohort_slim_v2 - skipping BPCells extraction."
  )
  new_expr_all <- NULL
}


# ── 5. Build unified tumour expression table ──────────────────────────────────

log_msg("Building unified tumour expression table...")

# Start with cohort_slim_v2 tumour cells
# Select cell_id, sample_name, and any immune genes already present
cols_from_slim <- c(
  "cell_id",
  "sample_name",
  intersect(all_immune_genes_needed, colnames(cohort_tumor))
)

tumor_expr <- cohort_tumor |>
  dplyr::select(all_of(cols_from_slim))

# Merge new genes if extracted
if (!is.null(new_expr_all) && nrow(new_expr_all) > 0) {
  new_cols <- c(
    "cell_id",
    "sample_name",
    setdiff(
      colnames(new_expr_all),
      c("cell_id", "sample_name", colnames(tumor_expr))
    )
  )

  if (length(new_cols) > 2) {
    tumor_expr <- tumor_expr |>
      dplyr::left_join(
        new_expr_all |> dplyr::select(all_of(new_cols)),
        by = c("cell_id", "sample_name")
      )
  }
}

# Final list of immune genes available
avail_immune_genes <- intersect(
  all_immune_genes_needed,
  colnames(tumor_expr)
)
avail_immune_genes <- avail_immune_genes[
  avail_immune_genes != "cell_id" &
    avail_immune_genes != "sample_name"
]

log_cat(sprintf(
  "Unified tumor table: %d cells × %d immune genes\n",
  nrow(tumor_expr),
  length(avail_immune_genes)
))
log_cat(sprintf(
  "Available: %s\n",
  paste(sort(avail_immune_genes), collapse = ", ")
))


# ── 6. Compute H-Score per gene per sample ────────────────────────────────────

log_msg("Computing H-Scores per immune marker gene...")

hscore_list <- lapply(avail_immune_genes, function(g) {
  tumor_expr |>
    dplyr::select(sample_name, expr = all_of(g)) |>
    filter(!is.na(expr)) |>
    group_by(sample_name) |>
    summarise(
      n_tumor = n(),
      n_positive = sum(expr > T_CELL),
      pct_pos = 100 * n_positive / n_tumor,
      mean_int = ifelse(n_positive > 0, mean(expr[expr > T_CELL]), 0),
      H_Score = pct_pos * mean_int,
      .groups = "drop"
    ) |>
    mutate(gene = g)
})

hscores <- dplyr::bind_rows(hscore_list)

write.csv(hscores, file.path(OUT_DIR, "immune_hscores.csv"), row.names = FALSE)

log_cat("\nH-Score summary (top 5 per gene):\n")
log_cat(
  capture.output(
    print(
      hscores |>
        group_by(gene) |>
        slice_max(H_Score, n = 3) |>
        dplyr::select(gene, sample_name, pct_pos, H_Score) |>
        mutate(across(where(is.numeric), ~ round(.x, 2))),
      n = 60
    )
  ),
  "\n"
)


# ── 7. Compute population-level aggregate score ───────────────────────────────

log_msg("Computing population aggregate scores...")

pop_scores_list <- lapply(names(IMMUNE_POPULATIONS), function(pop) {
  pop_def <- IMMUNE_POPULATIONS[[pop]]
  genes <- pop_def$genes
  weights <- setNames(pop_def$weights, genes)

  # Get H-Scores for available genes in this population
  pop_h <- hscores |>
    filter(gene %in% genes, gene %in% avail_immune_genes) |>
    dplyr::select(sample_name, gene, H_Score)

  if (nrow(pop_h) == 0) {
    return(data.frame(
      sample_name = sample_names,
      population = pop,
      pop_score = NA_real_,
      n_genes = 0L,
      stringsAsFactors = FALSE
    ))
  }

  # Weighted mean H-Score across available genes
  pop_h |>
    group_by(sample_name) |>
    summarise(
      pop_score = sum(H_Score * weights[gene], na.rm = TRUE) /
        sum(weights[gene[gene %in% names(weights)]], na.rm = TRUE),
      n_genes = n(),
      .groups = "drop"
    ) |>
    mutate(population = pop)
})

pop_scores <- dplyr::bind_rows(pop_scores_list) |>
  mutate(pop_score = round(pop_score, 3))

write.csv(
  pop_scores,
  file.path(OUT_DIR, "immune_population_scores.csv"),
  row.names = FALSE
)


# ── 8. Classify enrichment per population per sample ─────────────────────────
# Use K-means (k=2) within each population - same approach as Scripts 03/05
# Enriched = cluster with higher mean score

log_msg("Classifying immune enrichment via K-means...")

enrich_list <- lapply(names(IMMUNE_POPULATIONS), function(pop) {
  scores_vec <- pop_scores |>
    filter(population == pop, !is.na(pop_score)) |>
    arrange(sample_name)

  if (nrow(scores_vec) < 3) {
    scores_vec$enrichment <- "Indeterminate"
    return(scores_vec)
  }

  x <- matrix(scores_vec$pop_score, ncol = 1)

  # Gap statistic to choose k=2 vs k=1
  best_k <- 2L
  tryCatch(
    {
      gap <- cluster::clusGap(
        x,
        FUN = kmeans,
        K.max = min(3, nrow(x) - 1),
        B = 50,
        verbose = FALSE
      )
      best_k <- max(
        1L,
        min(
          as.integer(cluster::maxSE(
            gap$Tab[, "gap"],
            gap$Tab[, "SE.sim"],
            method = "Tibs2001SEmax"
          )),
          3L
        )
      )
    },
    error = function(e) NULL
  )

  if (best_k == 1) {
    # All samples in same cluster - check if any signal above baseline
    med_score <- median(scores_vec$pop_score)
    scores_vec$enrichment <- ifelse(
      scores_vec$pop_score > med_score,
      "Enriched",
      "Not enriched"
    )
    return(scores_vec)
  }

  set.seed(42)
  km <- kmeans(x, centers = best_k, nstart = 50)

  # Label: highest centroid = Enriched
  center_order <- rank(km$centers[, 1])
  labels <- if (best_k == 2) {
    c("Not enriched", "Enriched")[center_order]
  } else {
    c("Not enriched", "Low", "Enriched")[center_order]
  }

  scores_vec$enrichment <- labels[km$cluster]
  scores_vec
})

enrichment <- dplyr::bind_rows(enrich_list) |>
  dplyr::select(sample_name, population, pop_score, enrichment)

write.csv(
  enrichment,
  file.path(OUT_DIR, "immune_enrichment_class.csv"),
  row.names = FALSE
)

log_cat("\nImmune enrichment summary:\n")
log_cat(
  capture.output(
    print(
      enrichment |>
        dplyr::select(sample_name, population, pop_score, enrichment) |>
        pivot_wider(names_from = population, values_from = enrichment)
    )
  ),
  "\n"
)


# ── 9. Plots ──────────────────────────────────────────────────────────────────

log_msg(strrep("─", 60))
log_msg("Generating plots...")
log_msg(strrep("─", 60))

pop_colors <- sapply(IMMUNE_POPULATIONS, `[[`, "color")

# ── 9A. H-Score heatmap per gene ──────────────────────────────────────────────
log_msg("Plotting gene-level H-Score heatmap...")

# Build matrix: genes × samples
hscore_wide <- hscores |>
  dplyr::select(gene, sample_name, H_Score) |>
  pivot_wider(
    names_from = sample_name,
    values_from = H_Score,
    values_fill = 0
  ) |>
  column_to_rownames("gene") |>
  as.matrix()

# Gene group annotation for row
gene_to_pop <- lapply(names(IMMUNE_POPULATIONS), function(pop) {
  data.frame(
    gene = IMMUNE_POPULATIONS[[pop]]$genes,
    Population = pop,
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows()

gene_ann <- gene_to_pop |>
  filter(gene %in% rownames(hscore_wide)) |>
  distinct(gene, .keep_all = TRUE) |>
  column_to_rownames("gene")

ann_colors <- list(
  Population = pop_colors[unique(gene_ann$Population)]
)

png(
  file.path(PLOT_DIR, "immune_heatmap_hscores.png"),
  width = 14,
  height = 8,
  units = "in",
  res = 300,
  bg = "white"
)

pheatmap::pheatmap(
  mat = hscore_wide[rownames(gene_ann), , drop = FALSE],
  color = colorRampPalette(c(
    "grey97",
    "#FEE08B",
    "#FC8D59",
    "#D73027",
    "#7F0000"
  ))(100),
  annotation_row = gene_ann,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  clustering_method = "ward.D2",
  border_color = "white",
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 9,
  cellwidth = 35,
  cellheight = 20,
  main = "Immune Marker H-Score - Tumour Cells Only",
  gaps_row = cumsum(sapply(IMMUNE_POPULATIONS, function(p) {
    length(intersect(p$genes, rownames(hscore_wide)))
  }))
)

dev.off()
log_msg("H-Score heatmap saved.")

# ── 9B. Enrichment classification heatmap ─────────────────────────────────────
log_msg("Plotting enrichment classification heatmap...")

enrich_wide <- enrichment |>
  dplyr::select(sample_name, population, enrichment) |>
  pivot_wider(names_from = population, values_from = enrichment) |>
  column_to_rownames("sample_name")

ENRICH_COLS <- c(
  "Enriched" = "#D73027",
  "Low" = "#FEE08B",
  "Not enriched" = "#4393C3",
  "Indeterminate" = "grey70"
)

# Convert to numeric for pheatmap
enrich_num <- as.data.frame(lapply(enrich_wide, function(x) {
  recode(x, "Enriched" = 2, "Low" = 1, "Not enriched" = 0, "Indeterminate" = -1)
}))
rownames(enrich_num) <- rownames(enrich_wide)

png(
  file.path(PLOT_DIR, "immune_heatmap_enrichment.png"),
  width = 13,
  height = 7,
  units = "in",
  res = 300,
  bg = "white"
)

pheatmap::pheatmap(
  mat = t(as.matrix(enrich_num)),
  color = c("grey70", "#4393C3", "#FEE08B", "#D73027"),
  breaks = c(-1.5, -0.5, 0.5, 1.5, 2.5),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  clustering_method = "ward.D2",
  border_color = "white",
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 9,
  cellwidth = 35,
  cellheight = 28,
  main = "Immune Cell Enrichment - Tumour Compartment",
  legend_breaks = c(-1, 0, 1, 2),
  legend_labels = c("Indeterminate", "Not enriched", "Low", "Enriched")
)

dev.off()
log_msg("Enrichment heatmap saved.")

# ── 9C. Population score bar chart per sample ──────────────────────────────────
log_msg("Plotting population score bar chart...")

pop_scores_plot <- pop_scores |>
  filter(!is.na(pop_score)) |>
  mutate(
    population = factor(population, levels = names(IMMUNE_POPULATIONS)),
    sample_name = factor(sample_name, levels = sort(unique(sample_name)))
  )

p_bar <- ggplot(
  pop_scores_plot,
  aes(x = sample_name, y = pop_score, fill = population)
) +
  geom_col(position = position_dodge(0.8), width = 0.75, alpha = 0.9) +
  scale_fill_manual(values = pop_colors, name = "Immune population") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Immune Cell Population Scores - Tumour Compartment",
    subtitle = paste0(
      "Weighted Spatial H-Score (penetrance × intensity)\n",
      "Tumour cells only (CancerFinder3 annotation)"
    ),
    x = NULL,
    y = "Population H-Score"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(PLOT_DIR, "immune_population_barplot.png"),
  p_bar,
  width = 16,
  height = 7,
  dpi = 300,
  bg = "white"
)
log_msg("Population bar chart saved.")

# ── 9D. Immune landscape summary - bubble/dot plot ────────────────────────────
log_msg("Plotting immune landscape summary...")

# Merge scores and enrichment for bubble plot
landscape <- pop_scores |>
  left_join(
    enrichment |>
      dplyr::select(sample_name, population, enrichment),
    by = c("sample_name", "population")
  ) |>
  filter(!is.na(pop_score)) |>
  mutate(
    population = factor(population, levels = rev(names(IMMUNE_POPULATIONS))),
    sample_name = factor(sample_name, levels = sort(unique(sample_name))),
    enrichment = factor(
      enrichment,
      levels = c("Not enriched", "Low", "Enriched", "Indeterminate")
    )
  )

p_landscape <- ggplot(
  landscape,
  aes(
    x = sample_name,
    y = population,
    size = pmax(pop_score, 0.5), # min size for visibility
    fill = enrichment,
    color = enrichment
  )
) +
  geom_point(shape = 21, stroke = 0.8, alpha = 0.85) +
  scale_size_area(
    max_size = 14,
    name = "H-Score",
    guide = guide_legend(override.aes = list(fill = "grey50", color = "grey30"))
  ) +
  scale_fill_manual(values = ENRICH_COLS, name = "Enrichment", drop = FALSE) +
  scale_color_manual(
    values = c(
      "Enriched" = "#7F0000",
      "Low" = "#CC8800",
      "Not enriched" = "#1A5276",
      "Indeterminate" = "grey40"
    ),
    name = "Enrichment",
    drop = FALSE
  ) +
  labs(
    title = "Immune Landscape - CUP Cohort (Tumour Compartment)",
    subtitle = paste0(
      "Bubble size = H-Score | Fill = enrichment status\n",
      "Red = enriched | Yellow = low | Blue = not enriched"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    legend.position = "right"
  )

ggsave(
  file.path(PLOT_DIR, "immune_landscape_summary.png"),
  p_landscape,
  width = 14,
  height = 8,
  dpi = 300,
  bg = "white"
)
log_msg("Immune landscape summary saved.")


# ── 10. Final summary ─────────────────────────────────────────────────────────

log_msg(strrep("=", 60))
log_msg("SCRIPT 09 COMPLETE")
log_msg(sprintf("Immune genes analysed: %d", length(avail_immune_genes)))
log_msg(sprintf("Populations scored:    %d", length(IMMUNE_POPULATIONS)))
log_msg("Outputs:")
log_msg("  immune_marker_expression.rds")
log_msg("  immune_hscores.csv")
log_msg("  immune_population_scores.csv")
log_msg("  immune_enrichment_class.csv")
log_msg("  plots/immune_heatmap_hscores.png")
log_msg("  plots/immune_heatmap_enrichment.png")
log_msg("  plots/immune_population_barplot.png")
log_msg("  plots/immune_landscape_summary.png")
log_msg(strrep("=", 60))

close(.log_con)
