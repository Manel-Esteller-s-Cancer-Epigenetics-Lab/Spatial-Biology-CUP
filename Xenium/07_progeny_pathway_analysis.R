# =============================================================================
# 07_progeny_pathway_analysis.R
# -----------------------------------------------------------------------------
# PURPOSE : Pseudobulk PROGENy pathway activity inference on CUP/IC samples.
#           Replicates and extends the analysis from the CUP paper (Visium),
#           applied to Xenium single-cell data with CancerFinder3 annotations.
#
# APPROACH:
#   1. Aggregate raw counts per sample per compartment (pseudobulk)
#      from BPCells on-disk matrices — no Xenium reload needed
#   2. Filter low-expressed genes (≥50 counts across samples, matching paper)
#   3. CPM normalise (matching paper)
#   4. Run PROGENy via progeny::progeny(scale=TRUE) — matching paper exactly
#      14 pathways × 500 top responsive genes
#   5. Classify samples into ONC-CUP vs IMI-CUP via hierarchical clustering
#   6. Produce heatmaps for Tumor-only and Tumor+Stroma analyses
#
# ANALYSES:
#   A. Tumor cells only    — matches paper exactly
#   B. Tumor + Stroma separately — novel comparison showing TME pathway context
#
# INPUT  : cohort_slim.rds      (for Compartment_binary labels per cell)
#          per_sample/SAMPLE/bpcells/counts/  (raw count matrices)
#          sample_sheet.csv
#
# OUTPUT : Analysis/PROGENy/
#   ├── pseudobulk_counts_tumor.csv
#   ├── pseudobulk_counts_stroma.csv
#   ├── progeny_scores_tumor.csv
#   ├── progeny_scores_stroma.csv
#   ├── sample_classification_ONC_IMI.csv
#   ├── plots/
#   │   ├── progeny_heatmap_tumor.png
#   │   ├── progeny_heatmap_stroma.png
#   │   ├── progeny_heatmap_tumor_vs_stroma.png
#   │   ├── progeny_ONC_IMI_classification.png
#   │   └── progeny_pathway_comparison.png
#   └── 07_progeny.log
#
# USAGE  : Rscript 07_progeny_pathway_analysis.R
# =============================================================================

# ── 0. Configuration ──────────────────────────────────────────────────────────

ANALYSIS_ROOT <- "."
COHORT_RDS <- file.path(ANALYSIS_ROOT, "cohort_slim_v2.rds")
SAMPLE_SHEET <- file.path(ANALYSIS_ROOT, "sample_sheet.csv")
OUT_DIR <- file.path(ANALYSIS_ROOT, "PROGENy")
PLOT_DIR <- file.path(ANALYSIS_ROOT, "PLOTs", "07_progeny")
LOG_PATH <- file.path(ANALYSIS_ROOT, "logs", "07_progeny.log")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)


# PROGENy parameters — matching the paper
N_TOP_GENES <- 500 # top responsive genes per pathway
MIN_COUNTS <- 50 # minimum total counts per gene across samples
MIN_CELLS <- 50 # minimum cells contributing to pseudobulk per sample

# Pathways — exactly as defined in the paper
# ONC-CUP: oncogene-driven hyperproliferation and angiogenesis
ONC_PATHWAYS <- c("PI3K", "EGFR", "MAPK", "VEGF")
# IMI-CUP: chronic inflammation and immune-related markers
IMI_PATHWAYS <- c("JAK-STAT", "NFkB", "TNFa")
ALL_PATHWAYS <- c(ONC_PATHWAYS, IMI_PATHWAYS) # 7 pathways total

# Clustering parameters for ONC vs IMI classification
N_CLUSTERS <- 2 # paper found 2 groups


# ── 1. Libraries ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(progeny) # PROGENy — matching paper method
  library(BPCells) # on-disk matrix access
  library(Matrix) # sparse matrix ops
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(pheatmap) # clustered heatmaps
  library(RColorBrewer)
  library(scales)
  library(edgeR) # filterByExpr — matching paper's filtering approach
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

log_msg("Script 07 — PROGENy pseudobulk pathway analysis — started")


# ── 3. Load cohort slim (for compartment labels) ──────────────────────────────

log_msg("Loading cohort_slim.rds for compartment labels...")
cohort <- readRDS(COHORT_RDS) 

# Build per-sample cell ID lookup for Tumor and Stroma
# Key: sample_name → vector of cell_ids per compartment
tumor_cells <- cohort |>
  filter(Compartment_binary == "Tumor") |>
  dplyr::select(sample_name, cell_id)

stroma_cells <- cohort |>
  filter(Compartment_binary == "Stroma") |>
  dplyr::select(sample_name, cell_id)

sample_names <- sort(unique(cohort$sample_name))
log_cat(sprintf(
  "Samples: %d | Tumor cells: %d | Stroma cells: %d\n",
  length(sample_names),
  nrow(tumor_cells),
  nrow(stroma_cells)
))

samples <- read.csv(SAMPLE_SHEET, stringsAsFactors = FALSE) 

# ── 4. Build pseudobulk count matrices from BPCells ───────────────────────────
# For each sample, load BPCells on-disk matrix, subset to compartment cells,
# and sum counts across cells → one count vector per sample per compartment.

log_msg(strrep("─", 60))
log_msg("Building pseudobulk count matrices...")
log_msg(strrep("─", 60))

build_pseudobulk <- function(compartment_label, cell_lookup) {
  # Returns: genes × samples sparse matrix of summed raw counts

  pb_list <- vector("list", nrow(samples))

  for (i in seq_len(nrow(samples))) {
    sname <- samples$sample_name[i]
    bp_path <- file.path(
      ANALYSIS_ROOT,
      "per_sample",
      sname,
      "bpcells",
      "counts"
    )

    log_msg(sprintf(
      "  [%d/%d] %s — %s",
      i,
      nrow(samples),
      sname,
      compartment_label
    ))

    if (!dir.exists(bp_path)) {
      log_msg(sprintf("    BPCells dir not found — SKIPPING"))
      next
    }

    tryCatch(
      {
        # Open BPCells on-disk matrix (genes × cells)
        bp_mat <- BPCells::open_matrix_dir(bp_path)

        # Get cells for this compartment in this sample
        comp_cells <- cell_lookup |>
          filter(sample_name == sname) |>
          pull(cell_id)

        # Intersect with cells actually in the matrix
        cells_in_mat <- intersect(comp_cells, colnames(bp_mat))
        n_cells <- length(cells_in_mat)

        if (n_cells < MIN_CELLS) {
          log_msg(sprintf(
            "    Only %d cells — below MIN_CELLS (%d), SKIPPING",
            n_cells,
            MIN_CELLS
          ))
          next
        }

        log_cat(sprintf(
          "    %d cells in compartment | %d in matrix overlap\n",
          length(comp_cells),
          n_cells
        ))

        # Subset to compartment cells
        bp_subset <- bp_mat[, cells_in_mat, drop = FALSE]

        # Sum across cells → gene-level pseudobulk counts
        # BPCells::rowSums is memory-efficient: reads disk in chunks
        pb_counts <- BPCells::rowSums(bp_subset)

        pb_list[[i]] <- pb_counts
        names(pb_list)[i] <- sname

        rm(bp_mat, bp_subset)
        gc(verbose = FALSE)
      },
      error = function(e) {
        log_msg(sprintf("    ERROR: %s", conditionMessage(e)))
      }
    )
  }

  # Remove NULLs
  pb_list_clean <- Filter(Negate(is.null), pb_list)
  valid_names <- names(pb_list_clean)

  if (length(pb_list_clean) == 0) {
    stop(sprintf(
      "No pseudobulk data built for compartment: %s",
      compartment_label
    ))
  }

  # Align genes across samples (union, fill missing with 0)
  all_genes <- Reduce(union, lapply(pb_list_clean, names))

  pb_mat <- sapply(pb_list_clean, function(v) {
    out <- numeric(length(all_genes))
    names(out) <- all_genes
    matched <- intersect(names(v), all_genes)
    out[matched] <- v[matched]
    out
  })

  # pb_mat is now genes × samples dense matrix
  log_cat(sprintf(
    "  Pseudobulk matrix (%s): %d genes × %d samples\n",
    compartment_label,
    nrow(pb_mat),
    ncol(pb_mat)
  ))

  pb_mat
}

# Build for both compartments
pb_tumor <- build_pseudobulk("Tumor", tumor_cells)
pb_stroma <- build_pseudobulk("Stroma", stroma_cells)

# Save raw pseudobulk counts
write.csv(
  as.data.frame(pb_tumor),
  file.path(OUT_DIR, "pseudobulk_counts_tumor.csv")
)
write.csv(
  as.data.frame(pb_stroma),
  file.path(OUT_DIR, "pseudobulk_counts_stroma.csv")
)
log_msg("Raw pseudobulk count matrices saved.")


# ── 5. Filter and normalise — matching paper's approach ───────────────────────

filter_and_normalise <- function(count_mat, label) {
  log_msg(sprintf("Filtering and normalising (%s)...", label))
  log_cat(sprintf(
    "  Before filtering: %d genes × %d samples\n",
    nrow(count_mat),
    ncol(count_mat)
  ))

  # ── No filterByExpr — paper explicitly comments this out ──────────────────
  # The paper uses all genes present after summing counts
  # Only remove genes with zero counts across ALL samples (truly absent)
  gene_totals <- rowSums(count_mat)
  count_mat <- count_mat[gene_totals > 0, , drop = FALSE]

  log_cat(sprintf(
    "  After removing all-zero genes: %d genes × %d samples\n",
    nrow(count_mat),
    ncol(count_mat)
  ))

  # ── Normalisation — exactly matching paper (line 149-150) ─────────────────
  # Paper: cpm_gex <- t(t(gex + 1) / colSums(gex))
  #        cpm_gex <- log1p(cpm_gex * scale_factor)
  # Note: +1 pseudocount added BEFORE dividing by library size
  scale_factor <- 10000
  lib_sizes <- colSums(count_mat)
  cpm_mat <- t(t(count_mat + 1) / lib_sizes) # pseudocount + fraction
  log_cpm <- log1p(cpm_mat * scale_factor) # log1p(CPM × 10000)

  log_cat(sprintf(
    "  Library sizes range: %.1fM – %.1fM counts\n",
    min(lib_sizes) / 1e6,
    max(lib_sizes) / 1e6
  ))

  list(
    counts = count_mat,
    cpm = cpm_mat,
    log_cpm = log_cpm,
    lib_sizes = lib_sizes
  )
}

norm_tumor <- filter_and_normalise(pb_tumor, "Tumor")
norm_stroma <- filter_and_normalise(pb_stroma, "Stroma")


# ── 6. Retrieve PROGENy network ───────────────────────────────────────────────

log_msg("Running PROGENy — exactly matching paper method...")
log_msg("  progeny::progeny(scale=TRUE, top=500, perm=1)")

# ── 7. Run PROGENy — paper method ─────────────────────────────────────────────
# Paper: progeny(cpm_gex, scale=T, organism="Human", top=500, perm=1)
# progeny() takes genes × samples log-CPM matrix
# scale=TRUE: z-scores scores across samples (critical for between-sample comparison)
# perm=1: no permutation testing
# Returns: samples × pathways matrix

run_progeny_scores <- function(norm_data, label) {
  log_msg(sprintf("  Processing %s...", label))
  mat <- norm_data$log_cpm

  log_cat(sprintf("    Input: %d genes × %d samples\n", nrow(mat), ncol(mat)))

  pb_progeny <- progeny::progeny(
    expr = mat,
    scale = TRUE,
    organism = "Human",
    top = N_TOP_GENES,
    perm = 1
  )

  # Filter to our 7 pathways, transpose to pathways × samples
  pb_progeny <- pb_progeny[,
    colnames(pb_progeny) %in% ALL_PATHWAYS,
    drop = FALSE
  ]
  progeny_mat <- t(as.matrix(pb_progeny))
  progeny_mat <- progeny_mat[, sort(colnames(progeny_mat))]

  log_cat(sprintf(
    "    Output: %d pathways × %d samples\n",
    nrow(progeny_mat),
    ncol(progeny_mat)
  ))

  # Long format for downstream compatibility
  long_result <- as.data.frame(progeny_mat) |>
    tibble::rownames_to_column("source") |>
    tidyr::pivot_longer(-source, names_to = "condition", values_to = "score") |>
    dplyr::mutate(statistic = "progeny")

  list(matrix = progeny_mat, long = long_result)
}

res_tumor <- run_progeny_scores(norm_tumor, "Tumor")
res_stroma <- run_progeny_scores(norm_stroma, "Stroma")

# Assign for downstream compatibility
mlm_tumor <- res_tumor$long
mlm_stroma <- res_stroma$long
mat_tumor <- res_tumor$matrix
mat_stroma <- res_stroma$matrix

# Save scores
write.csv(
  as.data.frame(mat_tumor),
  file.path(OUT_DIR, "progeny_scores_tumor.csv")
)
write.csv(
  as.data.frame(mat_stroma),
  file.path(OUT_DIR, "progeny_scores_stroma.csv")
)
log_msg("PROGENy scores saved.")


# ── 8. Score matrices already built inside run_progeny_scores ────────────────

log_cat(sprintf(
  "\nTumor score matrix: %d pathways × %d samples\n",
  nrow(mat_tumor),
  ncol(mat_tumor)
))
log_cat(sprintf(
  "Stroma score matrix: %d pathways × %d samples\n",
  nrow(mat_stroma),
  ncol(mat_stroma)
))


# ── 9. ONC-CUP vs IMI-CUP classification ─────────────────────────────────────
# Replicate paper's approach: hierarchical clustering of PROGENy scores
# to classify samples into oncogene-driven (ONC) vs inflammatory-immune (IMI)

log_msg(strrep("─", 60))
log_msg("ONC-CUP vs IMI-CUP classification (hierarchical clustering)...")
log_msg(strrep("─", 60))

classify_ONC_IMI <- function(score_mat, n_clusters = 2, label = "") {
  # progeny() with scale=TRUE already z-scored across samples
  # Use raw scores directly — do NOT re-scale
  mat_scaled <- score_mat

  # Hierarchical clustering of samples (columns)
  # Ward.D2 linkage on Euclidean distance — robust for transcriptomic data
  hc_samples <- hclust(
    dist(t(mat_scaled), method = "euclidean"),
    method = "ward.D2"
  )

  # Cut tree into n_clusters
  clusters <- cutree(hc_samples, k = n_clusters)

  # Label clusters by their dominant biological signature:
  # ONC = high MAPK, EGFR, PI3K, VEGF (proliferation/angiogenesis)
  # IMI = high NFkB, JAK-STAT, TNFa, TRAIL (inflammation/immunity)
  # Use paper-defined pathway groupings from config
  onc_pathways <- ONC_PATHWAYS # PI3K, EGFR, MAPK, VEGF
  imi_pathways <- IMI_PATHWAYS # JAK-STAT, NFkB, TNFa

  # For each cluster, compute mean score across ONC and IMI pathways
  cluster_means <- sapply(seq_len(n_clusters), function(k) {
    samps <- names(clusters)[clusters == k]
    if (length(samps) == 0) {
      return(c(onc = 0, imi = 0))
    }
    sub_mat <- score_mat[, samps, drop = FALSE]

    onc_genes <- intersect(onc_pathways, rownames(sub_mat))
    imi_genes <- intersect(imi_pathways, rownames(sub_mat))

    c(
      onc = mean(sub_mat[onc_genes, , drop = FALSE], na.rm = TRUE),
      imi = mean(sub_mat[imi_genes, , drop = FALSE], na.rm = TRUE)
    )
  })

  # Assign cluster labels
  cluster_labels <- ifelse(
    cluster_means["onc", ] > cluster_means["imi", ],
    "ONC-CUP",
    "IMI-CUP"
  )

  sample_class <- data.frame(
    sample_name = names(clusters),
    cluster_id = clusters,
    phenotype = cluster_labels[clusters],
    stringsAsFactors = FALSE
  ) |>
    arrange(sample_name)

  log_cat(sprintf("\n%s Classification:\n", label))
  log_cat(capture.output(print(table(sample_class$phenotype))), "\n")
  log_cat(capture.output(print(sample_class)), "\n")

  list(
    classification = sample_class,
    hclust = hc_samples,
    scaled_mat = mat_scaled,
    clusters = clusters
  )
}

# Classify using Tumor scores (matches paper)
onc_imi_tumor <- classify_ONC_IMI(mat_tumor, N_CLUSTERS, "Tumor-only")
# Also classify using combined Tumor+Stroma average (novel)
mat_combined <- (mat_tumor[, intersect(
  colnames(mat_tumor),
  colnames(mat_stroma)
)] +
  mat_stroma[, intersect(colnames(mat_tumor), colnames(mat_stroma))]) /
  2
onc_imi_combined <- classify_ONC_IMI(
  mat_combined,
  N_CLUSTERS,
  "Tumor+Stroma combined"
)

# Save classifications
write.csv(
  onc_imi_tumor$classification |>
    left_join(
      onc_imi_combined$classification |>
        dplyr::select(sample_name, phenotype_combined = phenotype),
      by = "sample_name"
    ),
  file.path(OUT_DIR, "sample_classification_ONC_IMI.csv"),
  row.names = FALSE
)
log_msg("ONC/IMI classification saved.")


# ── 10. Plots ─────────────────────────────────────────────────────────────────

log_msg(strrep("─", 60))
log_msg("Generating plots...")
log_msg(strrep("─", 60))

# ── Shared colour palette ──────────────────────────────────────────────────────
# Diverging blue-white-red: negative (suppressed) → positive (active)
pathway_cols <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(11, "RdBu"))
)(100)

onc_imi_colors <- c("ONC-CUP" = "#D73027", "IMI-CUP" = "#4393C3")

# Helper: build annotation data frame for pheatmap
make_col_annotation <- function(
  classification_df,
  score_mat,
  extra_cols = NULL
) {
  ann <- classification_df |>
    filter(sample_name %in% colnames(score_mat)) |>
    dplyr::select(sample_name, phenotype) |>
    as.data.frame() |> # ensure plain data.frame, no tibble row names
    tibble::remove_rownames() |>
    tibble::column_to_rownames("sample_name")
  if (!is.null(extra_cols)) {
    ann <- cbind(ann, extra_cols)
  }
  ann
}

# ── Plot A: Tumor-only PROGENy heatmap with ONC/IMI annotation ────────────────
log_msg("Plotting tumor PROGENy heatmap...")

ann_tumor <- make_col_annotation(
  onc_imi_tumor$classification,
  mat_tumor
)

png(
  file.path(PLOT_DIR, "progeny_heatmap_tumor.png"),
  width = 14,
  height = 8,
  units = "in",
  res = 300,
  bg = "white"
)

# Paper uses pheatmap with cluster_cols=T, cluster_rows=T, fontsize=15
# We add ONC/IMI annotation on top
pheatmap::pheatmap(
  mat = mat_tumor, # already scaled (scale=TRUE in progeny())
  color = pathway_cols,
  annotation_col = ann_tumor,
  annotation_colors = list(phenotype = onc_imi_colors),
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  clustering_method = "ward.D2",
  border_color = "white",
  fontsize = 13,
  fontsize_row = 12,
  fontsize_col = 11,
  cellwidth = 40,
  cellheight = 25,
  main = "PROGENy Pathway Activity - Tumour Cells (Pseudobulk)",
  treeheight_row = 30,
  treeheight_col = 30
)

dev.off()
log_msg("Tumor heatmap saved.")

# ── Plot B: Stroma PROGENy heatmap ────────────────────────────────────────────
log_msg("Plotting stroma PROGENy heatmap...")

# No re-scaling — already scaled by progeny()
mat_stroma_scaled <- mat_stroma
hc_stroma <- hclust(
  dist(t(mat_stroma_scaled), method = "euclidean"),
  method = "ward.D2"
)

png(
  file.path(PLOT_DIR, "progeny_heatmap_stroma.png"),
  width = 14,
  height = 8,
  units = "in",
  res = 300,
  bg = "white"
)

pheatmap::pheatmap(
  mat = mat_stroma_scaled[,
    hc_stroma$labels[hc_stroma$order]
  ],
  color = pathway_cols,
  cluster_cols = hc_stroma,
  cluster_rows = TRUE,
  clustering_method = "ward.D2",
  border_color = "white",
  fontsize = 11,
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 40,
  cellheight = 22,
  main = "PROGENy Pathway Activity - Stromal Cells (Pseudobulk)"
)

dev.off()
log_msg("Stroma heatmap saved.")

# ── Plot C: Side-by-side Tumor vs Stroma for shared samples ───────────────────
log_msg("Plotting tumor vs stroma comparison heatmap...")

shared_samples <- intersect(colnames(mat_tumor), colnames(mat_stroma))

# Build combined long data frame for ggplot faceted heatmap
tumor_long <- as.data.frame(mat_tumor[, shared_samples]) |>
  rownames_to_column("pathway") |>
  pivot_longer(-pathway, names_to = "sample_name", values_to = "score") |>
  mutate(compartment = "Tumour")

stroma_long <- as.data.frame(mat_stroma[, shared_samples]) |>
  rownames_to_column("pathway") |>
  pivot_longer(-pathway, names_to = "sample_name", values_to = "score") |>
  mutate(compartment = "Stroma")

combined_long <- bind_rows(tumor_long, stroma_long) |>
  left_join(
    onc_imi_tumor$classification |>
      dplyr::select(sample_name, phenotype),
    by = "sample_name"
  ) |>
  mutate(
    compartment = factor(compartment, levels = c("Tumour", "Stroma")),
    # Order samples by ONC/IMI then by hierarchical clustering order
    sample_name = factor(
      sample_name,
      levels = onc_imi_tumor$hclust$labels[
        onc_imi_tumor$hclust$order
      ]
    ),
    # Order pathways by hierarchical clustering from tumor analysis
    pathway = factor(
      pathway,
      levels = rownames(
        onc_imi_tumor$scaled_mat
      )[
        hclust(
          dist(onc_imi_tumor$scaled_mat),
          method = "ward.D2"
        )$order
      ]
    )
  )

p_comparison <- ggplot(
  combined_long,
  aes(x = sample_name, y = pathway, fill = score)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = round(score, 1)), size = 2.2, color = "grey10") +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    limits = c(-4, 4),
    oob = squish,
    name = "Activity\nscore"
  ) +
  facet_grid(. ~ compartment, scales = "free_x", space = "free_x") +
  labs(
    title = "PROGENy Pathway Activity - Tumour vs Stroma Comparison",
    subtitle = paste0(
      "Pseudobulk MLM scores | Top 500 responsive genes per pathway\n",
      "Positive (red) = pathway active | Negative (blue) = pathway suppressed"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
    axis.text.y = element_text(face = "bold", size = 9),
    strip.text = element_text(face = "bold", size = 11, color = "white"),
    strip.background = element_rect(fill = "grey30", color = NA),
    legend.position = "right",
    panel.spacing = unit(0.5, "lines")
  )

ggsave(
  file.path(PLOT_DIR, "progeny_heatmap_tumor_vs_stroma.png"),
  p_comparison,
  width = 16,
  height = 8,
  dpi = 300,
  bg = "white"
)
log_msg("Tumor vs Stroma comparison heatmap saved.")

# ── Plot D: ONC-CUP vs IMI-CUP classification summary ────────────────────────
log_msg("Plotting ONC/IMI classification summary...")

# Merge tumor scores with classification
class_scores <- mlm_tumor |>
  filter(statistic == "progeny") |>
  dplyr::select(pathway = source, sample_name = condition, score) |>
  left_join(
    onc_imi_tumor$classification |>
      dplyr::select(sample_name, phenotype),
    by = "sample_name"
  ) |>
  mutate(
    sample_name = factor(
      sample_name,
      levels = onc_imi_tumor$hclust$labels[onc_imi_tumor$hclust$order]
    ),
    phenotype = factor(phenotype, levels = c("ONC-CUP", "IMI-CUP"))
  )

# Summary: mean pathway score per phenotype
phenotype_means <- class_scores |>
  group_by(phenotype, pathway) |>
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

p_onc_imi <- ggplot(
  phenotype_means,
  aes(x = reorder(pathway, mean_score), y = mean_score, fill = phenotype)
) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
  scale_fill_manual(values = onc_imi_colors, name = NULL) +
  coord_flip() +
  labs(
    title = "ONC-CUP vs IMI-CUP - Mean Pathway Activity",
    subtitle = paste0(
      "ONC-CUP (n=",
      sum(onc_imi_tumor$classification$phenotype == "ONC-CUP"),
      "): oncogene-driven hyperproliferation & angiogenesis\n",
      "IMI-CUP (n=",
      sum(onc_imi_tumor$classification$phenotype == "IMI-CUP"),
      "): chronic inflammation & immune-related markers"
    ),
    x = "Pathway",
    y = "Mean MLM activity score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
    legend.position = "top"
  )

ggsave(
  file.path(PLOT_DIR, "progeny_ONC_IMI_classification.png"),
  p_onc_imi,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)
log_msg("ONC/IMI classification plot saved.")

# ── Plot E: Per-sample ONC/IMI pathway profile bar chart ──────────────────────
log_msg("Plotting per-sample pathway profiles...")

p_per_sample <- ggplot(
  class_scores,
  aes(x = pathway, y = score, fill = score > 0)
) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_fill_manual(
    values = c("TRUE" = "#D73027", "FALSE" = "#4393C3"),
    guide = "none"
  ) +
  facet_wrap(
    ~ paste0(sample_name, "\n(", phenotype, ")"),
    ncol = 5,
    scales = "free_y"
  ) +
  coord_flip() +
  labs(
    title = "PROGENy Pathway Activity - Per Sample",
    subtitle = "Red = active | Blue = suppressed | Facet label shows ONC/IMI class",
    x = NULL,
    y = "MLM score"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8),
    strip.text = element_text(face = "bold", size = 7),
    axis.text.y = element_text(size = 7)
  )

ggsave(
  file.path(PLOT_DIR, "progeny_per_sample_profiles.png"),
  p_per_sample,
  width = 18,
  height = 12,
  dpi = 300,
  bg = "white"
)
log_msg("Per-sample profiles plot saved.")


# ── 11. Summary statistics ────────────────────────────────────────────────────

log_msg(strrep("=", 60))
log_msg("SUMMARY")
log_msg(strrep("=", 60))

log_cat("\nONC-CUP vs IMI-CUP (Tumor-only, matching paper):\n")
log_cat(
  capture.output(
    print(
      onc_imi_tumor$classification |>
        dplyr::select(sample_name, phenotype) |>
        arrange(phenotype, sample_name)
    )
  ),
  "\n"
)

log_cat("\nPathway scores (Tumor) - top active per sample:\n")
top_pathways <- mlm_tumor |>
  filter(statistic == "progeny") |>
  group_by(condition) |>
  slice_max(order_by = score, n = 3) |>
  dplyr::select(sample = condition, pathway = source, score) |>
  mutate(score = round(score, 3))
log_cat(capture.output(print(top_pathways, n = 50)), "\n")

log_msg(strrep("=", 60))
log_msg("SCRIPT 07 COMPLETE")
log_msg(sprintf("Outputs: %s", OUT_DIR))
log_msg("  pseudobulk_counts_tumor.csv")
log_msg("  pseudobulk_counts_stroma.csv")
log_msg("  progeny_scores_tumor.csv")
log_msg("  progeny_scores_stroma.csv")
log_msg("  sample_classification_ONC_IMI.csv")
log_msg("  plots/progeny_heatmap_tumor.png")
log_msg("  plots/progeny_heatmap_stroma.png")
log_msg("  plots/progeny_heatmap_tumor_vs_stroma.png")
log_msg("  plots/progeny_ONC_IMI_classification.png")
log_msg("  plots/progeny_per_sample_profiles.png")
log_msg(strrep("=", 60))

close(.log_con)
