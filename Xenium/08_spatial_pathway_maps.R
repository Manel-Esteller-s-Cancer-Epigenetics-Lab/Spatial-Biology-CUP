# =============================================================================
# 08_spatial_pathway_maps.R
# -----------------------------------------------------------------------------
# PURPOSE : Generate per-sample spatial maps of PROGENy pathway activity.
#           For each sample: one PNG with 7 spatial panels (one per pathway)
#           showing cells at their X/Y tissue coordinates coloured by
#           weighted pathway activity score.
#
# METHOD  : Weighted gene scoring (Option B)
#           For each cell × pathway:
#             Score = sum(w_g × expr_g) for all pathway target genes g
#           where w_g = PROGENy weight, expr_g = log-normalised expression.
#           This is computed directly from the BPCells on-disk matrices
#           per sample - no full matrix materialisation into RAM.
#
# PATHWAYS:
#   ONC-CUP: PI3K, EGFR, MAPK, VEGF
#   IMI-CUP: JAK-STAT, NFkB, TNFα
#
# OUTPUT : Analysis/PROGENy/spatial_maps/
#   └── <SAMPLE>_spatial_pathway_map.png   (one per sample, 7-panel grid)
#
# USAGE  : Rscript 08_spatial_pathway_maps.R
# =============================================================================

# ── 0. Configuration ──────────────────────────────────────────────────────────

ANALYSIS_ROOT <- "."
COHORT_RDS <- file.path(ANALYSIS_ROOT, "cohort_slim_v2.rds")
SAMPLE_SHEET <- file.path(ANALYSIS_ROOT, "sample_sheet.csv")
PROGENY_DIR <- file.path(ANALYSIS_ROOT, "PROGENy")
OUT_DIR <- file.path(ANALYSIS_ROOT, "PLOTs", "08_spatial_maps")
LOG_PATH <- file.path(ANALYSIS_ROOT, "logs", "08_spatial_maps.log")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Pathways to map - exactly as specified
ONC_PATHWAYS <- c("PI3K", "EGFR", "MAPK", "VEGF")
IMI_PATHWAYS <- c("JAK-STAT", "NFkB", "TNFa")
ALL_PATHWAYS <- c(ONC_PATHWAYS, IMI_PATHWAYS)

# PROGENy top genes per pathway
N_TOP_GENES <- 500

# Plot parameters
POINT_SIZE <- 0.15 # cell dot size - small since cells are dense
POINT_ALPHA <- 0.85
DPI <- 200 # lower than publication for speed; increase for final


# ── 1. Libraries ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(decoupleR)
  library(BPCells)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
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

log_msg("Script 08 - Spatial Pathway Activity Maps - started")


# ── 3. Load PROGENy network ───────────────────────────────────────────────────

log_msg("Fetching PROGENy network (human, top 500 genes)...")

progeny_net <- decoupleR::get_progeny(organism = "human", top = N_TOP_GENES)

# Filter to only the 7 pathways we need
progeny_net <- progeny_net |>
  filter(source %in% ALL_PATHWAYS)

log_cat(sprintf(
  "PROGENy network filtered to %d pathways: %s\n",
  length(unique(progeny_net$source)),
  paste(sort(unique(progeny_net$source)), collapse = ", ")
))
log_cat(sprintf("Total pathway-gene interactions: %d\n", nrow(progeny_net)))

# Build a named list: pathway → named weight vector (gene → weight)
# This enables fast matrix multiplication later
pathway_weights <- lapply(ALL_PATHWAYS, function(pw) {
  sub <- progeny_net |> filter(source == pw)
  setNames(sub$weight, sub$target)
})
names(pathway_weights) <- ALL_PATHWAYS


# ── 4. Load cohort slim for coordinates and compartment labels ────────────────

log_msg("Loading cohort_slim.rds...")
cohort <- readRDS(COHORT_RDS) 
# Keep only spatial coords, compartment, and cell_id
coords_all <- cohort |>
  dplyr::select(cell_id, sample_name, Compartment_binary, x, y) |>
  filter(!is.na(x), !is.na(y))

log_cat(sprintf(
  "Cells with coordinates: %d across %d samples\n",
  nrow(coords_all),
  length(unique(coords_all$sample_name))
))

samples <- read.csv(SAMPLE_SHEET, stringsAsFactors = FALSE) 


# ── 5. Colour palettes ────────────────────────────────────────────────────────

# ONC pathways: warm orange-red gradient
# IMI pathways: cool blue-purple gradient
# Both: grey for zero/near-zero, saturated for high activity

# Shared diverging palette centred at 0
# Negative scores (suppressed): blue | Zero: light grey | Positive: red/orange
make_score_scale <- function(pathway, limits = NULL) {
  if (pathway %in% ONC_PATHWAYS) {
    # Warm: blue (suppressed) → white → orange-red (active)
    scale_color_gradient2(
      low = "#313695",
      mid = "grey90",
      high = "#D73027",
      midpoint = 0,
      limits = limits,
      oob = squish,
      name = "Activity\nscore",
      guide = guide_colorbar(
        barwidth = 0.5,
        barheight = 4,
        title.position = "top"
      )
    )
  } else {
    # Cool: orange (suppressed) → white → blue-purple (active)
    scale_color_gradient2(
      low = "#D73027",
      mid = "grey90",
      high = "#4575B4",
      midpoint = 0,
      limits = limits,
      oob = squish,
      name = "Activity\nscore",
      guide = guide_colorbar(
        barwidth = 0.5,
        barheight = 4,
        title.position = "top"
      )
    )
  }
}

# Pathway panel title colours
PATHWAY_TITLE_COLS <- c(
  "PI3K" = "#B22222",
  "EGFR" = "#CC4400",
  "MAPK" = "#E06000",
  "VEGF" = "#8B0000",
  "JAK-STAT" = "#1A5276",
  "NFkB" = "#154360",
  "TNFa" = "#0B5345"
)

# Compartment outline colours for spatial context
COMP_COLS <- c(
  "Tumor" = "#B22222",
  "Stroma" = "#4393C3",
  "Unassigned" = "grey60"
)


# ── 6. Compute cell-level pathway scores per sample ───────────────────────────
# For each sample:
#   1. Open BPCells on-disk matrix (genes × cells)
#   2. Log-normalise (log1p CPM)
#   3. For each pathway: weighted sum of target gene expression
#   4. Merge with spatial coordinates
#   5. Plot 7-panel spatial map and save

log_msg(strrep("─", 60))
log_msg("Processing samples...")
log_msg(strrep("─", 60))

# Shared spatial theme (dark background)
theme_spatial <- function() {
  theme_void() +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      plot.title = element_text(
        color = "white",
        hjust = 0.5,
        size = 9,
        face = "bold"
      ),
      plot.subtitle = element_text(color = "grey60", hjust = 0.5, size = 6),
      legend.text = element_text(color = "white", size = 5),
      legend.title = element_text(color = "white", size = 6),
      legend.background = element_rect(fill = "black", color = NA),
      plot.margin = margin(3, 3, 3, 3)
    )
}

for (i in seq_len(nrow(samples))) {
  sname <- samples$sample_name[i]
  bp_path <- file.path(ANALYSIS_ROOT, "per_sample", sname, "bpcells", "counts")
  out_png <- file.path(OUT_DIR, paste0(sname, "_spatial_pathway_map.png"))

  log_msg(sprintf("[%d/%d] %s", i, nrow(samples), sname))

  # Checkpoint - skip if already done
  if (file.exists(out_png)) {
    log_msg("  Spatial map exists - skipping.")
    next
  }

  if (!dir.exists(bp_path)) {
    log_msg("  BPCells dir not found - SKIPPING")
    next
  }

  tryCatch(
    {
      # ── 6.1 Get spatial coordinates for this sample ─────────────────────────
      coords_s <- coords_all |> filter(sample_name == sname)
      n_cells <- nrow(coords_s)
      log_cat(sprintf("  Cells with coords: %d\n", n_cells))

      # ── 6.2 Open BPCells matrix and log-normalise ───────────────────────────
      bp_mat <- BPCells::open_matrix_dir(bp_path)

      # Align cells: intersect BPCells columns with cells that have coordinates
      shared_cells <- intersect(colnames(bp_mat), coords_s$cell_id)
      log_cat(sprintf(
        "  Cells in BPCells ∩ coords: %d\n",
        length(shared_cells)
      ))

      if (length(shared_cells) < 100) {
        log_msg("  Too few shared cells - SKIPPING")
        next
      }

      # Subset BPCells to shared cells only
      bp_sub <- bp_mat[, shared_cells, drop = FALSE]

      # Library sizes for CPM normalisation
      lib_sizes <- BPCells::colSums(bp_sub)
      lib_sizes[lib_sizes == 0] <- 1 # guard division by zero

      # ── 6.3 Compute weighted pathway score per cell ─────────────────────────
      # For each pathway:
      #   1. Find intersection of pathway genes with matrix rownames
      #   2. Extract those gene rows from BPCells (streams from disk)
      #   3. Log-CPM normalise the slice
      #   4. Multiply by weights and sum across genes → per-cell score

      pathway_scores <- lapply(ALL_PATHWAYS, function(pw) {
        w_vec <- pathway_weights[[pw]]

        # Intersect pathway genes with genes in this sample's matrix
        genes_avail <- intersect(names(w_vec), rownames(bp_sub))

        if (length(genes_avail) < 10) {
          log_msg(sprintf(
            "    %s: only %d genes available - using NA",
            pw,
            length(genes_avail)
          ))
          return(rep(NA_real_, length(shared_cells)))
        }

        # Extract gene slice from BPCells (only these rows read from disk)
        gene_slice <- bp_sub[genes_avail, , drop = FALSE]

        # Materialise as sparse matrix for arithmetic
        gene_sparse <- as(gene_slice, "dgCMatrix")

        # CPM normalise: divide each column by lib_size/1e6
        # Using sweep on sparse matrix
        gene_cpm <- sweep(
          gene_sparse,
          2,
          lib_sizes[shared_cells] / 1e6,
          FUN = "/"
        )

        # Log1p transform
        gene_logcpm <- log1p(gene_cpm)

        # Weighted sum: w^T × expr_matrix → score per cell
        # weights vector aligned to genes_avail
        w_aligned <- w_vec[genes_avail]

        # Matrix multiplication: (1 × genes) %*% (genes × cells) = (1 × cells)
        # scores <- as.numeric(Matrix::crossprod(
        #   Matrix::sparseVector(
        #     w_aligned,
        #     seq_along(w_aligned),
        #    length(w_aligned)
        #   ),
        #   gene_logcpm
        # ))
        w_aligned <- as.numeric(w_vec[genes_avail])
        scores <- as.numeric(Matrix::t(w_aligned) %*% gene_logcpm)
        # Clean up
        rm(gene_slice, gene_sparse, gene_cpm, gene_logcpm)
        gc(verbose = FALSE)

        scores
      })
      names(pathway_scores) <- ALL_PATHWAYS

      rm(bp_mat, bp_sub)
      gc(verbose = FALSE)

      # ── 6.4 Build plotting data frame ───────────────────────────────────────
      # Merge scores with spatial coordinates
      scores_df <- coords_s |>
        filter(cell_id %in% shared_cells) |>
        arrange(match(cell_id, shared_cells)) # align row order

      for (pw in ALL_PATHWAYS) {
        scores_df[[pw]] <- pathway_scores[[pw]]
      }

      # Compute per-pathway symmetric colour limits
      # Use 99th percentile of absolute scores for robust limits
      pw_limits <- lapply(ALL_PATHWAYS, function(pw) {
        vals <- scores_df[[pw]]
        vals <- vals[!is.na(vals)]
        lim <- quantile(abs(vals), 0.99, na.rm = TRUE)
        c(-lim, lim)
      })
      names(pw_limits) <- ALL_PATHWAYS

      log_cat(sprintf("  Pathway score ranges computed.\n"))

      # ── 6.5 Build 7 spatial panels ──────────────────────────────────────────
      panel_list <- lapply(ALL_PATHWAYS, function(pw) {
        # Sort cells so high-scoring cells plot on top (visible)
        plot_df <- scores_df |>
          filter(!is.na(.data[[pw]])) |>
          arrange(abs(.data[[pw]])) # low scores first, high scores on top

        lims <- pw_limits[[pw]]

        # Group label
        group_label <- if (pw %in% ONC_PATHWAYS) "ONC" else "IMI"

        ggplot(
          plot_df,
          aes(
            x = y,
            y = x, # Xenium: x=row, y=col → flip for display
            color = .data[[pw]]
          )
        ) +
          geom_point(
            size = POINT_SIZE,
            alpha = POINT_ALPHA,
            stroke = 0,
            shape = 16
          ) +
          make_score_scale(pw, limits = lims) +
          coord_fixed() +
          ggtitle(
            label = pw,
            subtitle = group_label
          ) +
          theme_spatial() +
          theme(
            plot.title = element_text(
              color = PATHWAY_TITLE_COLS[[pw]],
              size = 10,
              face = "bold",
              hjust = 0.5
            ),
            plot.subtitle = element_text(
              color = ifelse(pw %in% ONC_PATHWAYS, "#E06000", "#1A5276"),
              size = 7,
              hjust = 0.5
            )
          )
      })
      names(panel_list) <- ALL_PATHWAYS

      # ── 6.6 Assemble 7-panel grid and save ──────────────────────────────────
      # Layout: ONC on top row (4 panels), IMI on bottom row (3 panels)
      p_onc <- wrap_plots(panel_list[ONC_PATHWAYS], nrow = 1) +
        plot_annotation(
          theme = theme(
            plot.background = element_rect(fill = "black", color = NA)
          )
        )

      p_imi <- wrap_plots(panel_list[IMI_PATHWAYS], nrow = 1) +
        plot_annotation(
          theme = theme(
            plot.background = element_rect(fill = "black", color = NA)
          )
        )

      p_combined <- (p_onc / p_imi) +
        plot_annotation(
          title = paste0(sname, " - Spatial Pathway Activity"),
          subtitle = paste0(
            "Top row: ONC-CUP pathways (PI3K | EGFR | MAPK | VEGF)\n",
            "Bottom row: IMI-CUP pathways (JAK-STAT | NFkB | TNFα)\n",
            "Warm red = active | Cool blue = suppressed | Grey = neutral"
          ),
          theme = theme(
            plot.background = element_rect(fill = "black", color = NA),
            plot.title = element_text(
              color = "white",
              hjust = 0.5,
              size = 14,
              face = "bold"
            ),
            plot.subtitle = element_text(
              color = "grey60",
              hjust = 0.5,
              size = 8
            )
          )
        )

      # Dynamic figure dimensions based on tissue aspect ratio
      # Xenium coords are flipped: x→y, y→x, so correct aspect accordingly
      # x_range <- diff(range(scores_df$x, na.rm = TRUE))
      # y_range <- diff(range(scores_df$y, na.rm = TRUE))

      # After coord flip: plot width corresponds to y_range, height to x_range
      # aspect <- x_range / y_range # height-to-width of each panel after flip

      # Cap panel width: reasonable range for a spatial plot
      # panel_w <- max(3, min(4.5, 4 * (y_range / x_range)))
      # fig_w <- panel_w * 4 + 1 # 4 ONC panels + legend
      # fig_h <- panel_w * aspect * 2 + 1.5 # 2 rows + title/subtitle

      fig_w <- 16 # fixed width: 4 panels × ~5 in + margins
      fig_h <- 12 # fixed height: 2 rows × ~5 in + title

      ggsave(
        out_png,
        p_combined,
        width = fig_w,
        height = fig_h,
        dpi = DPI,
        bg = "black",
        limitsize = FALSE
      )

      log_msg(sprintf(
        "  Saved: %s (%.1f × %.1f in)",
        basename(out_png),
        fig_w,
        fig_h
      ))

      rm(scores_df, pathway_scores, panel_list, p_onc, p_imi, p_combined)
      gc(verbose = FALSE)
    },
    error = function(e) {
      log_msg(sprintf("  ERROR: %s", conditionMessage(e)))
    }
  )
} # end sample loop


# ── 7. Summary ────────────────────────────────────────────────────────────────

maps_produced <- list.files(OUT_DIR, pattern = "_spatial_pathway_map\\.png$")

log_msg(strrep("=", 60))
log_msg("SCRIPT 08 COMPLETE")
log_msg(sprintf(
  "Spatial maps produced: %d / %d samples",
  length(maps_produced),
  nrow(samples)
))
log_cat(paste(" ", maps_produced, collapse = "\n"), "\n")
log_msg(sprintf("Output directory: %s", OUT_DIR))
log_msg(strrep("=", 60))

close(.log_con)
