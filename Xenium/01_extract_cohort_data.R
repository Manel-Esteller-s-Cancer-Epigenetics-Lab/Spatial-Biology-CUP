# =============================================================================
# 01_extract_cohort_data.R
# -----------------------------------------------------------------------------
# PURPOSE : Loop over all 16 Xenium CUP/IC samples sequentially.
#           For each sample:
#             1. Load raw Xenium data via BPCells-backed on-disk matrix
#             2. QC filtering + LogNormalise
#             3. Export raw count matrix → MTX → h5ad (Python)
#             4. Run CancerFinder3 inference (Python, conda env "scf")
#             5. Load CF3 predictions → assign Compartment_binary
#             6. Sketch-based PCA + UMAP (sketch cells only, no ProjectData)
#             7. Save per-sample UMAP PNG for QC inspection
#             8. Extract slim per-cell data.frame (marker expression,
#                compartment, spatial coords, sketch UMAP for sketch cells)
#             9. Discard Seurat object; force gc() before next sample
#
# MEMORY STRATEGY :
#   - BPCells keeps the full count matrix on disk — only accessed slices
#     are pulled into RAM at any time
#   - Raw counts are written to MTX directly from the BPCells layer,
#     never materialised as a dense in-memory matrix
#   - ScaleData operates only on HVGs (N_HVG genes), not all 5000
#   - ProjectData is omitted — sketch UMAP is sufficient for QC
#   - Aggressive rm() + gc() after every major step
#   - SKETCH_MAX caps sketch size regardless of cell count
#
# LOGGING :
#   - All messages mirrored to terminal (stdout/stderr)
#   - Per-sample log: <ANALYSIS_ROOT>/per_sample/<SAMPLE>/<SAMPLE>_processing.log
#   - Master log:     <ANALYSIS_ROOT>/logs/master_run.log
#
# OUTPUT  : <ANALYSIS_ROOT>/cohort_slim.rds
#           <ANALYSIS_ROOT>/per_sample/<SAMPLE>/
#               ├── input/
#               │   ├── matrix.mtx, barcodes.txt, features.txt, metadata.csv
#               │   ├── <SAMPLE>.h5ad
#               │   └── <SAMPLE>_transposed.h5ad
#               ├── bpcells/              on-disk BPCells matrix cache
#               ├── <SAMPLE>_Prediction.csv
#               ├── <SAMPLE>_UMAP.png
#               ├── <SAMPLE>_slim.rds     checkpoint (skip if present)
#               └── <SAMPLE>_processing.log
#
# USAGE   : Rscript 01_extract_cohort_data.R
#           Resume-safe: completed samples are skipped automatically.
# =============================================================================

# ── 0. Global configuration ──────────────────────────────────────────────────

ANALYSIS_ROOT <- "."
SAMPLE_SHEET <- file.path(ANALYSIS_ROOT, "sample_sheet.csv")

PY_MTX2H5AD <- "path/to/CancerFinder3/mtx_to_h5ad.py"  # download from https://github.com/JiangBioLab/CancerFinder
PY_INFER    <- "path/to/CancerFinder3/infer.py"
CF3_MODEL   <- "path/to/CancerFinder3/models/st_pretrain_article.pkl"
CF3_THRESHOLD <- 0.85
CONDA_ENV <- "scf"

QC_MIN_UMI <- 15
QC_MIN_GENES <- 10

SKETCH_FRAC <- 0.10 # 10 % of cells
SKETCH_MIN <- 5000 # lower bound
SKETCH_MAX <- 50000 # hard upper bound — keeps ScaleData RAM-safe

N_HVG <- 2000
N_PCS <- 50
N_NEIGHBORS <- 20
RESOLUTION <- 0.5

GENES_HISTOLOGY <- c("KRT20")
GENES_ACTIONABLE <- c("MET", "ERBB2", "KIT", "FGFR2", "CCNE1", "MYC")
GENES_ICI <- c("CD274", "LAG3", "CTLA4", "TIGIT", "PDCD1")
ALL_MARKER_GENES <- unique(c(GENES_HISTOLOGY, GENES_ACTIONABLE, GENES_ICI))


# ── 1. Libraries ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells) # must load before first Seurat v5 assay is created
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# Activate Seurat v5 assay so BPCells on-disk backend is used automatically
options(Seurat.object.assay.version = "v5")


# ── 2. Helper functions ───────────────────────────────────────────────────────

make_sample_dirs <- function(sample_name) {
  base <- file.path(ANALYSIS_ROOT, "per_sample", sample_name)
  input <- file.path(base, "input")
  bpcells <- file.path(base, "bpcells")
  logs <- file.path(ANALYSIS_ROOT, "logs")
  for (d in c(input, bpcells, logs)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  list(base = base, input = input, bpcells = bpcells, logs = logs)
}

# ── Dual-output logging (terminal + file) ─────────────────────────────────────
.log_con <- NULL

log_open <- function(path) {
  if (!is.null(.log_con)) {
    try(close(.log_con), silent = TRUE)
  }
  con <- file(path, open = "wt")
  assign(".log_con", con, envir = .GlobalEnv)
  log_msg(sprintf("Log opened  : %s", path))
  log_msg(sprintf("Run started : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}

log_close <- function() {
  if (!is.null(.log_con)) {
    log_msg(sprintf(
      "Run ended   : %s",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
    try(close(.log_con), silent = TRUE)
    assign(".log_con", NULL, envir = .GlobalEnv)
  }
}

log_msg <- function(...) {
  txt <- paste0(format(Sys.time(), "[%H:%M:%S] "), paste(..., sep = ""))
  message(txt) # terminal stderr
  if (!is.null(.log_con)) {
    tryCatch(writeLines(txt, con = .log_con), error = function(e) NULL)
  }
}

log_cat <- function(...) {
  txt <- paste(..., sep = "")
  cat(txt) # terminal stdout
  if (!is.null(.log_con)) {
    tryCatch(
      writeLines(trimws(txt, "right"), con = .log_con),
      error = function(e) NULL
    )
  }
}

# ── Shell command with live output captured to log ────────────────────────────
conda_python <- function(env = CONDA_ENV) {
  paste("conda run -n", env, "python -u")
}

run_cmd <- function(cmd, step_label) {
  log_msg("[CMD] ", step_label)
  log_msg("  → ", cmd)
  out <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
  status <- attr(out, "status")
  if (length(out) > 0) {
    log_cat(paste(out, collapse = "\n"), "\n")
  }
  if (!is.null(status) && status != 0) {
    stop(sprintf("Command failed (exit %d): %s", status, cmd))
  }
  invisible(out)
}

# ── RAM usage reporter ────────────────────────────────────────────────────────
report_mem <- function(label) {
  mem_mb <- sum(gc(verbose = FALSE)[, 2]) * 8 / 1024
  log_msg(sprintf("[MEM] %-40s R heap ~%.0f MB", label, mem_mb))
}


# -- 3. Master log ------------------------------------------------------------
# R does not allow split=TRUE with type="message" sink, so the master log is
# handled by writing to a persistent open connection inside log_msg/log_cat.

dir.create(
  file.path(ANALYSIS_ROOT, "logs"),
  recursive = TRUE,
  showWarnings = FALSE
)
master_log_path <- file.path(ANALYSIS_ROOT, "logs", "master_run.log")
.master_con <- file(master_log_path, open = "at")
assign(".master_con", .master_con, envir = .GlobalEnv)

# Redefine log_msg and log_cat to also write to .master_con
log_msg <- function(...) {
  txt <- paste0(format(Sys.time(), "[%H:%M:%S] "), paste(..., sep = ""))
  message(txt)
  if (!is.null(.log_con)) {
    tryCatch(writeLines(txt, con = .log_con), error = function(e) NULL)
  }
  tryCatch(writeLines(txt, con = .master_con), error = function(e) NULL)
}

log_cat <- function(...) {
  txt <- paste(..., sep = "")
  cat(txt)
  if (!is.null(.log_con)) {
    tryCatch(
      writeLines(trimws(txt, "right"), con = .log_con),
      error = function(e) NULL
    )
  }
  tryCatch(
    writeLines(trimws(txt, "right"), con = .master_con),
    error = function(e) NULL
  )
}

writeLines(
  c(
    "",
    strrep("=", 70),
    paste("MASTER RUN START:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    strrep("=", 70)
  ),
  con = .master_con
)
cat(sprintf("Master log: %s\n", master_log_path))


# ── 4. Load sample sheet ──────────────────────────────────────────────────────

if (!file.exists(SAMPLE_SHEET)) {
  stop("Sample sheet not found: ", SAMPLE_SHEET)
}
samples <- read.csv(SAMPLE_SHEET, stringsAsFactors = FALSE)
stopifnot(all(c("sample_name", "data_dir") %in% colnames(samples)))
message(sprintf("Sample sheet loaded: %d samples", nrow(samples)))
print(samples[, c("sample_name", "data_dir")])


# ── 5. Per-sample processing loop ────────────────────────────────────────────

slim_list <- vector("list", nrow(samples))

for (i in seq_len(nrow(samples))) {
  sname <- samples$sample_name[i]
  data_dir <- samples$data_dir[i]
  dirs <- make_sample_dirs(sname)
  slim_rds <- file.path(dirs$base, paste0(sname, "_slim.rds"))

  log_open(file.path(dirs$base, paste0(sname, "_processing.log")))
  log_msg(strrep("=", 70))
  log_msg(sprintf("[%d/%d] %s", i, nrow(samples), sname))
  log_msg(strrep("=", 70))

  # ── Checkpoint ────────────────────────────────────────────────────────────
  if (file.exists(slim_rds)) {
    log_msg("Slim RDS found — loading checkpoint and skipping.")
    slim_list[[i]] <- readRDS(slim_rds)
    log_close()
    next
  }

  # Wrap in tryCatch so one bad sample never kills the loop
  tryCatch(
    {
      # ── 5.1 Validate paths ──────────────────────────────────────────────────
      if (!dir.exists(data_dir)) {
        stop("data_dir not found: ", data_dir)
      }

      # ── 5.2 Load Xenium → BPCells on-disk backend ───────────────────────────
      log_msg("Loading Xenium data...")
      xen_raw <- LoadXenium(data.dir = data_dir, assay = "Xenium")

      log_cat(sprintf(
        "[%s] %d cells × %d genes | median UMI %.0f | median genes %.0f\n",
        sname,
        ncol(xen_raw),
        nrow(xen_raw),
        median(xen_raw$nCount_Xenium),
        median(xen_raw$nFeature_Xenium)
      ))

      # Write counts to disk as BPCells MatrixDir
      bpcells_path <- file.path(dirs$bpcells, "counts")
      counts_bp <- BPCells::write_matrix_dir(
        mat = GetAssayData(xen_raw, assay = "Xenium", layer = "counts"),
        dir = bpcells_path,
        overwrite = TRUE
      )
      # Swap in-memory sparse matrix for the on-disk reference
      xen_raw[["Xenium"]] <- CreateAssay5Object(counts = counts_bp)
      rm(counts_bp)
      gc(verbose = FALSE)
      report_mem("after BPCells write")

      # ── 5.3 QC filter ───────────────────────────────────────────────────────
      n_before <- ncol(xen_raw)
      xen <- subset(
        xen_raw,
        nCount_Xenium >= QC_MIN_UMI &
          nFeature_Xenium >= QC_MIN_GENES
      )
      rm(xen_raw)
      gc(verbose = FALSE)
      log_msg(sprintf(
        "QC: %d → %d cells (removed %d)",
        n_before,
        ncol(xen),
        n_before - ncol(xen)
      ))
      report_mem("after QC")

      # ── 5.4 Normalise ───────────────────────────────────────────────────────
      xen <- NormalizeData(
        xen,
        normalization.method = "LogNormalize",
        assay = "Xenium",
        verbose = FALSE
      )
      report_mem("after normalisation")

      # ── 5.5 Export count matrix (sparse, from BPCells layer) ─────────────────
      log_msg("Exporting count matrix for CancerFinder3...")
      # BPCells objects cannot be coerced with as(...,"CsparseMatrix").
      # BPCells::convert_matrix_type() materialises the on-disk matrix into a
      # standard dgCMatrix (genes x cells) that writeMM can consume directly.
      mat_sparse <- BPCells::convert_matrix_type(
        GetAssayData(xen, assay = "Xenium", layer = "counts"),
        type = "uint32_t" # lossless integer — required by writeMM
      ) |>
        as("dgCMatrix")
      writeMM(mat_sparse, file.path(dirs$input, "matrix.mtx"))
      writeLines(colnames(mat_sparse), file.path(dirs$input, "barcodes.txt"))
      writeLines(rownames(mat_sparse), file.path(dirs$input, "features.txt"))
      write.csv(
        xen@meta.data,
        file.path(dirs$input, "metadata.csv"),
        row.names = TRUE
      )
      log_msg(sprintf(
        "Exported: %d cells × %d genes",
        ncol(mat_sparse),
        nrow(mat_sparse)
      ))
      rm(mat_sparse)
      gc(verbose = FALSE)
      report_mem("after MTX export")

      # ── 5.6 Python: MTX → h5ad ──────────────────────────────────────────────
      h5ad_path <- file.path(dirs$input, paste0(sname, ".h5ad"))
      h5ad_T_path <- file.path(dirs$input, paste0(sname, "_transposed.h5ad"))

      run_cmd(
        paste(
          conda_python(),
          PY_MTX2H5AD,
          "--mtx_dir",
          shQuote(dirs$input),
          "--out_h5ad",
          shQuote(h5ad_path),
          "--out_h5ad_T",
          shQuote(h5ad_T_path)
        ),
        "mtx_to_h5ad.py"
      )

      # ── 5.7 Python: CancerFinder3 inference ─────────────────────────────────
      cf3_csv <- file.path(dirs$base, paste0(sname, "_Prediction.csv"))

      run_cmd(
        paste(
          conda_python(),
          PY_INFER,
          paste0("--ckp=", shQuote(CF3_MODEL)),
          paste0("--matrix=", shQuote(h5ad_T_path)),
          paste0("--out=", shQuote(cf3_csv)),
          "--threshold",
          CF3_THRESHOLD
        ),
        "infer.py (CancerFinder3)"
      )

      # ── 5.8 CF3 predictions → Compartment_binary ────────────────────────────
      if (!file.exists(cf3_csv)) {
        stop("CancerFinder3 output not found: ", cf3_csv)
      }

      cf_preds <- read.csv(cf3_csv, stringsAsFactors = FALSE) |>
        dplyr::mutate(CF_label = ifelse(predict >= 0.5, "Tumor", "Stroma"))

      if (!all(c("sample", "predict") %in% colnames(cf_preds))) {
        stop(
          "CF3 CSV missing expected columns. Found: ",
          paste(colnames(cf_preds), collapse = ", ")
        )
      }

      cells_s <- colnames(xen)
      n_matched <- sum(cells_s %in% cf_preds$sample)
      match_pct <- 100 * n_matched / length(cells_s)
      log_cat(sprintf(
        "CF3 barcode match: %d / %d (%.1f%%)\n",
        n_matched,
        length(cells_s),
        match_pct
      ))
      if (match_pct < 90) {
        log_msg(sprintf(
          "WARNING: low barcode match %.1f%% — check barcode format",
          match_pct
        ))
      }

      cf_score_v <- setNames(cf_preds$predict, cf_preds$sample)
      cf_label_v <- setNames(cf_preds$CF_label, cf_preds$sample)
      xen <- AddMetaData(xen, cf_score_v[cells_s], col.name = "CF_score")
      xen <- AddMetaData(xen, cf_label_v[cells_s], col.name = "CF_label")
      xen$Compartment_binary <- dplyr::case_when(
        xen$CF_label == "Tumor" ~ "Tumor",
        xen$CF_label == "Stroma" ~ "Stroma",
        TRUE ~ "Unassigned"
      )
      log_cat("Compartment distribution:\n")
      log_cat(
        paste(
          capture.output(print(table(xen$Compartment_binary))),
          collapse = "\n"
        ),
        "\n"
      )

      rm(cf_preds, cf_score_v, cf_label_v, cells_s)
      gc(verbose = FALSE)
      report_mem("after CF3 annotation")

      # ── 5.9 Sketch → PCA → UMAP (sketch cells only) ──────────────────────────
      sketch_size <- max(SKETCH_MIN, floor(SKETCH_FRAC * ncol(xen)))
      sketch_size <- min(sketch_size, SKETCH_MAX, ncol(xen))
      log_msg(sprintf(
        "Sketching: %d / %d cells (%.1f%%)",
        sketch_size,
        ncol(xen),
        100 * sketch_size / ncol(xen)
      ))

      DefaultAssay(xen) <- "Xenium"
      # LeverageScore densifies the BPCells matrix internally and fails on large
      # samples. "uniform" random sampling is RAM-safe and sufficient for QC UMAP.
      xen <- SketchData(
        object = xen,
        ncells = sketch_size,
        method = "Uniform",
        sketched.assay = "sketch",
        verbose = FALSE
      )
      report_mem("after SketchData")

      DefaultAssay(xen) <- "sketch"
      xen <- FindVariableFeatures(
        xen,
        selection.method = "vst",
        nfeatures = N_HVG,
        verbose = FALSE
      )

      # Scale ONLY HVGs — critical for RAM: avoids densifying all 5K genes
      xen <- ScaleData(xen, features = VariableFeatures(xen), verbose = FALSE)
      report_mem("after ScaleData (HVGs only)")

      xen <- RunPCA(xen, npcs = N_PCS, verbose = FALSE)
      report_mem("after PCA")

      xen <- FindNeighbors(
        xen,
        dims = 1:N_PCS,
        k.param = N_NEIGHBORS,
        verbose = FALSE
      )
      xen <- FindClusters(xen, resolution = RESOLUTION, verbose = FALSE)

      # return.model = FALSE: no projection needed, saves memory
      xen <- RunUMAP(
        xen,
        dims = 1:N_PCS,
        reduction = "pca",
        reduction.name = "umap.sketch",
        return.model = FALSE,
        verbose = FALSE
      )
      report_mem("after UMAP")
      log_msg(sprintf(
        "Sketch clusters: %d",
        length(unique(xen$seurat_clusters))
      ))

      # ── 5.10 UMAP QC PNG ─────────────────────────────────────────────────────
      log_msg("Saving UMAP QC plot...")
      umap_png <- file.path(dirs$base, paste0(sname, "_UMAP.png"))

      p_clust <- DimPlot(
        xen,
        reduction = "umap.sketch",
        group.by = "seurat_clusters",
        pt.size = 0.4,
        label = TRUE,
        repel = TRUE,
        label.size = 3
      ) +
        ggtitle(paste0(sname, " — clusters")) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
        NoLegend()

      p_comp <- DimPlot(
        xen,
        reduction = "umap.sketch",
        group.by = "Compartment_binary",
        cols = c(
          "Tumor" = "#B22222",
          "Stroma" = "#4393C3",
          "Unassigned" = "grey60"
        ),
        pt.size = 0.4
      ) +
        ggtitle(paste0(sname, " — CancerFinder3")) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))

      p_cf <- ggplot(
        data.frame(
          CF_score = xen$CF_score,
          Compartment = xen$Compartment_binary
        ),
        aes(x = Compartment, y = CF_score, fill = Compartment)
      ) +
        geom_violin(trim = FALSE, alpha = 0.75, linewidth = 0.3) +
        geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white") +
        scale_fill_manual(
          values = c(
            "Tumor" = "#B22222",
            "Stroma" = "#4393C3",
            "Unassigned" = "grey60"
          )
        ) +
        labs(
          title = paste0(sname, " — CF3 score"),
          x = NULL,
          y = "CancerFinder3 score"
        ) +
        theme_classic(base_size = 11) +
        theme(
          legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5)
        )

      p_combined <- (p_clust | p_comp | p_cf) +
        plot_annotation(
          title = paste0("QC — ", sname, "  [sketch cells only]"),
          subtitle = sprintf(
            "Total cells: %d | Sketch: %d | Tumor: %d (%.1f%%) | Stroma: %d (%.1f%%)",
            ncol(xen),
            sketch_size,
            sum(xen$Compartment_binary == "Tumor"),
            100 * mean(xen$Compartment_binary == "Tumor"),
            sum(xen$Compartment_binary == "Stroma"),
            100 * mean(xen$Compartment_binary == "Stroma")
          )
        )

      ggsave(
        umap_png,
        p_combined,
        width = 18,
        height = 6,
        dpi = 150,
        bg = "white"
      )
      log_msg(sprintf("UMAP saved: %s", umap_png))
      rm(p_clust, p_comp, p_cf, p_combined)
      gc(verbose = FALSE)

      # ── 5.11 Extract slim data.frame ─────────────────────────────────────────
      # Pull only the marker genes from the BPCells layer (tiny slice of 5K)
      log_msg("Extracting slim data frame...")

      DefaultAssay(xen) <- "Xenium"
      genes_present <- ALL_MARKER_GENES[ALL_MARKER_GENES %in% rownames(xen)]
      genes_absent <- ALL_MARKER_GENES[!ALL_MARKER_GENES %in% rownames(xen)]
      if (length(genes_absent) > 0) {
        log_msg(sprintf(
          "Genes absent (NA): %s",
          paste(genes_absent, collapse = ", ")
        ))
      }

      expr_slice <- GetAssayData(xen, assay = "Xenium", layer = "data")[
        genes_present,
        ,
        drop = FALSE
      ]
      expr_df <- as.data.frame(t(as.matrix(expr_slice)))
      expr_df$cell_id <- rownames(expr_df)
      rm(expr_slice)
      gc(verbose = FALSE)

      for (g in genes_absent) {
        expr_df[[g]] <- NA_real_
      }

      meta_df <- data.frame(
        cell_id = rownames(xen@meta.data),
        sample_name = sname,
        Compartment_binary = xen$Compartment_binary,
        CF_score = xen$CF_score,
        stringsAsFactors = FALSE
      )

      coords <- GetTissueCoordinates(xen, image = "fov", which = "centroids") |>
        dplyr::rename(cell_id = cell) |>
        dplyr::mutate(cell_id = as.character(cell_id))

      # Sketch UMAP: only sketch cells have coordinates; non-sketch cells → NA
      umap_embed <- as.data.frame(Embeddings(xen, reduction = "umap.sketch"))
      colnames(umap_embed) <- c("UMAP_1", "UMAP_2")
      umap_embed$cell_id <- rownames(umap_embed)

      slim <- meta_df |>
        dplyr::left_join(expr_df, by = "cell_id") |>
        dplyr::left_join(coords[, c("cell_id", "x", "y")], by = "cell_id") |>
        dplyr::left_join(umap_embed, by = "cell_id")

      canonical_cols <- c(
        "cell_id",
        "sample_name",
        "Compartment_binary",
        "CF_score",
        "x",
        "y",
        "UMAP_1",
        "UMAP_2",
        GENES_HISTOLOGY,
        GENES_ACTIONABLE,
        GENES_ICI
      )
      for (col in canonical_cols) {
        if (!col %in% colnames(slim)) slim[[col]] <- NA
      }
      slim <- slim[, canonical_cols]

      log_cat(sprintf(
        "Slim table: %d cells × %d columns\n",
        nrow(slim),
        ncol(slim)
      ))
      log_cat(sprintf(
        "Cells with UMAP coords (sketch): %d / %d\n",
        sum(!is.na(slim$UMAP_1)),
        nrow(slim)
      ))

      saveRDS(slim, slim_rds)
      log_msg(sprintf("Slim RDS saved: %s", slim_rds))
      slim_list[[i]] <- slim

      # ── 5.12 Free memory ─────────────────────────────────────────────────────
      rm(xen, expr_df, meta_df, coords, umap_embed, slim)
      gc(verbose = FALSE)
      report_mem(sprintf("end of %s — memory freed", sname))
      log_msg("DONE.")
    },
    error = function(e) {
      log_msg(sprintf("ERROR: %s", conditionMessage(e)))
      log_msg("Sample skipped — continuing to next sample.")
      suppressWarnings(
        for (obj in c(
          "xen",
          "xen_raw",
          "mat_sparse",
          "counts_bp",
          "cf_preds",
          "expr_df",
          "meta_df",
          "slim"
        )) {
          if (exists(obj)) rm(list = obj)
        }
      )
      gc(verbose = FALSE)
    }
  )

  log_close()
} # ── end sample loop ──────────────────────────────────────────────────────────


# ── 6. Combine slim tables ────────────────────────────────────────────────────

message("\n", strrep("=", 70))
message("Combining slim tables into cohort_slim.rds...")

slim_list_clean <- Filter(Negate(is.null), slim_list)
if (length(slim_list_clean) == 0) {
  stop("No slim tables produced.")
}

cohort_slim <- dplyr::bind_rows(slim_list_clean)

cat(sprintf(
  "\nCohort summary:\n  Samples : %d\n  Cells   : %d\n  Columns : %d\n",
  length(unique(cohort_slim$sample_name)),
  nrow(cohort_slim),
  ncol(cohort_slim)
))
cat("\nCells per sample:\n")
print(table(cohort_slim$sample_name))
cat("\nCompartment distribution:\n")
print(table(cohort_slim$sample_name, cohort_slim$Compartment_binary))


# ── 7. Save ───────────────────────────────────────────────────────────────────

cohort_rds <- file.path(ANALYSIS_ROOT, "cohort_slim.rds")
saveRDS(cohort_slim, cohort_rds)
message(sprintf("Saved: %s", cohort_rds))
message("Script 01 complete. Run 02_cohort_visualisation.qmd next.")

close(.master_con) # close master log connection
