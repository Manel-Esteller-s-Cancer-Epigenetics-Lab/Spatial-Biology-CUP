# Load necessary libraries
library(tidyverse)
library(Seurat)
library(patchwork)

projectPath <- "02_ST_CUPS"

source(file.path(projectPath, "misc", "project_utilities.R"))

add_UCD_predictions <- function(args) {
  if (length(args) != 2) {
    stop("2 arguments required: sampleSubtype and sampleId.")
  }
  
  sampleSubtype <- args[1]
  sampleId <- args[2]
  baseDir <- "01_analysis_by_sample"
  step <- "04_sfc_pred"
  
  dirsAndLog <- setupDirectoriesAndLog(projectPath, baseDir, sampleSubtype, sampleId, step)
  logFile <- dirsAndLog$logFile
  dataPath <- file.path(projectPath, "processed_data", baseDir, "03.3_add_UCD_pred", sampleSubtype, sampleId, "data.RDS")
  savePlots <- dirsAndLog$savePlots

  # Load the Seurat object
  if (!file.exists(dataPath)) {
    logMessage("Data file does not exist", logFile)
  }
  data <- readRDS(dataPath)
  
  thresholds <- c("0.50", "0.55", "0.60", "0.65", "0.70", "0.75", "0.80", "0.85", "0.90", "0.95", "1.00")
  
  for (threshold in thresholds) {
    predPath <- file.path(projectPath, "results", baseDir, "07.1_run_scf", sampleSubtype, sampleId, paste0("out_", threshold, ".csv"))
    if (!file.exists(predPath)) {
      next
    }
    
                                                                       
    predData <- read.csv(predPath)
    colnames(predData)[2] <- paste0("scf_", threshold)
    data <- data %>% AddMetaData(metadata = predData)

    logMessage(paste("scf predictions added:", threshold), logFile)

    pt.size.factor <- get_pt_size_factor(sampleId)
    plot <- SpatialDimPlot(data, group.by = paste0("scf_", threshold), stroke = 0, pt.size.factor = pt.size.factor, image.alpha = F) + scale_fill_manual(values = c("0" = "lightblue", "1" = "red"))
    ggsave(plot, file = here::here(savePlots, paste0("scf_threshold_inference_",threshold, ".pdf"))

  }
  
  # Save the updated Seurat object
  saveRDS(data, file.path(dirsAndLog$saveProcessedData, "data.RDS"))
}

add_UCD_predictions(commandArgs(trailingOnly = TRUE))
