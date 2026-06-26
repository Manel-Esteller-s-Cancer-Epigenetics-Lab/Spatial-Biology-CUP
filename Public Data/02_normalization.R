library(tidyverse)
library(Seurat)

projectPath <- "02_ST_CUPS"

source(file.path(projectPath, "misc", "project_utilities.R"))

normalizeData <- function(args) {
    
    sampleSubtype <- args[1]
    sampleId <- args[2]
    baseDir  <- "01_analysis_by_sample"
    step  <- "02_normalization" 
  
    dirsAndLog <- setupDirectoriesAndLog(projectPath, baseDir, sampleSubtype, sampleId, step)
    logFile <- dirsAndLog$logFile
    saveProcessedData <- dirsAndLog$saveProcessedData

    dataPath <- file.path(projectPath, "processed_data", baseDir, "01_preprocessing_qc", sampleSubtype, sampleId, "data.RDS")
    if (!file.exists(dataPath)) {
        stop("Data file does not exist.")
    }
    data <- readRDS(dataPath)

    data <- SCTransform(data, assay = "Spatial", vst.flavor = "v2")

    saveRDS(data, file = file.path(saveProcessedData, "data.RDS"))
}

normalizeData(commandArgs(trailingOnly = TRUE))

