library(tidyverse)
library(Seurat)

projectPath <- "02_ST_CUPS"  
source(file.path(projectPath, "misc", "project_utilities.R"))

loadData <- function(sampleSubtype, sampleId, logFile) {
     dataDir <- file.path(projectPath, "raw_data", sampleSubtype, sampleId)
    
    if (!dir.exists(dataDir)) {
        stop("Data not found: ", dataDir)
    }
        
    # Handle samples fron Masood Nat Genetics 2024 as they are directly data.RDS (previously processed using reformat_samples_LNM_Masood-NatGen2024.R)
    if(sampleId %in% c("IU_PDA_LNM12", "IU_PDA_LNM10", "IU_PDA_LNM8", "IU_PDA_LNM7", "IU_PDA_LNM6")){
        data <- readRDS(file.path(dataDir, "data.RDS"))
	data <- VisiumV2_to_VisiumV1(se_obj_V2 = data, path_to_tissue_positions = file.path(dataDir, "spatial", "tissue_positions.csv"))
    } else{

        # Handle primary samples
        if (sampleSubtype == "primaries") {
            data <- readRDS(file.path(dataDir, "data.RDS"))
        } else {
            # Determine the correct file to load based on the sampleId
            matrixFile <- if (sampleId %in% c("GSM5732357", "GSM5732358", "GSM5732359", "GSM5732360")) {
                "raw_feature_bc_matrix.h5"
            } else {
                "filtered_feature_bc_matrix.h5"
            }
            
            # Load the spatial data with Seurat
            data <- Seurat::Load10X_Spatial(
                data.dir = dataDir,
                filename = matrixFile,
                assay = "Spatial",
                slice = sampleId,
                filter.matrix = TRUE,
                to.upper = FALSE,
                image = NULL
            )

        # For this samples, which only has matrixFile raw_feature_bc_matrix.h5 available, subset the barcodes so they match the slice@coordinates barcodes (which are fewer, as in the slice@coordinates they are filtered)
        if (sampleId %in% c("GSM5732357", "GSM5732358", "GSM5732359", "GSM5732360")){           
        	slice <- data@images[[1]]
          
        	# Extract barcodes from slice@coordinates
        	barcodes_in_coordinates <- rownames(slice@coordinates)
          
        	# Extract barcodes from Seurat object (colnames of the count matrix)
        	barcodes_in_assay <- colnames(data)
          
        	# Filter the barcodes to keep only those present in both Seurat object and slice coordinates
        	common_barcodes <- intersect(barcodes_in_assay, barcodes_in_coordinates)
        	data <- subset(data, cells = common_barcodes)
        }
            
            # Check if the low-res image is missing and replace it with hi-res scale factors
            spatialDir <- file.path(projectPath, "raw_data", sampleSubtype, sampleId, "spatial")
            if (!"tissue_hires_image.png" %in% list.files(spatialDir)) {
                data@images[[1]]@scale.factors$lowres <- data@images[[1]]@scale.factors$hires
            }
        }
    }

        return(data)
}


# Perform QC
performQC <- function(data, sampleId, savePlots, logFile) {
    logMessage("Performing quality control", logFile)
    logMessage(paste("Number of spots before filtering: ", ncol(data)), logFile)
    logMessage(paste("Number of genes before filtering: ", nrow(data)), logFile)
    data <- subset(data, subset = nFeature_Spatial > 100)
    
    if(sampleId %in% c("GSM5420749_pt15", "GSM5420750_pt16", "GSM5420751_pt19", "GSM5420752_pt24", "GSM5420753_pt26", "GSM5420754_pt27")){
        data <- subset(data, subset = percent.mt < 50)
    }else{
        data <- subset(data, subset = percent.mt < 20)
    }

    logMessage(paste("Number of spots after filtering: ", ncol(data)), logFile)
    gene_counts <- rowSums(data@assays$Spatial@layers$counts)
    filtered_genes <- names(gene_counts[gene_counts >= 5])
    data <- data[filtered_genes, ]
    logMessage(paste("Number of genes after filtering: ", nrow(data)), logFile)
    logMessage("Filtering performed: 1) nFeature_Spatial > 100; 2) gene_counts >= 5", logFile)
    return(data)
}


plotQC <- function(data, sampleId, savePlots, logFile, step) {
    logMessage(paste("Ploting QC metrics ", step, " filtering"), logFile)
    pt.size.factor <- get_pt_size_factor(sampleId)
    variables <- c("nCount_Spatial","nFeature_Spatial", "percent.mt") %>% setNames(c("nCount_Spatial","nFeature_Spatial", "percent.mt"))
    qc_vln_before <- variables %>% purrr::map(~VlnPlot(data, features = .x, pt.size = 0.1) + NoLegend())
    qc_sf_before <- variables %>% purrr::map(~SpatialFeaturePlot(data, features = .x, pt.size.factor = pt.size.factor, stroke = 0.0))
    qc_before <- qc_vln_before %>% purrr::map2(qc_sf_before, ~patchwork::wrap_plots(.x, .y))
    qc_before %>% purrr::imap(~ggsave(.x, file = file.path(savePlots, paste0("QC_", step, "_", .y, ".pdf")), width = 8, height = 6))
    pdf(file = file.path(savePlots,paste0(step, "_nFeature_density.pdf")))
    plot(density(data$nFeature_Spatial))
    abline(v = 100)
    dev.off()
}

saveResults <- function(data, savePath, logFile) {
    logMessage("Saving processed data", logFile)
    saveRDS(data, file = file.path(savePath, "data.RDS"))
}

preprocessing_visium <- function(args) {

    sampleSubtype <- args[1]
    sampleId <- args[2]

    dirsAndLog <- setupDirectoriesAndLog(projectPath, "01_analysis_by_sample", sampleSubtype, sampleId, "01_preprocessing_qc")
    logFile <- dirsAndLog$logFile
    savePlots <- dirsAndLog$savePlots
    saveProcessedData <- dirsAndLog$saveProcessedData

    data <- loadData(sampleSubtype, sampleId, logFile)
    data <- subset(data, subset = nCount_Spatial > 0)
    data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
    plotQC(data, sampleId, savePlots, logFile, "before")
    data <- performQC(data, sampleId, savePlots, logFile)
    plotQC(data, sampleId, savePlots, logFile, "after")
    saveResults(data, saveProcessedData, logFile)
    
}

preprocessing_visium(commandArgs(trailingOnly = TRUE))

