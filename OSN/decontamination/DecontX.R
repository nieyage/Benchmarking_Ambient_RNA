library(Seurat)
library(decontX)
library(tibble)
library(dplyr)
library(SingleCellExperiment)

samples <- list(
  "NC_1" = list(
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/NC_1_filtered.rds",
    counts_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_1_merged/GSM5608794/outs/filtered_feature_bc_matrix/",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/09_DecontX/NC_1/"
  ),
  "NC_2" = list(
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/NC_2_filtered.rds",
    counts_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_2_merged/GSM5608795/outs/filtered_feature_bc_matrix/",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/09_DecontX/NC_2/"
  ),
  "NC_3" = list(
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/NC_3_filtered.rds",
    counts_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_3_merged/GSM5608796/outs/filtered_feature_bc_matrix/",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/09_DecontX/NC_3/"
  ),
  "NC_4" = list(
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/NC_4_filtered.rds",
    counts_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_4_merged/GSM5608797/outs/filtered_feature_bc_matrix/",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/09_DecontX/NC_4/"
  ),
  "SA_3" = list(
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/SA_3_filtered.rds",
    counts_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752985/cellranger/GSM4752985/outs/filtered_feature_bc_matrix/",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/09_DecontX/SA_3/"
  ),
  "SA_4" = list(
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/SA_4_filtered.rds",
    counts_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752986/cellranger/GSM4752986/outs/filtered_feature_bc_matrix/",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/09_DecontX/SA_4/"
  ),
  "SA_5" = list(
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/SA_5_filtered.rds",
    counts_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752989/cellranger/GSM4752989/outs/filtered_feature_bc_matrix/",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/09_DecontX/SA_5/"
  ),
  "SA_6" = list(
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/SA_6_filtered.rds",
    counts_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752990/cellranger/GSM4752990/outs/filtered_feature_bc_matrix/",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/09_DecontX/SA_6/"
  )
)

run_contamination_analysis <- function(seurat_obj, sample_name) {
  cat("Running contamination analysis for sample:", sample_name, "\n")
  
  nonOSNs <- colnames(seurat_obj)[!seurat_obj$celltype %in% c("Mature OSNs")]
  mOSNs <- colnames(seurat_obj)[seurat_obj$celltype == "Mature OSNs"]
  
  nonOSNs_counts <- LayerData(seurat_obj, assay = "RNA", layer = "counts")[, nonOSNs]
  mOSNs_counts <- LayerData(seurat_obj, assay = "RNA", layer = "counts")[, mOSNs]
  
  nonOSNs_decontX <- decontX(nonOSNs_counts)
  seurat_obj$Contamination_nonOSNs <- NA
  seurat_obj$Contamination_nonOSNs[nonOSNs] <- nonOSNs_decontX$contamination
  
  mOSNs_decontX <- decontX(mOSNs_counts)
  seurat_obj$Contamination_mOSNs <- NA
  seurat_obj$Contamination_mOSNs[mOSNs] <- mOSNs_decontX$contamination
  
  contamination_df <- seurat_obj[[]] %>%
    rownames_to_column(var = "cell_id") %>%
    select(cell_id, celltype, Contamination_nonOSNs, Contamination_mOSNs) %>%
    mutate(Contamination = ifelse(!is.na(Contamination_nonOSNs), Contamination_nonOSNs, Contamination_mOSNs))
  
  nonOSNs_stats <- seurat_obj[[]] %>%
    rownames_to_column(var = "cell_id") %>%
    filter(cell_id %in% nonOSNs) %>%
    summarise(group = "nonOSNs", n_cells = n(), mean_contam = mean(Contamination_nonOSNs, na.rm = TRUE), median_contam = median(Contamination_nonOSNs, na.rm = TRUE))
  
  mOSNs_stats <- seurat_obj[[]] %>%
    rownames_to_column(var = "cell_id") %>%
    filter(cell_id %in% mOSNs) %>%
    summarise(group = "mOSNs", n_cells = n(), mean_contam = mean(Contamination_mOSNs, na.rm = TRUE), median_contam = median(Contamination_mOSNs, na.rm = TRUE))
  
  return(list(seurat_obj = seurat_obj, contamination_df = contamination_df, group_summary = bind_rows(nonOSNs_stats, mOSNs_stats)))
}

run_decontamination_correction <- function(counts_path, output_dir, sample_name) {
  cat("Running decontamination correction for sample:", sample_name, "\n")
  
  counts <- Read10X(counts_path)
  sce <- SingleCellExperiment(list(counts = counts))
  sce <- decontX(sce)
  
  return(decontXcounts(sce))
}

for (sample_name in names(samples)) {
  cat("\nProcessing sample:", sample_name, "\n")
  info <- samples[[sample_name]]
  
  if (!dir.exists(info$output_dir)) dir.create(info$output_dir, recursive = TRUE)
  
  tryCatch({
    cat("Loading Seurat object...\n")
    seurat_obj <- readRDS(info$seurat_path)
    
    cat("Running contamination analysis...\n")
    results <- run_contamination_analysis(seurat_obj, sample_name)
    
    write.csv(results$contamination_df, file.path(info$output_dir, paste0(sample_name, "_cell_contamination_details.csv")), row.names = FALSE)
    write.csv(results$group_summary, file.path(info$output_dir, paste0(sample_name, "_group_contamination_summary.csv")), row.names = FALSE)
    saveRDS(results$seurat_obj, file.path(info$output_dir, paste0(sample_name, "_with_contamination.rds")))
    
    cat("Running decontamination correction...\n")
    corrected_counts <- run_decontamination_correction(info$counts_path, info$output_dir, sample_name)
    saveRDS(corrected_counts, file.path(info$output_dir, paste0(sample_name, "_decontX_corrected.rds")))
    
    cat("Completed:", sample_name, "\n")
    
  }, error = function(e) {
    cat("Error in", sample_name, ":", e$message, "\n")
  })
}
