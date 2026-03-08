library(Seurat)
library(scCDC)
library(dplyr)
library(readr)

samples <- list(
  "NC_1" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_1_merged/GSM5608794/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/NC_1_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/05_scCDC_contamination_ratio/11samples_1126/NC_1/"
  ),
  "NC_2" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_2_merged/GSM5608795/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/NC_2_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/05_scCDC_contamination_ratio/11samples_1126/NC_2/"
  ),
  "NC_3" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_3_merged/GSM5608796/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/NC_3_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/05_scCDC_contamination_ratio/11samples_1126/NC_3/"
  ),
  "NC_4" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_4_merged/GSM5608797/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/NC_4_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/05_scCDC_contamination_ratio/11samples_1126/NC_4/"
  ),
  "SA_3" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752985/cellranger/GSM4752985/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/SA_3_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/05_scCDC_contamination_ratio/11samples_1126/SA_3/"
  ),
  "SA_4" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752986/cellranger/GSM4752986/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/SA_4_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/05_scCDC_contamination_ratio/11samples_1126/SA_4/"
  ),
  "SA_5" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752989/cellranger/GSM4752989/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/SA_5_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/05_scCDC_contamination_ratio/11samples_1126/SA_5/"
  ),
  "SA_6" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752990/cellranger/GSM4752990/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/SA_6_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/05_scCDC_contamination_ratio/11samples_1126/SA_6/"
  )
)

process_single_sample <- function(seurat_path, output_dir, sample_name) {
  cat("Processing:", sample_name, "\n")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  object <- readRDS(seurat_path)
  
  if ("RNA" %in% Assays(object)) {
    layers <- Layers(object, assay = "RNA")
    if ("counts" %in% layers) {
      counts_matrix <- LayerData(object, assay = "RNA", layer = "counts")
    } else if ("data" %in% layers) {
      counts_matrix <- LayerData(object, assay = "RNA", layer = "data")
    } else {
      counts_matrix <- GetAssayData(object, assay = "RNA", slot = "counts")
      if (is.null(counts_matrix)) counts_matrix <- GetAssayData(object, assay = "RNA", slot = "data")
    }
  } else {
    assay_name <- Assays(object)[1]
    counts_matrix <- GetAssayData(object, assay = assay_name, slot = "counts")
    if (is.null(counts_matrix)) counts_matrix <- GetAssayData(object, assay = assay_name, slot = "data")
  }
  
  v4_object <- CreateSeuratObject(counts = counts_matrix, meta.data = object[[]])
  
  if ("RNA" %in% Assays(object) && "data" %in% Layers(object, assay = "RNA")) {
    v4_object[["RNA"]]$data <- LayerData(object, assay = "RNA", layer = "data")
  }
  
  if ("celltype" %in% names(v4_object@meta.data)) {
    Idents(v4_object) <- v4_object$celltype
  } else if ("seurat_clusters" %in% names(v4_object@meta.data)) {
    Idents(v4_object) <- v4_object$seurat_clusters
  }
  
  GCGs <- ContaminationDetection(v4_object)
  contaminated_genes <- rownames(GCGs)
  
  overall_contamination <- ContaminationQuantification(v4_object, contaminated_genes)
  
  v4_object_corrected <- ContaminationCorrection(v4_object, contaminated_genes)
  DefaultAssay(v4_object_corrected) <- "Corrected"
  
  original_counts <- GetAssayData(v4_object, layer = "counts")
  total_umi_per_cell <- Matrix::colSums(original_counts)
  contam_umi_per_cell <- Matrix::colSums(original_counts[contaminated_genes, , drop = FALSE])
  
  all_cells_results <- data.frame(
    cell_id = names(contamination_ratio),
    total_umi = total_umi_per_cell,
    contam_umi = contam_umi_per_cell,
    contamination_ratio = contam_umi_per_cell / total_umi_per_cell,
    celltype = if ("celltype" %in% names(v4_object@meta.data)) {
      v4_object@meta.data[names(contamination_ratio), "celltype"]
    } else if ("seurat_clusters" %in% names(v4_object@meta.data)) {
      paste0("Cluster_", v4_object@meta.data[names(contamination_ratio), "seurat_clusters"])
    } else {
      "Unknown"
    },
    stringsAsFactors = FALSE
  )
  
  celltype_stats <- all_cells_results %>%
    group_by(celltype) %>%
    summarise(
      n_cells = n(),
      mean_contamination = mean(contamination_ratio),
      median_contamination = median(contamination_ratio),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_contamination))
  
  saveRDS(v4_object_corrected, file.path(output_dir, paste0(sample_name, "_scCDC_corrected.rds")))
  write_tsv(all_cells_results, file.path(output_dir, paste0(sample_name, "_per_cell_contamination_rates.tsv")))
  write_tsv(data.frame(gene = contaminated_genes), file.path(output_dir, paste0(sample_name, "_contaminated_genes.tsv")))
  write_tsv(data.frame(GCGs), file.path(output_dir, paste0(sample_name, "_GCGs_detailed.tsv")))
  write_tsv(celltype_stats, file.path(output_dir, paste0(sample_name, "_contamination_by_celltype.tsv")))
  
  overall_stats <- data.frame(
    sample = sample_name,
    total_cells = ncol(v4_object),
    n_contaminated_genes = length(contaminated_genes),
    overall_contamination_ratio = overall_contamination,
    mean_per_cell_contamination = mean(all_cells_results$contamination_ratio)
  )
  write_tsv(overall_stats, file.path(output_dir, paste0(sample_name, "_overall_contamination_stats.tsv")))
  
  return(list(success = TRUE, n_cells = ncol(v4_object), n_genes = length(contaminated_genes), overall_contamination = overall_contamination))
}

process_all_samples <- function(samples) {
  cat("Processing", length(samples), "samples...\n")
  start_time <- Sys.time()
  
  results <- list()
  successful <- c()
  failed <- c()
  
  for (name in names(samples)) {
    cat("\n---", name, "---\n")
    tryCatch({
      res <- process_single_sample(samples[[name]]$seurat_path, samples[[name]]$output_dir, name)
      if (res$success) {
        successful <- c(successful, name)
        results[[name]] <- res
        cat("SUCCESS:", name, "- Cells:", res$n_cells, "Contam genes:", res$n_genes, 
            "Overall contam:", round(res$overall_contamination, 4), "\n")
      }
    }, error = function(e) {
      cat("FAILED:", name, "-", e$message, "\n")
      failed <- c(failed, name)
      results[[name]] <- list(success = FALSE)
    })
  }
}

final_results <- process_all_samples(samples)
```
