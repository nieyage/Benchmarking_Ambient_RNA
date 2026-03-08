library(SoupX)
library(Seurat)
library(Matrix)
library(dplyr)

samples <- list(
  "NC_1" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_1_merged/GSM5608794/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/NC_1_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/08_SoupX/NC_1/"
  ),
  "NC_2" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_2_merged/GSM5608795/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/NC_2_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/08_SoupX/NC_2/"
  ),
  "NC_3" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_3_merged/GSM5608796/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/NC_3_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/08_SoupX/NC_3/"
  ),
  "NC_4" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_4_merged/GSM5608797/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/NC_4_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/08_SoupX/NC_4/"
  ),
  "SA_3" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752985/cellranger/GSM4752985/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/SA_3_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/08_SoupX/SA_3/"
  ),
  "SA_4" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752986/cellranger/GSM4752986/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/SA_4_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/08_SoupX/SA_4/"
  ),
  "SA_5" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752989/cellranger/GSM4752989/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/SA_5_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/08_SoupX/SA_5/"
  ),
  "SA_6" = list(
    raw_data_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752990/cellranger/GSM4752990/outs/",
    seurat_path = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/19_11sample_new/split_samples/SA_6_filtered.rds",
    output_dir = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/08_SoupX/SA_6/"
  )
)

process_all_samples <- function() {
  for (sample_name in names(samples)) {
    cat("Processing", sample_name, "...\n")
    
    info <- samples[[sample_name]]
    if (!dir.exists(info$output_dir)) dir.create(info$output_dir, recursive = TRUE)
    
    tryCatch({
      # SoupX decontamination
      data_dir <- info$raw_data_path
      tod <- Read10X(file.path(data_dir, "raw_feature_bc_matrix"))
      toc <- Read10X(file.path(data_dir, "filtered_feature_bc_matrix"))
      
      common <- intersect(rownames(tod), rownames(toc))
      tod <- as(tod[common, ], "dgCMatrix")
      toc <- as(toc[common, ], "dgCMatrix")
      
      # Clustering
      so <- CreateSeuratObject(counts = toc) %>%
        NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>%
        RunPCA() %>% FindNeighbors() %>% FindClusters(resolution = 0.5)
      
      # SoupChannel
      sc <- SoupChannel(tod, toc, calcSoupProfile = TRUE) %>%
        setClusters(Idents(so)) %>%
        autoEstCont()
      
      corrected <- adjustCounts(sc, roundToInt = TRUE)
      saveRDS(corrected, file.path(info$output_dir, paste0(sample_name, "_soupx_corrected_counts.rds")))
      cat("  Saved corrected counts\n")
      
      # Calculate pollution rates
      original_obj <- readRDS(info$seurat_path)
      
      corrected_mat <- if (inherits(corrected, "Matrix")) corrected else GetAssayData(corrected, slot = "counts")
      
      if ("RNA" %in% names(original_obj@assays)) {
        original_mat <- GetAssayData(original_obj, assay = "RNA", slot = "counts")
      } else {
        original_mat <- GetAssayData(original_obj, slot = "counts")
      }
      
      # Match barcodes
      common_cells <- intersect(colnames(corrected_mat), colnames(original_mat))
      if (length(common_cells) == 0) {
        orig_bc <- colnames(original_mat)
        corr_bc <- colnames(corrected_mat)
        
        if (all(grepl("-1$", corr_bc)) && !all(grepl("-1$", orig_bc))) {
          colnames(original_mat) <- paste0(orig_bc, "-1")
        } else if (!all(grepl("-1$", corr_bc)) && all(grepl("-1$", orig_bc))) {
          colnames(original_mat) <- gsub("-1$", "", orig_bc)
        }
        common_cells <- intersect(colnames(corrected_mat), colnames(original_mat))
      }
      
      if (length(common_cells) > 0) {
        corrected_mat <- corrected_mat[, common_cells, drop = FALSE]
        original_mat <- original_mat[, common_cells, drop = FALSE]
        
        original_umi <- colSums(original_mat)
        corrected_umi <- colSums(corrected_mat)
        
        celltype <- if ("celltype" %in% names(original_obj@meta.data)) {
          original_obj@meta.data[common_cells, "celltype"]
        } else {
          as.character(Idents(original_obj)[common_cells])
        }
        
        df <- data.frame(
          cell_id = common_cells,
          celltype = celltype,
          original_umi = original_umi,
          corrected_umi = corrected_umi,
          pollution_rate = 1 - (corrected_umi / original_umi),
          sample = sample_name
        )
        
        write.csv(df, file.path(info$output_dir, paste0(sample_name, "_percell_pollution_rates.csv")), row.names = FALSE)
        
        cat("  Success:", sample_name, "- Cells:", length(common_cells),
            "Mean pollution:", round(mean(df$pollution_rate), 4), "\n\n")
      } else {
        cat("  No common cells found\n\n")
      }
      
    }, error = function(e) {
      cat("  Error:", e$message, "\n\n")
    })
  }
}

process_all_samples()
```
