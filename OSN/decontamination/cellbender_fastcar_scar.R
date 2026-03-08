library(Seurat)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(readr)
library(ggplot2)
library(Matrix)
library(qlcMatrix)
library(FastCAR)
library(reticulate)
library(SiftCell)
library(rhdf5)

proj_dir_base <- "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking"
samples <- list(
  list(name="NC1", dir="01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_1_merged/GSM5608794", rds="20_8samples/8sample_rds/NC_1.rds"),
  list(name="NC2", dir="01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_2_merged/GSM5608795", rds="20_8samples/8sample_rds/NC_2.rds"),
  list(name="NC3", dir="01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_3_merged/GSM5608796", rds="20_8samples/8sample_rds/NC_3.rds"),
  list(name="NC4", dir="01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_4_merged/GSM5608797", rds="20_8samples/8sample_rds/NC_4.rds"),
  list(name="SA3", dir="01_rawdata/GSE157100/GSM4752985/cellranger/GSM4752985", rds="20_8samples/8sample_rds/SA_3.rds"),
  list(name="SA4", dir="01_rawdata/GSE157100/GSM4752986/cellranger/GSM4752986", rds="20_8samples/8sample_rds/SA_4.rds"),
  list(name="SA5", dir="01_rawdata/GSE157100/GSM4752989/cellranger/GSM4752989", rds="20_8samples/8sample_rds/SA_5.rds"),
  list(name="SA6", dir="01_rawdata/GSE157100/GSM4752990/cellranger/GSM4752990", rds="20_8samples/8sample_rds/SA_6.rds")

)


for (s in samples) {
  cat("Processing sample:", s$name, "\n")
    
  cat("Load data","\n")
  filtered_path <- file.path(proj_dir_base, s$dir, "outs/filtered_feature_bc_matrix")
  raw_path <- file.path(proj_dir_base, s$dir, "outs/raw_feature_bc_matrix")
  rds_path <- file.path(proj_dir_base, s$rds)
    
  filtered_GEX <- Read10X(filtered_path)
  raw_GEX <- Read10X(raw_path)
  cellranger_seuratObject <- readRDS(rds_path)
  
  ###########################
  # load cells and separate #
  ###########################

  # set ident as celltype
  cellranger_seuratObject <- SetIdent(cellranger_seuratObject, value = "celltype")
  unique(Idents(cellranger_seuratObject))

  # withdraw OSN and non-OSN seurat object
  non_osn_cells <- subset(cellranger_seuratObject, idents = "Mature OSNs", invert = TRUE) 
  osn_cells <- subset(cellranger_seuratObject, idents = "Mature OSNs") 

  # get count matrix
  # non_osn_cells_counts <- GetAssayData(non_osn_cells, slot = "counts")
  # osn_cells_counts <- GetAssayData(osn_cells, slot = "counts")

  # if seurat v5 
  # gene x cell
  non_osn_cells_counts <- GetAssayData(non_osn_cells[["RNA"]], layer = "counts")
  osn_cells_counts <- GetAssayData(osn_cells[["RNA"]], layer = "counts")

  #===================================cellbender====================
  cat("CellBender","\n")

  # conda activate cb
  # 150 epochs, and use default settings

  # 1. run cellbender using bash cli
  # please see cellbender.sh

  # 2. load results
  # the h5 file also includes the corrected matrix

  result_dir <- file.path("/data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/",s$name,"/cellbender_output_file.h5")

  h5ls(result_dir) 
  metadata <- h5read(result_dir, "/metadata") 
  all_barcode <- metadata$barcodes_analyzed
  droplet_latents <- h5read(result_dir, "/droplet_latents") 
  contamination_ratio <- data.frame(background_fraction = droplet_latents$background_fraction,
                                        cell_probability = droplet_latents$cell_probability)
  rownames(contamination_ratio) <- all_barcode
  contamination_ratio_filtered <- contamination_ratio[contamination_ratio$cell_probability > 0.5,]

  # 3. calculate contam_ratio per-cell
  get_contamination_rate <- function(cells, contamination_df) {
  barcodes <- colnames(cells)
    vals <- contamination_df[barcodes, "background_fraction"]
    names(vals) <- barcodes
    return(vals)
  }
  cellbender_non_osn_contamination_rate_percell <- get_contamination_rate(non_osn_cells, contamination_ratio_filtered)
  cellbender_osn_contamination_rate_percell <- get_contamination_rate(osn_cells, contamination_ratio_filtered)

  df_non_osn <- data.frame(
    cell_id = colnames(non_osn_cells),
    cellbender = cellbender_non_osn_contamination_rate_percell,
    celltype = "non-OSN"
  )
  df_osn <- data.frame(
    cell_id = colnames(osn_cells),
    cellbender = cellbender_osn_contamination_rate_percell,
    celltype = "OSN"
  )
  df_merged <- rbind(df_non_osn, df_osn)
  write.csv(df_merged, file.path("/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/",s$name,"/cellbender_contamination_percell.csv"), row.names = FALSE)
  cat("Done!","\n")

  #===============================FastCAR============================
  cat("FastCAR","\n")
  ambProfile = describe.ambient.RNA.sequence(fullCellMatrix = raw_GEX, start = 10, stop = 500, by = 10, contaminationChanceCutoff = 0.05)
  emptyDropletCutoff = recommend.empty.cutoff(ambProfile) 
  contaminationChanceCutoff = 0.05
  ambientProfile = determine.background.to.remove(raw_GEX, emptyDropletCutoff, contaminationChanceCutoff)

  # find ambient genes
  ambient_genes <- names(ambientProfile[ambientProfile > 0])
  ambient_genes

  # Select different calculation methods based on the quantity of ambient_genes
  if(length(ambient_genes) <= 1) {
        fastcar_non_osn_contamination_rate_percell <- non_osn_cells_counts[ambient_genes, ] / colSums(non_osn_cells_counts)
        fastcar_osn_contamination_rate_percell <- osn_cells_counts[ambient_genes, ] / colSums(osn_cells_counts)
  } else {
        fastcar_non_osn_contamination_rate_percell <- colSums(non_osn_cells_counts[ambient_genes, ]) / colSums(non_osn_cells_counts)
        fastcar_osn_contamination_rate_percell <- colSums(osn_cells_counts[ambient_genes, ]) / colSums(osn_cells_counts)
  }

  # save contamination_rate_percell
  df_non_osn <- data.frame(barcodes = names(fastcar_non_osn_contamination_rate_percell), fastcar = as.numeric(fastcar_non_osn_contamination_rate_percell), celltype = "non-OSN",stringsAsFactors = FALSE)
  df_osn <- data.frame(barcodes = names(fastcar_osn_contamination_rate_percell), fastcar = as.numeric(fastcar_osn_contamination_rate_percell), celltype = "OSN", stringsAsFactors = FALSE)
  df_fastcar <- rbind(df_osn, df_non_osn)
  write.csv(df_fastcar, file.path("/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/", s$name, "/fastcar_contamination_rate_percell.csv"), row.names = FALSE)

  # save corrected matrix to h5ad
  fastcar_cellMatrix = remove.background(filtered_GEX, ambientProfile)
  fastcar_seurat_obj <- CreateSeuratObject(fastcar_cellMatrix)
  sceasy::convertFormat(fastcar_seurat_obj, from = "seurat", to = "anndata", outFile = file.path("/data/R04/zhangchao/joint_nuclei/03_analysis/FastCAR/", paste0(s$name, "_corrected_mat.h5ad")))
  cat("Done!","\n")
   
   
  #=================scar================================
  # run scar first

  cat("scAR","\n")

  use_condaenv("/data/R04/zhangchao/anaconda3/envs/scar", required = TRUE)
  py_config()

  file_path <- file.path("/data/R04/zhangchao/joint_nuclei/03_analysis/scar/output/", s$name, "/filtered_feature_bc_matrix_denoised_mRNA.h5ad")

  cmd <- paste(
    "import scanpy as sc",
    sprintf("adata = sc.read_h5ad('%s')", file_path),
    "all_contam = adata.obs['noise_ratio']",
    sep = "\n"
  )
  # cat(cmd)
  py_run_string(cmd)


  # py$all_contam

  scar_non_osn_contamination_rate_percell <- py$all_contam[colnames(non_osn_cells)]
  scar_osn_contamination_rate_percell <- py$all_contam[colnames(osn_cells)]

  df_non_osn <- data.frame(
    cell_id = colnames(non_osn_cells),
    scar = scar_non_osn_contamination_rate_percell,
    celltype = "non-OSN"
  )
  df_osn <- data.frame(
    cell_id = colnames(osn_cells),
    scar = scar_osn_contamination_rate_percell,
    celltype = "OSN"
  )
  df_merged <- rbind(df_non_osn, df_osn)
  write.csv(df_merged, file.path("/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/",s$name,"/scar_contamination_percell.csv"), row.names = FALSE)
  cat("Done!","\n")
}
