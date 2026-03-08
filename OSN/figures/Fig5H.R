library(Seurat)
library(tidyverse)
library(readxl)
library(pheatmap)
library(Matrix)

outdir <- "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/OR_heatmap/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

specific_or_genes <- c(
    "Olfr1123", "Olfr1140", "Olfr117", "Olfr1297", "Olfr1348", 
    "Olfr1384", "Olfr1437", "Olfr1440", "Olfr17", "Olfr195", 
    "Olfr222", "Olfr24", "Olfr279", "Olfr403", "Olfr523", "Olfr536", 
    "Olfr686", "Olfr71", "Olfr728", "Olfr794", "Olfr827"
)

rds_list <- list(
  "Original" = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/integrated_8samples_filtered.rds",
  "CellClear" = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/cellclear_umap/final_integrated_8samples_cellclear.rds",
  "SCAR" = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/scar_umap/integrated_8samples_clustered.rds"
)

fix_bc <- function(bc) substr(bc, nchar(bc)-15, nchar(bc))
get_merged_matrix <- function(seurat_obj, celltype_value="Mature OSNs"){
  data_layers <- grep("^data\\.", Layers(seurat_obj[["RNA"]]), value=TRUE)
  
  merged_data <- do.call(cbind, lapply(data_layers, function(lay){
    LayerData(seurat_obj, assay="RNA", layer=lay)
  }))
  
  mature_idx <- which(seurat_obj$celltype == celltype_value)
  merged_data_mature <- merged_data[, mature_idx, drop=FALSE]
  
  short_bc <- fix_bc(colnames(merged_data_mature))
  
  list(expr=merged_data_mature, short_bc=short_bc)
}

get_OR_matrix <- function(seurat_obj, bc_list, genes){
  data_layers <- grep("^data\\.", Layers(seurat_obj[["RNA"]]), value=TRUE)
  merged_data <- do.call(cbind, lapply(data_layers, function(lay){
    LayerData(seurat_obj, assay="RNA", layer=lay)
  }))
  
  short_bc <- fix_bc(colnames(merged_data))
  keep_idx <- which(short_bc %in% bc_list)
  
  expr <- merged_data[genes, keep_idx, drop=FALSE]
  expr <- expr[, order(match(short_bc[keep_idx], bc_list))]
  
  return(expr)
}

orig <- readRDS(rds_list$Original)
merged_orig <- get_merged_matrix(orig, celltype_value="Mature OSNs")
mat_orig <- merged_orig$expr
mature_bc <- merged_orig$short_bc

available_genes <- intersect(specific_or_genes, rownames(mat_orig))
cat("Available OR genes in Original:", length(available_genes), "\n")
cat("Missing genes:", setdiff(specific_or_genes, rownames(mat_orig)), "\n")

mat_orig_filtered <- mat_orig[available_genes, , drop=FALSE]
cat("Filtered matrix dimensions:", dim(mat_orig_filtered), "\n")

cellclear <- readRDS(rds_list$CellClear)
mat_cellclear <- get_OR_matrix(cellclear, mature_bc, available_genes)

scar <- readRDS(rds_list$SCAR)
mat_scar <- get_OR_matrix(scar, mature_bc, available_genes)

writeLines(available_genes, file.path(outdir, "selected_or_genes.txt"))
cat("Saved gene list to:", file.path(outdir, "selected_or_genes.txt"), "\n")

# Original heatmap
pdf(file.path(outdir, "OR_heatmap_Original_specific.pdf"), width=12, height=6)
pheatmap(mat_orig_filtered, 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         main=paste("Original - Specific OR Genes (", length(available_genes), "genes, ", 
                   ncol(mat_orig_filtered), "cells)", sep=""),
         show_colnames=FALSE,
         fontsize_row=8)
dev.off()

# CellClear heatmap
pdf(file.path(outdir, "OR_heatmap_CellClear_specific.pdf"), width=12, height=6)
pheatmap(mat_cellclear, 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         main=paste("CellClear - Specific OR Genes (", length(available_genes), "genes, ", 
                   ncol(mat_cellclear), "cells)", sep=""),
         show_colnames=FALSE,
         fontsize_row=8)
dev.off()

# SCAR heatmap
pdf(file.path(outdir, "OR_heatmap_SCAR_specific.pdf"), width=12, height=6)
pheatmap(mat_scar, 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         main=paste("SCAR - Specific OR Genes (", length(available_genes), "genes, ", 
                   ncol(mat_scar), "cells)", sep=""),
         show_colnames=FALSE,
         fontsize_row=8)
dev.off()
