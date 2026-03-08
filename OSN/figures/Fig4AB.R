### Fig4A
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(patchwork)
set.seed(100)

# Create output directory 
out_dir <- "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir), recursive = TRUE)
}

NC_1_path <- "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_1_merged/GSM5608794/outs/"
NC_2_path <- "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_2_merged/GSM5608795/outs/"
NC_3_path <- "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_3_merged/GSM5608796/outs/"
NC_4_path <- "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_4_merged/GSM5608797/outs/"
SA_3_path <- "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752985/cellranger/GSM4752985/outs/"
SA_4_path <- "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752986/cellranger/GSM4752986/outs/"
SA_5_path <- "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752989/cellranger/GSM4752989/outs/"
SA_6_path <- "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752990/cellranger/GSM4752990/outs/"

obj_NC_1 <- CreateSeuratObject(
  counts = Read10X_h5(paste0(NC_1_path, "filtered_feature_bc_matrix.h5")), 
  project = "NC_1", 
  assay = "RNA"
)

obj_NC_2 <- CreateSeuratObject(
  counts = Read10X_h5(paste0(NC_2_path, "filtered_feature_bc_matrix.h5")), 
  project = "NC_2", 
  assay = "RNA"
)

obj_NC_3 <- CreateSeuratObject(
  counts = Read10X_h5(paste0(NC_3_path, "filtered_feature_bc_matrix.h5")), 
  project = "NC_3", 
  assay = "RNA"
)

obj_NC_4 <- CreateSeuratObject(
  counts = Read10X_h5(paste0(NC_4_path, "filtered_feature_bc_matrix.h5")), 
  project = "NC_4", 
  assay = "RNA"
)

obj_SA_3 <- CreateSeuratObject(
  counts = Read10X_h5(paste0(SA_3_path, "filtered_feature_bc_matrix.h5")), 
  project = "SA_3", 
  assay = "RNA"
)

obj_SA_4 <- CreateSeuratObject(
  counts = Read10X_h5(paste0(SA_4_path, "filtered_feature_bc_matrix.h5")), 
  project = "SA_4", 
  assay = "RNA"
)

obj_SA_5 <- CreateSeuratObject(
  counts = Read10X_h5(paste0(SA_5_path, "filtered_feature_bc_matrix.h5")), 
  project = "SA_5", 
  assay = "RNA"
)

obj_SA_6 <- CreateSeuratObject(
  counts = Read10X_h5(paste0(SA_6_path, "filtered_feature_bc_matrix.h5")), 
  project = "SA_6", 
  assay = "RNA"
)

objList <- list(obj_NC_1, obj_NC_2, obj_NC_3, obj_NC_4, 
                obj_SA_3, obj_SA_4, obj_SA_5, obj_SA_6)

samples <- c("NC_1", "NC_2", "NC_3", "NC_4", 
             "SA_3", "SA_4", "SA_5", "SA_6")

# Quality control for all samples
for (i in 1:length(objList)){
  objList[[i]][["percent.mt"]] <- PercentageFeatureSet(objList[[i]], pattern = "^mt-")
  pdf(paste0(out_dir, samples[i], "_qc_plot.pdf"), width = 9)
  p1 <- VlnPlot(objList[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p2 <- ggplot(data = objList[[i]][["nFeature_RNA"]], aes(x = nFeature_RNA)) + geom_density()
  p3 <- ggplot(data = objList[[i]][["nCount_RNA"]], aes(x = nCount_RNA)) + geom_density()
  p4 <- ggplot(data = objList[[i]][["percent.mt"]], aes(x = percent.mt)) + geom_density()
  p5 <- FeatureScatter(objList[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  p6 <- FeatureScatter(objList[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  dev.off()
}

# Filter cells for each sample
objList[[1]] <- subset(objList[[1]], subset = nFeature_RNA >= 1300 & nFeature_RNA <= 7500 & percent.mt < 15)   # NC1
objList[[2]] <- subset(objList[[2]], subset = nFeature_RNA >= 1000 & nFeature_RNA <= 7000 & percent.mt < 15)  # NC2
objList[[3]] <- subset(objList[[3]], subset = nFeature_RNA >= 2000 & nFeature_RNA <= 7500 & percent.mt < 25)  # NC3
objList[[4]] <- subset(objList[[4]], subset = nFeature_RNA >= 1200 & nFeature_RNA <= 6500 & percent.mt < 15)  # NC4
objList[[5]] <- subset(objList[[5]], subset = nFeature_RNA >= 1500 & nFeature_RNA <= 10000 & percent.mt < 25)  # SA3
objList[[6]] <- subset(objList[[6]], subset = nFeature_RNA >= 1500 & nFeature_RNA <= 10000 & percent.mt < 25)   # SA4
objList[[7]] <- subset(objList[[7]], subset = nFeature_RNA >= 1200 & nFeature_RNA <= 11000 & percent.mt < 25)  # SA5
objList[[8]] <- subset(objList[[8]], subset = nFeature_RNA >= 1500 & nFeature_RNA <= 10000 & percent.mt < 25)  # SA6

# Merge objects
merged_obj <- merge(objList[[1]], y = objList[-1], add.cell.ids = samples)

# Split by sample for integration
obj.ls <- SplitObject(merged_obj, split.by = "orig.ident")

# Normalize and find variable features
for (i in 1:length(obj.ls)){
  obj.ls[[i]] <- NormalizeData(obj.ls[[i]], verbose = FALSE)
  obj.ls[[i]] <- FindVariableFeatures(obj.ls[[i]], selection.method = "vst", nfeatures = 3000, verbose = FALSE)
}

# Find integration anchors and integrate
anchors <- FindIntegrationAnchors(object.list = obj.ls, anchor.features = 3000, dims = 1:30)
integrated_scRNA <- IntegrateData(anchorset = anchors, dims = 1:30)

# Switch to integrated assay
DefaultAssay(integrated_scRNA) <- "integrated"

# Scale data
integrated_scRNA <- ScaleData(integrated_scRNA, features = rownames(integrated_scRNA))

# Perform PCA
integrated_scRNA <- RunPCA(integrated_scRNA, npcs = 50, verbose = FALSE)
pdf(paste0(out_dir, "integrated_8samples_pc.pdf"))
ElbowPlot(integrated_scRNA, ndims = 50)
dev.off()

markers <- c(
  "Omp", "Gng13", "Stoml3","Cnga2", # Mature Olfactory Sensory Neurons (mOSNs)
  "Nqo1", "Ncam2", # Dorsal/Ventral OSNs
  "Gap43", "Gng8", # Immature Olfactory Sensory Neurons (imOSNs)
  "Neurog1", "Neurod1", "Sox11", # Immediate Neuronal Precursors (INPs)
  "Ascl1", "Kit", # Globose Basal Cells (GBCs)
  "Krt5", "Krt14", "Trp63", "Cxcl14", # Horizontal Basal Cells (HBCs)
  "Cyp2g1", "Cyp1a2", "Cbr2", "Ermn", "Sox2", # Sustentacular Cells
  "Trpm5", "Sh2d7", # Brush-like MVs
  "Ascl3", "Cftr", # Ionocyte-like MVs
  "Aqp5", "Sox9", "Sox10", "Muc5b", "Muc5ac", # Bowman's Gland
  "Gfap", "Atp1a2", "Fabp7", "S100b", "Plp1", # Olfactory Ensheathing Glia
  "Lum", "Dcn", # Fibroblasts
  "Tagln", "Myh11", # Vascular Smooth Muscle Cells
  "Cd3d", "Cd3e", # T Cells
  "Cd79a", "Cd19", "Cd37", # B Cells
  "S100a8", "S100a9", "Csf3r", # Neutrophils
  "Tpsb2", "Tpsab1", # Mast Cells
  "C1qa", "C1qb", "Ms4a7", "Adgre1", # Macrophages
  "H2-Ab1", "Itgax", # Dendritic Cells
  "Ccr2", "S100a4", "Lyz2", # Monocytes
  "Krt19", # Respiratory Cells
  "Foxj1", "Cfap126", # Respiratory Ciliated Cells
  "Hbb-bs", "Hbb-bt", # Erythrocytes
  "Mcpt8","Ccl4" #Basophils
)

for (nPCs in seq(20, 35, 5)){
  DefaultAssay(integrated_scRNA) <- "integrated"
  integrated_scRNA <- FindNeighbors(integrated_scRNA, dims = 1:nPCs)
  integrated_scRNA <- FindClusters(integrated_scRNA, resolution = 0.2)
  integrated_scRNA <- RunUMAP(integrated_scRNA, dims = 1:nPCs)
  pdf(paste0(out_dir, "integrated_8samples_PC", nPCs, "_resolution0.2_UMAP.pdf"))
  p1 <- DimPlot(integrated_scRNA, reduction = "umap", label = TRUE)
  p2 <- DimPlot(integrated_scRNA, reduction = "umap", group.by = "orig.ident")
  print(p1)
  print(p2)
  dev.off()
  DefaultAssay(integrated_scRNA) <- "RNA"
  pdf(paste0(out_dir, "integrated_8samples_PC", nPCs, "_markers_UMAP.pdf"))
  p <- FeaturePlot(integrated_scRNA, features = markers, combine = FALSE, cols = c("lightgrey", "red"), order = TRUE)
  print(p)
  dev.off()
}

# Select optimal PCs
DefaultAssay(integrated_scRNA) <- "integrated"
integrated_scRNA <- FindNeighbors(integrated_scRNA, dims = 1:25)
integrated_scRNA <- RunUMAP(integrated_scRNA, dims = 1:25)

# Clustering
integrated_scRNA <- FindClusters(integrated_scRNA, resolution = 1.5)

# Visualization
pdf(paste0(out_dir, "umap_clusters_8samples.pdf"), width = 9, height = 7)
DimPlot(integrated_scRNA, reduction = "umap", label = TRUE)
dev.off()

png(paste0(out_dir, "umap_clusters_8samples.png"), width = 9, height = 7, res = 300)
DimPlot(integrated_scRNA, reduction = "umap", label = TRUE)
dev.off()

# Marker gene dot plot
DefaultAssay(integrated_scRNA) <- "RNA"

dot_plot <- DotPlot(integrated_scRNA,
                   features = markers,
                   cols = c("lightgrey", "blue"),
                   dot.scale = 6) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)) +
  ggtitle("Marker Gene Expression DotPlot - 8 Samples") +
  labs(x = "Genes", y = "Clusters")

ggsave(paste0(out_dir, "marker_dotplot_8samples.pdf"), 
       plot = dot_plot, 
       width = 22, height = 15)

ggsave(paste0(out_dir, "marker_dotplot_8samples.png"), 
       plot = dot_plot, 
       width = 22, height = 15, dpi = 300)

# Cell type annotation
DefaultAssay(integrated_scRNA) <- "integrated"

new.cluster.ids <- c(
  "0" = "Mature OSNs",
  "1" = "Mature OSNs", 
  "2" = "Mature OSNs",
  "3" = "Mature OSNs",
  "4" = "Mature OSNs",
  "5" = "Mature OSNs",
  "6" = "Mature OSNs",
  "7" = "Ionocyte-like MVs",
  "8" = "Immature OSNs",
  "9" = "Mature OSNs",
  "10" = "Mature OSNs",
  "11" = "Mature OSNs",
  "12" = "Mature OSNs",
  "13" = "Mature OSNs",
  "14" = "Mature OSNs",
  "15" = "HBCs",
  "16" = "Immature OSNs",
  "17" = "Mature OSNs",
  "18" = "Immature OSNs",
  "19" = "Ionocyte-like MVs",
  "20" = "INPs",
  "21" = "HBCs",
  "22" = "Mature OSNs",
  "23" = "Mature OSNs",
  "24" = "GBCs",
  "25" = "Bowman's Gland",
  "26" = "Ionocyte-like MVs",
  "27" = "Ionocyte-like MVs",
  "28" = "Sustentacular Cells",
  "29" = "Immature OSNs",
  "30" = "Neutrophils",
  "31" = "Monocytes",
  "32" = "Monocytes",
  "33" = "Macrophages",
  "34" = "Neutrophils",
  "35" = "B Cells",
  "36" = "Brush-like MVs",
  "37" = "Respiratory Cells",
  "38" = "Sustentacular Cells",
  "39" = "Erythrocytes",
  "40" = "Olfactory Ensheathing Glia",
  "41" = "dk",
  "42" = "dk",
  "43" = "dk",
  "44" = "B Cells",
  "45" = "HBCs",
  "46" = "dk"
)

integrated_scRNA <- RenameIdents(integrated_scRNA, new.cluster.ids)

celltype_df <- data.frame(
  celltype = as.vector(Idents(integrated_scRNA)), 
  row.names = names(Idents(integrated_scRNA))
)
integrated_scRNA <- AddMetaData(integrated_scRNA, celltype_df)

cat("Cell type distribution before removing dk cells:\n")
print(table(integrated_scRNA$celltype))
integrated_scRNA_clean <- subset(integrated_scRNA, subset = celltype != "dk")
cat("\nCell type distribution after removing dk cells:\n")
print(table(integrated_scRNA_clean$celltype))

celltype_order <- c(
  "HBCs", "GBCs", "INPs", 
  "Immature OSNs", "Mature OSNs",
  "Ionocyte-like MVs", "Bowman's Gland", "Brush-like MVs",
  "Neutrophils", "Macrophages", "Monocytes", "B Cells", "Erythrocytes",
  "Sustentacular Cells", "Respiratory Cells", "Olfactory Ensheathing Glia"
)

existing_celltypes <- unique(integrated_scRNA_clean$celltype)
celltype_order <- celltype_order[celltype_order %in% existing_celltypes]

cat("Final cell type order:\n")
print(celltype_order)

integrated_scRNA_clean$celltype <- factor(integrated_scRNA_clean$celltype, levels = celltype_order)

library(RColorBrewer)

n_celltypes <- length(celltype_order)

color_vector <- c(
  brewer.pal(9, "YlGn")[c(7,6,5,4,3)],
  brewer.pal(9, "RdPu")[c(3,4,5)],
  brewer.pal(9, "Blues")[c(3,4,6,7,8)],
  brewer.pal(9, "Purples")[c(4,6,8)]
)

newpalette <- setNames(color_vector, celltype_order)

cat("Color mapping:\n")
print(newpalette)

p <- DimPlot(integrated_scRNA_clean, 
        reduction = "umap", 
        group.by = "celltype",
        label = TRUE,
        repel = TRUE,
        label.size = 3,
        cols = newpalette,
        order = celltype_order) +
  ggtitle("Integrated scRNA-seq - Cell Type Annotation") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

pdf(paste0(out_dir, "umap_celltype_annotated.pdf"), width = 10, height = 8)
print(p)
dev.off()

png(paste0(out_dir, "umap_celltype_annotated.png"), 
    width = 10, height = 8, res = 300, units = "in")
print(p)
dev.off()

final_obj_clean <- integrated_scRNA_clean

Idents(final_obj_clean) <- "celltype"

saveRDS(final_obj_clean, paste0(out_dir, "integrated_8samples_filtered.rds"))

### Fig4B
pdf(paste0(out_dir, "umap_omp_feature.pdf"), width = 7, height = 7)
FeaturePlot(integrated_scRNA_clean, 
            features = "Omp",
            reduction = "umap",
            cols = c("lightgrey", "red"),
            pt.size = 0.5,
            order = FALSE) +
  ggtitle("Omp") +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
dev.off()
