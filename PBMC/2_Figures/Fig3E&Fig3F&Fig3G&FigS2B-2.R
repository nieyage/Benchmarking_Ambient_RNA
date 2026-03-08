################################ 01 Separate RNA data ################################
cd /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/01
conda activate new_r4_base
R

library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
library(DropletUtils)


mix_10X <- Read10X('/md01/nieyg/project/AmbientRNA_Benchmarking/01_data/PBMC/joint_FA/joint_FA/outs/filtered_feature_bc_matrix')
pbmc <- CreateSeuratObject(counts = mix_10X$'Gene Expression', project = "mix")
DropletUtils:::write10xCounts("human_mix_filtered", 
                              pbmc@assays$RNA$counts,
                              barcodes = colnames(pbmc))


mix_10X_raw <- Read10X('/md01/nieyg/project/AmbientRNA_Benchmarking/01_data/PBMC/joint_FA/joint_FA/outs/raw_feature_bc_matrix')
pbmc <- CreateSeuratObject(counts = mix_10X_raw$'Gene Expression', project = "mix")
DropletUtils:::write10xCounts("human_mix_raw", 
                              pbmc@assays$RNA$counts,
                              barcodes = colnames(pbmc))







################################ 02 scAR ################################

cd /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/02
conda activate scAR
# conda create -n scAR python=3.10
# conda activate scAR
# pip install matplotlib==3.1.3  # Specify this matplotlib version to avoid errors
# pip install scanpy
# pip install git+https://github.com/Novartis/scAR.git


python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scar import model, setup_anndata

import warnings
warnings.simplefilter("ignore")



### data ###

human_mix_raw = sc.read_10x_mtx('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/01/human_mix_raw')
human_mix_raw.var_names_make_unique()

human_mix = sc.read_10x_mtx('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/01/human_mix_filtered')
human_mix.var_names_make_unique()

human_mix_raw.var['feature_types'] = 'Gene Expression'
human_mix.var['feature_types'] = 'Gene Expression'

### Estimated environmental overview ###
import matplotlib
matplotlib.use("Agg")   

import matplotlib.pyplot as plt
from scar import setup_anndata

setup_anndata(
    adata=human_mix,
    raw_adata=human_mix_raw,
    prob=0.995,
    kneeplot=True
)

plt.savefig("kneeplot.png", dpi=300, bbox_inches="tight")
plt.close()

human_mix.uns["ambient_profile_Gene Expression"].head()



### train ###
human_mix_scar = model(raw_count=human_mix, # In the case of Anndata object, scar will automatically use the estimated ambient_profile present in adata.uns.
                      ambient_profile=human_mix.uns['ambient_profile_Gene Expression'],
                      feature_type='mRNA',
                      sparsity=1,
                      device='cuda' # CPU, CUDA and MPS are supported.
                     )

human_mix_scar.train(epochs=200,
                    batch_size=64,
                    verbose=True
                   )

human_mix_scar.inference()
denoised_count = pd.DataFrame(human_mix_scar.native_counts.toarray(), index=human_mix.obs_names, columns=human_mix.var_names)
denoised_count.head()


denoised_count.to_csv("denoised_count.csv")






################################ 03 seurat ################################
cd /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/05_combined_vcf
conda activate new_r4_base
R

library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
library(rhdf5)

ifnb_scAR <- read.csv("/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/01/denoised_count.csv",row.names=1)
ifnb_scAR_seurat = CreateSeuratObject(counts =  t(ifnb_scAR))

saveRDS(ifnb_scAR_seurat,'/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/humanmix_scAR_seurat.rds')
















cd /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/03_anno
conda activate new_r4_base
R

library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
library(rhdf5)


pbmc <- readRDS('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/humanmix_scAR_seurat.rds')


# add sample information
cell_id_table <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/MitoSort/Demultiplex_output_4/result_pvalue.txt')
cell_id_table <- as.data.frame(cell_id_table)
rownames(cell_id_table) <- cell_id_table$Barcode
cell_id_table <- cell_id_table[cell_id_table$Demultiplex %in% c('Sample0','Sample1','Sample2','Sample3'),]

selected_cell <- intersect(rownames(cell_id_table),colnames(pbmc))
pbmc$sub <- ifelse(colnames(pbmc) %in% selected_cell,'yes','no')
Idents(pbmc) <- pbmc$sub
pbmc <- subset(pbmc, idents=c('yes'), invert=FALSE)
pbmc

cell_id_table_order <- cell_id_table[colnames(pbmc),]
table(rownames(cell_id_table_order)==colnames(pbmc))
pbmc$sample <- cell_id_table_order$Demultiplex
table(pbmc$sample)

# delete cell 
Idents(pbmc)<-colnames(pbmc) 

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 20000 & percent.mt < 40)
pbmc <- NormalizeData(pbmc,verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, verbose = FALSE)
pbmc <- ScaleData(pbmc, features =rownames(pbmc),verbose = FALSE)
pbmc <- RunPCA(pbmc, npcs = 50, verbose = FALSE)

# PCA selection #
pct <- pbmc[["pca"]]@stdev / sum(pbmc[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1
co2 <- sort(
  which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), 
  decreasing = TRUE
)[1] + 1
co2
pcs <- min(co1, co2)
pcs
# 13

pbmc <- FindNeighbors(object = pbmc, dims = 1:13)
pbmc <- FindClusters(object = pbmc, resolution = 1.5)
pbmc <- RunTSNE(object = pbmc, dims = 1:13)
pbmc <- RunUMAP(object = pbmc, , dims = 1:13,spread=3)

pdf("./unsupervised_umap_tsne.pdf",width=15)
DimPlot(pbmc, reduction = "tsne", group.by = c("sample","seurat_clusters"),label=T)
DimPlot(pbmc, reduction = "umap", group.by = c("sample","seurat_clusters"),label=T)
dev.off()

pdf(paste0("Cluster_qc_plot.pdf"),width=15)
plot1 <- VlnPlot(pbmc, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0,)
print(plot1)
dev.off()

################## Anotation ####################
monocytes= c("LYZ",'CD14','CSF3R','S100A8','S100A9',"VCAN","ITGAM","ITGAX","CD68","S100A12","MTSS1","SERPINA1","MS4A7","SIGLEC10","APLP2","MPEG1",'FCGR3A','CDKN1C','CX3CR1')                   
DC= c("FLT3","FCER1A","CLEC10A","CLEC9A","CD1C","AXL","SIGLEC6",'LILRA4','MZB1',"IRF8","IRF7","SPI8","DERL3","SPI1","PTPRS")
T= c("CD3D", "CD3E", "CD3G") 
NK= c("NKG7","GNLY","FGFBP2","SPON2","NCAM1","KLRC1","XCL2","KLRF1","GZMB","FCGR3A")                      
B = c("MS4A1","CD27","IGHD","IGHM","TCL1A","FCRL5")                      
plasma = c("JCHAIN","TNFRSF17","ITM2C", 'MZB1', 'IGKC')                  
platelet = c("ITGA2B","PF4","TUBB1","PPBP","GP1BB","NRGN")
Progenitor = c("CD34","PRSS57","SOX4")
CD8T=c("CD8A","CD8B")
CD4T=c("CD4")
B_naive=c("IGHD", 'FCER2', 'TCL1A', 'IL4R')
B_memory=c("CD27",'AIM2', 'TNFRSF13B')
B_detail=c("CD20","CD138","CD27","CD38","CD79A")
B_germinal_center=c('S1PI2', 'LRMP', 'SUGCT', 'MME', 'MKI67','AICDA')
Naive_CD4T=c('CCR7','SELL','TCF7','IL7R','GPR183','LEF1','LTB')
Naive_CD8T=c('CCR7','SELL','S100A8','CST3','AC020916.1','IL7R','LEF1','TCF7')
Naive_T=c('CCR7','SELL','TCF7','IL7R','FOXO1','LEF1','LTB','GZMK')
Naive_CD4_T_2 = c("SELL", "CCR7", "CD95")
Naive_CD8_T = c("LEF1", "CD8A", "IL2RB")
Memory_CD4_T_2 = c("IL7R", "CD52", "CD4")
Memory_CD8_T = c("IFNG", "EMOES")  

markerGenes  <- list(monocytes,DC,T,NK,B,plasma,platelet,Progenitor,CD8T,CD4T,B_naive,B_memory,B_detail,B_germinal_center,Naive_CD4T,
                     Naive_CD8T,Naive_T,Naive_CD4_T_2,Naive_CD8_T,Memory_CD4_T_2,Memory_CD8_T)

markerGenes_mouse <- markerGenes 
label<- list(rep("monocytes",length(monocytes)),
          rep("DC",length(DC)),
          rep("T",length(T)),
          rep("NK",length(NK)),
          rep("B",length(B)),
          rep("plasma",length(plasma)),
          rep("platelet",length(platelet)),
          rep("Progenitor",length(Progenitor)),
          rep("CD8T",length(CD8T)),
          rep("CD4T",length(CD4T)),
          rep("B_naive",length(B_naive)),
          rep("B_memory",length(B_memory)),
          rep("B_detail",length(B_detail)),
         rep("B_germinal_center",length(B_germinal_center)),
         rep("Naive_CD4T",length(Naive_CD4T)),
         rep("Naive_CD8T",length(Naive_CD8T)),
         rep("Naive_T",length(Naive_T)),
         rep("Naive_CD4_T_2",length(Naive_CD4_T_2)),
         rep("Naive_CD8_T",length(Naive_CD8_T)),
         rep("Memory_CD4_T_2",length(Memory_CD4_T_2)),
         rep("Memory_CD8_T",length(Memory_CD8_T))
         )
celltype <- c('monocytes','DC','T','NK','B','plasma','platelet','Progenitor','CD8T','CD4T','B_naive','B_memory',"B_detail","B_germinal_center",'Naive_CD4T','Naive_CD8T',
              'Naive_T','Naive_CD4_T_2','Naive_CD8_T','Memory_CD4_T_2','Memory_CD8_T')

DefaultAssay(pbmc)<-"RNA"
for (i in 1:length(celltype)) {
 print(celltype[i])
 pdf(paste0("./PBMC_mix_cluster-annotation-all_celltype-basic/",celltype[i],".pdf"))
 p<-DotPlot(pbmc, features = markerGenes_mouse[[i]],dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
 print(p)
 print(p&scale_x_discrete(labels=label[[i]]))
 dev.off()
}

##########

##########
######################################

Idents(pbmc) <- pbmc$RNA_snn_res.1.5

pbmc <- RenameIdents(
  object = pbmc,
  '0' = 'Memory_CD8_T',
  '1' = 'Naive_CD8_T',
  '2' = 'NK',
  '3' = 'NK',
  '4' = 'Mono',
  '5' = 'B_cell',
  '6' = 'CD4_T',
  '7' = 'CD4_T',
  '8' = 'Memory_CD8_T',
  '9' = 'CD4_T',
  '10' = 'Mono',
  '11' = 'Memory_CD8_T',
  '12' = 'CD4_T',
  '13' = 'Mono',
  '14' = 'Naive_CD8_T',
  '15' = 'CD4_T',
  '16' = 'Mono',
  '17' = 'DC',
  '18' = 'Mono',
  '19' = 'Memory_CD8_T',
  '20' = 'B_cell',
  '21' = 'DC',
  '22' = 'plasma',
  '23' = 'Memory_CD8_T'
  )

pbmc@meta.data$Annotation<-Idents(pbmc)
table(pbmc$Annotation,pbmc$sample)
pbmc$Annotation<-factor(pbmc$Annotation,levels=c('Mono','NK',"DC","CD4_T","Naive_CD8_T","Memory_CD8_T","B_cell","plasma"))
Idents(pbmc)<-pbmc$Annotation

myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

pdf("./Annotated_maincelltype_UMAP.pdf",width=6,height=5)
DimPlot(pbmc, repel = TRUE, cols=myUmapcolors, reduction = "umap",group.by = "Annotation",label = TRUE)
DimPlot(pbmc, reduction = "umap", group.by = "sample",label = TRUE)
dev.off()



top_gene <- c('VCAN','MS4A1','RNF175','YPEL5','SIGLEC10','KLRF1','SELL','CD3E','CD3D','SPI1')
pdf('./top_gene_FeaturePlot_all.pdf')
for (i in 1:length(top_gene)) {
  print(i)
  print(top_gene[i])
  p1 <- FeaturePlot(pbmc, features = top_gene[i],
                    cols = c("lightgrey", "#FF0A07"))
    print(p1)
}
dev.off()


pdf("./sample_split_umap.pdf",width=25)
DimPlot(pbmc,
        reduction = "umap", 
        group.by = c("sample"),
        label=F,
        split.by='sample',
        label.size = 10,
        cols=c("#499676", "#C26636", "#6F6BA6", "#C53A80"))
dev.off()

saveRDS(pbmc,'scAR_humanmix.rds')



##########

all_markers <- c(
  # Monocytes
  "LYZ", "CD14", "CSF3R", "S100A8", "S100A9", "VCAN", "ITGAM", "ITGAX",
  "CD68", "S100A12", "MTSS1", "SERPINA1", "MS4A7", "SIGLEC10", "APLP2", "MPEG1",
  "CDKN1C", "CX3CR1",
  # NK cells
  "NKG7", "GNLY", "FGFBP2", "SPON2", "NCAM1", "KLRC1", "XCL2", "KLRF1", "GZMB", "FCGR3A",
  # Dendritic Cells (DC)
  "FLT3", "FCER1A", "CLEC10A", "CLEC9A", "CD1C", "AXL", "SIGLEC6", "LILRA4",
  "IRF8", "IRF7", "SPI8", "DERL3", "SPI1", "PTPRS",
  # T cells
  "CD3D", "CD3E", "CD3G",
  # CD4+ T cells
  "CD4",
  # CD8+ T cells
  "CD8A", "CD8B",
  # Naive T
  'CCR7','SELL','TCF7','FOXO1','LEF1',
  # Memory T
  'IFNG','CD52','IL7R',
  # B cells
  "MS4A1", "IGHD", "IGHM",  "FCRL5",
  # plasma
  "JCHAIN","TNFRSF17","ITM2C", 'MZB1', 'IGKC'
)

marker_genes <- all_markers

pdf("cluster-annotation-all_celltype_5.pdf",width=10,height=8)
p<-DotPlot(pbmc, features = marker_genes, dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
dev.off()

saveRDS(pbmc,'pbmc_human_mix_new.rds')
write.csv(markers,'./all_markers_new.csv')





#####################################################################################################################
#####################################################################################################################
pbmc_raw <- readRDS('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/01_celltype_anno/pbmc_human_mix_new.rds')
pbmc_raw$Annotation <- as.character(pbmc_raw$Annotation)
pbmc_raw$Annotation_merged <- pbmc_raw$Annotation
pbmc_raw$Annotation_merged[pbmc_raw$Annotation_merged %in% 
                             c("Naive_CD4_T", "Memory_CD4_T")] <- "CD4_T"
pbmc_raw$Annotation_merged[pbmc_raw$Annotation_merged %in% 
                             c("naive_B", "memory_B")] <- "B_cell"
table(pbmc_raw$Annotation_merged)
pbmc_raw$Annotation <- pbmc_raw$Annotation_merged

pdf("./pbmc_raw/pbmc_raw_Annotated_maincelltype_UMAP.pdf",width=6,height=5)
DimPlot(pbmc_raw, repel = TRUE, cols=myUmapcolors, reduction = "umap",group.by = "Annotation",label = TRUE)
DimPlot(pbmc_raw, reduction = "umap", group.by = "sample",label = TRUE)
dev.off()


seurat_list <- list(pbmc,pbmc_raw)
seurat_list_name <- c('scAR','Rawdata')

all_expression_ratio <- data.frame()
all_umi_ratio <- data.frame()

for (i in 1:length(top_gene)) {
 current_gene <- top_gene[i]
 print(current_gene)
 for (j in 1:length(seurat_list_name)) {
   current_seurat <- seurat_list[[j]]
   print(seurat_list_name[j])
   current_seurat_name <- seurat_list_name[j]
   rawdata_vcan <- current_seurat@assays$RNA$counts[current_gene,]
   rawdata_numi <- current_seurat$nCount_RNA
   for (m in 1:length(unique(pbmc$Annotation))) {
     current_celltype <- unique(pbmc$Annotation)[m]
     print(current_celltype)
     current_cell <- colnames(current_seurat)[current_seurat$Annotation==current_celltype]

     expression_ratio <- table(rawdata_vcan[current_cell] > 0)[2] / length(rawdata_vcan[current_cell])
     expression_ratio <- ifelse(is.na(expression_ratio),0,expression_ratio)
     names(expression_ratio) <- NULL
     print(expression_ratio)
     
     umi_ratio <- rawdata_vcan/rawdata_numi
     names(umi_ratio) <- names(rawdata_vcan) 
     umi_ratio <- mean(umi_ratio[current_cell])
     umi_ratio <- ifelse(is.na(umi_ratio), 0, umi_ratio)
     names(umi_ratio) <- NULL
     print(umi_ratio)



     current_df_expression_ratio<- data.frame(gene=current_gene,
                                              tool=current_seurat_name,
                                              cell=current_celltype,
                                              expression_ratio=expression_ratio)
     current_df_umi_ratio<- data.frame(gene=current_gene,
                                              tool=current_seurat_name,
                                              cell=current_celltype,
                                              umi_ratio=umi_ratio)
     all_expression_ratio <- rbind(all_expression_ratio,current_df_expression_ratio)
     all_umi_ratio <- rbind(all_umi_ratio,current_df_umi_ratio)
    
   }
 }
}

write.csv(all_expression_ratio,'all_expression_ratio.csv')
write.csv(all_umi_ratio,'all_umi_ratio.csv')











#####################################################################################################################
#####################################################################################################################
cd /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/05
conda activate new_r4_base
R

library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)

ifnb <- readRDS('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/03_anno/scAR_humanmix.rds')
ifnb_old <- readRDS('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/01_celltype_anno/pbmc_human_mix_new.rds')
ifnb_old$Annotation <- as.character(ifnb_old$Annotation)
ifnb_old$Annotation_merged <- ifnb_old$Annotation
ifnb_old$Annotation_merged[ifnb_old$Annotation_merged %in% 
                             c("Naive_CD4_T", "Memory_CD4_T")] <- "CD4_T"
ifnb_old$Annotation_merged[ifnb_old$Annotation_merged %in% 
                             c("naive_B", "memory_B")] <- "B_cell"
table(ifnb_old$Annotation_merged)
ifnb_old$Annotation <- ifnb_old$Annotation_merged


old_meta <- ifnb_old@meta.data
ifnb$old_anno <- old_meta[colnames(ifnb),]$Annotation
ifnb$Annotation_diff <- ifelse(ifnb$Annotation == ifnb$old_anno,'same','different')


umap_info <- ifnb@reductions$umap@cell.embeddings
umap_info <- as.data.frame(umap_info)

data <- ifnb@meta.data
data <- data[,c('Annotation','Annotation_diff')]
data$umap1 <- umap_info$umap_1
data$umap2 <- umap_info$umap_2
head(data)

myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

pdf('umap_diff_anno.pdf')

# UMAP 
p <- ggplot(data, aes(x = umap1, y = umap2, color = Annotation)) +
  geom_point(shape = 16, size = 0.5) +  
  geom_point(
    data = subset(data, Annotation_diff == "different"),
    shape = 21, size = 0.7, stroke = 0.1, color = "black"
  ) +
  scale_color_manual(values = c(myUmapcolors)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_classic() +
  labs(x = "UMAP1", y = "UMAP2", color = "Annotation")

print(p)

dev.off()


#################################################################################################################################################

cd /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/04
conda activate new_r4_base
R

library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)
library(dplyr)
library(ggbreak)


cell_expression_ratio <-  read.csv('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/03_anno/all_expression_ratio.csv',header=T, row.names=1)
cell_expression_ratio <- cell_expression_ratio[cell_expression_ratio$gene == 'MS4A1',]


cell_expression_ratio$cell<-factor(cell_expression_ratio$cell,levels=c('Mono','NK',"DC","CD4_T","Naive_CD8_T","Memory_CD8_T","B_cell","plasma"))
cell_expression_ratio$tool <- factor(cell_expression_ratio$tool,levels=c('Rawdata','scAR'))

pdf('cell_expression_ratio.pdf', width = 5, height = 4)

p <- ggplot(cell_expression_ratio, aes(x = cell, y = expression_ratio, fill = tool)) +
  geom_col(position = "dodge", width = 0.7) +
  labs(title = "Cell Expression Ratio by Tool",
       x = "Tools",
       y = "Percentage",
       fill = "Cell Type") +
  scale_fill_manual(values = c("Rawdata" = "#BFBFBF", "scAR" = "#0052A1")) +
  scale_y_continuous(labels = percent_format(accuracy = 0.01)) +   
  theme_classic() +                                             
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 # scale_y_break(c(0.1, 0.3), scales = 1)  


print(p)
dev.off()




#################################################################################################################################################
cell_umi_ratio <-  read.csv('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/03_anno/all_umi_ratio.csv',header=T, row.names=1)
cell_umi_ratio <- cell_umi_ratio[cell_umi_ratio$gene == 'MS4A1',]

cell_umi_ratio$cell<-factor(cell_umi_ratio$cell,levels=c('Mono','NK',"DC","CD4_T","Naive_CD8_T","Memory_CD8_T","B_cell","plasma"))
cell_umi_ratio$tool <- factor(cell_umi_ratio$tool,levels=c('Rawdata','scAR'))

pdf('cell_umi_ratio.pdf', width = 5, height = 4)


p <- ggplot(cell_umi_ratio, aes(x = cell, y = umi_ratio, fill = tool)) +
  geom_col(position = "dodge", width = 0.7) +
  labs(title = "Cell umi Ratio by Tool",
       x = "Tools",
       y = "Percentage",
       fill = "Cell Type") +
  scale_fill_manual(values = c("Rawdata" = "#BFBFBF", "scAR" = "#0052A1")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.0001)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
 # scale_y_break(c(0.000050, 0.0014), scales = 1)  


print(p)

dev.off()





