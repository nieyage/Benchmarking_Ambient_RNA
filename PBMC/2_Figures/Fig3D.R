cd /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/09_contamination_gene_ratio
conda activate new_r4_base
R

library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
library(rhdf5)

current_gene <- 'MS4A1'
#################### Rawdata ####################
ifnb <- readRDS('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/01_celltype_anno/pbmc_human_mix_new.rds') 
rawdata_vcan <- ifnb@assays$RNA$counts[current_gene,]
rawdata_numi <- ifnb$nCount_RNA

table(names(rawdata_vcan) == names(rawdata_numi))

B_cell <- colnames(ifnb)[ifnb$Annotation %in% c('naive_B','memory_B')]
non_B_cell <- colnames(ifnb)[! ifnb$Annotation %in% c('naive_B','memory_B')]
table(rawdata_vcan[B_cell] > 0)[2] / length(rawdata_vcan[B_cell])
# 0.8198433
table(rawdata_vcan[non_B_cell] > 0)[2] / length(rawdata_vcan[non_B_cell])
# 0.01472168 

umi_ratio <- rawdata_vcan/rawdata_numi
names(umi_ratio) <- names(rawdata_vcan) 
mean(umi_ratio[B_cell])
# 0.001577039
mean(umi_ratio[non_B_cell])
# 1.25706e-05

save(B_cell, non_B_cell, file = "B.RData")
#################### cellbender ####################

data.file <- '/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/02_cellbender/pbmc_human_mix_output_filtered_seurat.h5'
data.data <- Read10X_h5(filename = data.file, use.names = TRUE)

# create Seurat object
obj <- CreateSeuratObject(counts = data.data$'Gene Expression')
obj

PBMC_10k_pbmc <- readRDS('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/01_celltype_anno/pbmc_human_mix.rds')
target_cell <- colnames(obj)[colnames(obj) %in% colnames(PBMC_10k_pbmc)]
obj_sub <- subset(obj, cells = target_cell)
obj_sub

rawdata_vcan <- obj_sub@assays$RNA$counts[current_gene,]
rawdata_numi <- obj_sub$nCount_RNA


current_B_cell <- colnames(obj_sub)[colnames(obj_sub) %in% B_cell]
current_non_B_cell <- setdiff(colnames(obj_sub),current_B_cell)
length(current_B_cell)
length(current_non_B_cell)

table(rawdata_vcan[current_B_cell] > 0)[2] / length(rawdata_vcan[current_B_cell])
# 0.8211169
table(rawdata_vcan[current_non_B_cell] > 0)[2] / length(rawdata_vcan[current_non_B_cell])
# 0.00126968 

umi_ratio <- rawdata_vcan/rawdata_numi
names(umi_ratio) <- names(rawdata_vcan) 
mean(umi_ratio[current_B_cell])
# 0.001922826
mean(umi_ratio[current_non_B_cell])
# 7.295317e-07

#################### cellclear ####################
ifnb_cellclear <- Read10X('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/09_contamination_gene_ratio/cellclear_matrix/')
ifnb_cellclear <- CreateSeuratObject(counts = ifnb_cellclear)

target_cell <- colnames(ifnb_cellclear)[colnames(ifnb_cellclear) %in% colnames(PBMC_10k_pbmc)] #16935 
obj_sub <- subset(ifnb_cellclear, cells = target_cell)
obj_sub

rawdata_vcan <- obj_sub@assays$RNA$counts[current_gene,]
rawdata_numi <- obj_sub$nCount_RNA


current_B_cell <- colnames(obj_sub)[colnames(obj_sub) %in% B_cell]
current_non_B_cell <- setdiff(colnames(obj_sub),current_B_cell)
length(current_B_cell)
length(current_non_B_cell)

table(rawdata_vcan[current_B_cell] > 0)[2] / length(rawdata_vcan[current_B_cell])
# 0.8198433
table(rawdata_vcan[current_non_B_cell] > 0)[2] / length(rawdata_vcan[current_non_B_cell])
# 0.01456987

umi_ratio <- rawdata_vcan/rawdata_numi
names(umi_ratio) <- names(rawdata_vcan) 
mean(umi_ratio[current_B_cell])
# 0.00227507
mean(umi_ratio[current_non_B_cell])
# 2.020072e-05

#################### SoupX ####################
ifnb_soupx <- readRDS('/data/R04/yanjw7/workspace/250801_soupX/joint_cell/with_FA_human_mix/pbmc_soupx_cleaned.rds')

target_cell <- colnames(ifnb_soupx)[colnames(ifnb_soupx) %in% colnames(PBMC_10k_pbmc)] #16976 
obj_sub <- subset(ifnb_soupx, cells = target_cell)
obj_sub

rawdata_vcan <- obj_sub@assays$RNA$counts[current_gene,]
rawdata_numi <- obj_sub$nCount_RNA

current_B_cell <- colnames(obj_sub)[colnames(obj_sub) %in% B_cell]
current_non_B_cell <- setdiff(colnames(obj_sub),current_B_cell)
length(current_B_cell)
length(current_non_B_cell)

table(rawdata_vcan[current_B_cell] > 0)[2] / length(rawdata_vcan[current_B_cell])
# 0.8198433
table(rawdata_vcan[current_non_B_cell] > 0)[2] / length(rawdata_vcan[current_non_B_cell])
# 0.006697416

umi_ratio <- rawdata_vcan/rawdata_numi
names(umi_ratio) <- names(rawdata_vcan) 
mean(umi_ratio[current_B_cell])
# 0.001566769
mean(umi_ratio[current_non_B_cell])
# 2.573656e-06


#################### Fastcar ####################
ifnb_fastcar <- readRDS('/data/R03/luoxy237/workspace/2025_summer/ambientRNA_pbmc_contamination/FastCAR/07_human_mix/human_mix.rds')
ifnb_fastcar <- CreateSeuratObject(counts = ifnb_fastcar)

target_cell <- colnames(ifnb_fastcar)[colnames(ifnb_fastcar) %in% colnames(PBMC_10k_pbmc)] #16976 
obj_sub <- subset(ifnb_fastcar, cells = target_cell)
obj_sub

rawdata_vcan <- obj_sub@assays$RNA$counts[current_gene,]
rawdata_numi <- obj_sub$nCount_RNA

current_B_cell <- colnames(obj_sub)[colnames(obj_sub) %in% B_cell]
current_non_B_cell <- setdiff(colnames(obj_sub),current_B_cell)
length(current_B_cell)
length(current_non_B_cell)

table(rawdata_vcan[current_B_cell] > 0)[2] / length(rawdata_vcan[current_B_cell])
# 0.8198433
table(rawdata_vcan[current_non_B_cell] > 0)[2] / length(rawdata_vcan[current_non_B_cell])
# 0.01472168

umi_ratio <- rawdata_vcan/rawdata_numi
names(umi_ratio) <- names(rawdata_vcan) 
mean(umi_ratio[current_B_cell])
# 0.001730944
mean(umi_ratio[current_non_B_cell])
# 1.375158e-05


#################### DecontX ####################
ifnb_decontx <- readRDS('/data/R04/wuchx37/data_wcx/pbmc_cellranger/tools_decontX/human_mix_cleaned_counts.rds')
ifnb_decontx <- CreateSeuratObject(round(ifnb_decontx))

target_cell <- colnames(ifnb_decontx)[colnames(ifnb_decontx) %in% colnames(PBMC_10k_pbmc)] #16976 
obj_sub <- subset(ifnb_decontx, cells = target_cell)
obj_sub

rawdata_vcan <- obj_sub@assays$RNA$counts[current_gene,]
rawdata_numi <- obj_sub$nCount_RNA

current_B_cell <- colnames(obj_sub)[colnames(obj_sub) %in% B_cell]
current_non_B_cell <- setdiff(colnames(obj_sub),current_B_cell)
length(current_B_cell)
length(current_non_B_cell)

table(rawdata_vcan[current_B_cell] > 0)[2] / length(rawdata_vcan[current_B_cell])
# 0.8154917
table(rawdata_vcan[current_non_B_cell] > 0)[2] / length(rawdata_vcan[current_non_B_cell])
# 0.001453213

umi_ratio <- rawdata_vcan/rawdata_numi
names(umi_ratio) <- names(rawdata_vcan) 
mean(umi_ratio[current_B_cell])
# 0.001587327
mean(umi_ratio[current_non_B_cell])
# 9.900056e-07


#################### siftcell ####################
ifnb_siftcell <- readRDS('/data/R04/wuchx37/data_wcx/pbmc_cellranger/tools_siftcell/Human_mix/filter_humanmix.rds')
ifnb_siftcell <- CreateSeuratObject(round(ifnb_siftcell@assays$RNA$counts), min.features = 200)


target_cell <- colnames(ifnb_siftcell)[colnames(ifnb_siftcell) %in% colnames(PBMC_10k_pbmc)] #16762 
obj_sub <- subset(ifnb_siftcell, cells = target_cell)
obj_sub

rawdata_vcan <- obj_sub@assays$RNA$counts[current_gene,]
rawdata_numi <- obj_sub$nCount_RNA

current_B_cell <- colnames(obj_sub)[colnames(obj_sub) %in% B_cell]
current_non_B_cell <- setdiff(colnames(obj_sub),current_B_cell)
length(current_B_cell)
length(current_non_B_cell)

table(rawdata_vcan[current_B_cell] > 0)[2] / length(rawdata_vcan[current_B_cell])
# 0.8236332
table(rawdata_vcan[current_non_B_cell] > 0)[2] / length(rawdata_vcan[current_non_B_cell])
# 0.0145892

umi_ratio <- rawdata_vcan/rawdata_numi
names(umi_ratio) <- names(rawdata_vcan) 
mean(umi_ratio[current_B_cell])
# 0.001590304
mean(umi_ratio[current_non_B_cell])
# 1.243533e-05


##################### scAR ####################
ifnb_scAR <- readRDS('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/humanmix_scAR_seurat.rds')
target_cell <- colnames(ifnb_scAR)[colnames(ifnb_scAR) %in% colnames(PBMC_10k_pbmc)] #16976
obj_sub <- subset(ifnb_scAR, cells = target_cell)
obj_sub

rawdata_vcan <- obj_sub@assays$RNA$counts['MS4A1',]
rawdata_numi <- obj_sub$nCount_RNA

current_B_cell <- colnames(obj_sub)[colnames(obj_sub) %in% B_cell]
current_non_B_cell <- setdiff(colnames(obj_sub),current_B_cell)
length(current_B_cell)
length(current_non_B_cell)


table(rawdata_vcan[current_B_cell] > 0)[2] / length(rawdata_vcan[current_B_cell])
# 0.9895561
table(rawdata_vcan[current_non_B_cell] > 0)[2] / length(rawdata_vcan[current_non_B_cell])
# 0.0008213812

umi_ratio <- rawdata_vcan/rawdata_numi
names(umi_ratio) <- names(rawdata_vcan) 
mean(umi_ratio[current_B_cell])
# 0.002205346
mean(umi_ratio[current_non_B_cell])
# 8.392954e-07



#################### scCDC ####################
cd /public/home/chenxy/project/Mrs.nie/contain/12_PBMC_human_mix/soupx
conda activate r4-base 
R

library(Seurat)
library(scCDC)

load('B.RData')
current_gene <- 'MS4A1'


seuratobject <- readRDS("./pbmc_human_mix_new.rds")
# detect global contamination causing genes(GCGs)
GCGs <- ContaminationDetection(seuratobject)           
head(GCGs)

# Contamination quantification 
if ("Annotation" %in% colnames(seuratobject@meta.data)) {annotations <- seuratobject@meta.data$Annotation
Idents(seuratobject) <- annotations} else {stop("no column named 'Annotation' in the Seurat object")}
new_ident <- gsub("_", "-", Idents(seuratobject))
Idents(seuratobject) <- new_ident

mislet_cont_ratio <- ContaminationQuantification(seuratobject,rownames(GCGs))
mislet_cont_ratio

# remove the contamination
seuratobj_corrected <- ContaminationCorrection(seuratobject,rownames(GCGs))

# extract the stain removal count matrix
corrected_count_matrix = data.frame(seuratobj_corrected@assays$Corrected$counts)
head(corrected_count_matrix)

ifnb_scCDC <- CreateSeuratObject(corrected_count_matrix) #16976

colnames(ifnb_scCDC) <- gsub("\\.", "-", colnames(ifnb_scCDC))
head(colnames(ifnb_scCDC))

rawdata_vcan <- ifnb_scCDC@assays$RNA$counts[current_gene,]
rawdata_numi <- ifnb_scCDC$nCount_RNA

current_B_cell <- colnames(ifnb_scCDC)[colnames(ifnb_scCDC) %in% B_cell]
current_non_B_cell <- setdiff(colnames(ifnb_scCDC),current_B_cell)
length(current_B_cell)
length(current_non_B_cell)

table(rawdata_vcan[current_B_cell] > 0)[2] / length(rawdata_vcan[current_B_cell])
# 0.8198433
table(rawdata_vcan[current_non_B_cell] > 0)[2] / length(rawdata_vcan[current_non_B_cell])
# 0.01472168

umi_ratio <- rawdata_vcan/rawdata_numi
names(umi_ratio) <- names(rawdata_vcan) 
mean(umi_ratio[current_B_cell])
# 0.001750829
mean(umi_ratio[current_non_B_cell])
# 1.398658e-05


######################################################################################################################

cd /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/09_contamination_gene_ratio/MS4A1
conda activate new_r4_base
R

library(ggplot2)
library(dplyr)
library(scales)
library(ggbreak)

cell_expression_ratio <-  read.csv('cell_expression_ratio.csv',header=T)
cell_expression_ratio$cell_ratio <- as.numeric(cell_expression_ratio$cell_ratio)
cell_expression_ratio$tool <- factor(cell_expression_ratio$tool,levels=c('Rawdata','DecontX','SoupX','scAR','Cellbender','Fastcar','scCDC','Cellclear'))


myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

pdf('cell_expression_ratio.pdf',width=5,height=4)

p <- ggplot(cell_expression_ratio, aes(x = tool, y = cell_ratio, fill = cell)) +
  geom_col(position = "dodge", width = 0.7) +  
  labs(title = "Cell Expression Ratio by Tool",
       x = "Tools",
       y = "Ratio",
       fill = "Cell Type") +
  scale_fill_manual(values = c("B_cell" = "#AB3282", "Others" = "#BFBFBF")) +  
  scale_y_continuous(labels = percent_format(accuracy = 1)) +   
  theme_classic() +                                             
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  scale_y_break(c(0.03, 0.6), scales = 1)    

print(p)
dev.off()




cell_umi_ratio <-  read.csv('cell_umi_ratio.csv',header=T)
cell_umi_ratio$UMI_ratio <- as.numeric(cell_umi_ratio$UMI_ratio)
cell_umi_ratio$tool <- factor(cell_umi_ratio$tool,levels=c('Rawdata','DecontX','SoupX','scAR','Cellbender','Fastcar','scCDC','Cellclear'))
cell_umi_ratio

pdf('cell_umi_ratio.pdf',width=5,height=4)

p <- ggplot(cell_umi_ratio, aes(x = tool, y = UMI_ratio, fill = cell)) +
  geom_col(position = "dodge", width = 0.7) + 
  labs(title = "Cell Expression Ratio by Tool",
       x = "Tools",
       y = "Ratio",
       fill = "Cell Type") +
  scale_fill_manual(values = c("B_cell" = "#AB3282", "Others" = "#BFBFBF")) +  
  scale_y_continuous(labels = percent_format(accuracy = 0.0001)) +   
  theme_classic() +                                             
  theme(axis.text.x = element_text(angle = 45, hjust = 1))   
  scale_y_break(c(0.00005, 0.001), scales = 1)  

print(p)
dev.off()


