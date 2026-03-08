cd /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/01_celltype_anno
conda activate new_r4_base
R

library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)

mix_10X <- Read10X('/md01/nieyg/project/AmbientRNA_Benchmarking/01_data/PBMC/joint_FA/joint_FA/outs/filtered_feature_bc_matrix')
pbmc <- CreateSeuratObject(counts = mix_10X$'Gene Expression', project = "mix", min.cells = 3, min.features = 200)

# add sample information
cell_id_table <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/MitoSort/Demultiplex_output_4/result_pvalue.txt')
cell_id_table <- as.data.frame(cell_id_table)
rownames(cell_id_table) <- cell_id_table$Barcode
cell_id_table_order <- cell_id_table[colnames(pbmc),]
pbmc$sample <- cell_id_table_order$Demultiplex

# delete cell 
Idents(pbmc)<-pbmc$sample 
pbmc <- subset(pbmc, idents=c('Unassign','Doublet'), invert=TRUE)

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
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2
pcs <- min(co1, co2)
pcs
# 12

pbmc <- FindNeighbors(object = pbmc, dims = 1:12)
pbmc <- FindClusters(object = pbmc, resolution = 0.8)
pbmc <- RunTSNE(object = pbmc, dims = 1:12)
pbmc <- RunUMAP(object = pbmc, dims = 1:12)


## unsupervised umap tsne
pdf("./unsupervised_umap_tsne.pdf",width=15)
DimPlot(pbmc, reduction = "tsne", group.by = c("sample", "seurat_clusters"),label=T)
DimPlot(pbmc, reduction = "umap", group.by = c("sample", "seurat_clusters"),label=T)
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


#cluster annotation
#Gene Proportion Chart
table(pbmc@active.ident)

pdf(paste0("subtype_qc_plot.pdf"),width=25)
plot1 <- VlnPlot(pbmc, features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0,)
print(plot1)
dev.off()

# make the trans dist tree 
object <- pbmc
DefaultAssay(object)<-"RNA"
embeddings <- Embeddings(object = object, reduction = "pca")[,1:12]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
library(ggtree);

pdf("ORG_cluster-tree-cosine.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()

 object = pbmc
 DefaultAssay(object)<-'RNA'
 Idents(object)<- object$RNA_snn_res.0.8
 markers <- FindAllMarkers(object,only.pos = TRUE)
 markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5) %>%
    slice_head(n = 50) %>%
    ungroup() -> top20
 pdf("seurat_cluster_heatmap_top50.pdf",width=15)
 heatmap <- DoHeatmap(object,features = top20$gene) + NoLegend()
 print(heatmap)
 dev.off()

################## enrichment ####################
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

 object = pbmc
 DefaultAssay(object)<-'RNA'
 object <- ScaleData(object, features =rownames(object),verbose = FALSE)
 markers <- FindAllMarkers(object,only.pos = TRUE,slot="scale.data")

## mouse
## Extract markers that do not contain mitochondria or ribosomes
ClusterMarker_noRibo <- markers[!grepl("^Rp[sl]", markers$gene, ignore.case=F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case=F),]

ClusterMarker_noRibo_noMito %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::filter(pct.1 > 0.1) %>%
    slice_head(n = 100) %>%
    ungroup() -> TopMarkers

Idents(object) <- object$RNA_snn_res.0.8    
pdf("seurat_cluter_heatmap_top.pdf",width=15)
 heatmap <- DoHeatmap(object,features = TopMarkers$gene) + NoLegend()
 print(heatmap)
 dev.off()

main_celltype <- 1:16
 for (i in c(7,12)) {
    current_celltype <- main_celltype[i]
    print(current_celltype)
    current_gene <- TopMarkers$gene[TopMarkers$cluster == current_celltype]
    print(current_gene)

    pdf(paste0('./enrichment_seurat_cluster/',current_celltype,"-GO-BP.pdf"),width = 10)
    gene.df <- bitr(current_gene, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Hs.eg.db)
    ego <- enrichGO(gene.df$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Hs.eg.db,
                    ont = "BP", ###BP,MF,CC
                    pAdjustMethod  = "BH",
                    pvalueCutoff = 0.3,
                    qvalueCutoff = 0.3,
                   readable = TRUE)
    ego
    barplot <- barplot(ego, showCategory=20,label_format = 80)
    write.csv(ego,paste0('./enrichment_seurat_cluster/',current_celltype,"-GO-BP.csv"))
    print(barplot)
    dev.off()

    pdf(paste0('./enrichment_seurat_cluster/',current_celltype,"-KEGG.pdf"),width = 10)
    ego <- enrichKEGG(gene = gene.df$ENTREZID,
                   keyType = "kegg",
                   organism  = 'hsa',
                   pvalueCutoff  = 0.3,
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.3)

     ego
     barplot <- barplot(ego, showCategory=20,label_format = 100)
     ego<-setReadable(ego,OrgDb=org.Hs.eg.db,keyType = "ENTREZID")
     write.csv(ego,paste0('./enrichment_seurat_cluster/',current_celltype,"-kegg.csv"))
     print(barplot)
     dev.off()
 }



################## Annotation ####################
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
 ## add new_label
 print(p&scale_x_discrete(labels=label[[i]]))
 dev.off()
}


######################################
Idents(pbmc) <- pbmc$RNA_snn_res.0.8
pbmc <- RenameIdents(
  object = pbmc,
  '0' = 'Memory_CD4_T',
  '1' = 'Naive_CD4_T',
  '2' = 'Memory_CD8_T',
  '3' = 'Mono',
  '4' = 'NK',
  '5' = 'Naive_CD8_T',
  '6' = 'NK',
  '7' = 'Memory_CD8_T',
  '8' = 'naive_B',
  '9' = 'Mono',
  '10' = 'Mono',
  '11' = 'memory_B',
  '12' = 'Memory_CD8_T',
  '13' = 'Mono',
  '14' = 'DC',
  '15' = 'plasma',
  '16' = 'DC'
  )

pbmc@meta.data$Annotation<-Idents(pbmc)
table(pbmc$Annotation,pbmc$sample)
pbmc$Annotation<-factor(pbmc$Annotation,levels=c('Mono','NK',"DC","Naive_CD4_T","Memory_CD4_T","Naive_CD8_T","Memory_CD8_T","naive_B","memory_B","plasma"))
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


all_markers <- c(
  # Progenitor cells
  "CD34", "PRSS57", "SOX4",
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
  # naive
  'FCER2', 'TCL1A', 'IL4R',
  # memory
  "CD27",'AIM2', 'TNFRSF13B',
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






