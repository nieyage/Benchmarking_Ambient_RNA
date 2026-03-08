cd /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/07_celltype_ratio
conda activate new_r4_base
R

library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

pbmc <- readRDS('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/01_celltype_anno/pbmc_human_mix_new.rds')



###################### histogram ######################
mydata <- table(pbmc$sample,pbmc$Annotation)
mydata <- t(mydata)
CD4_T <- colSums(mydata[c("Naive_CD4_T","Memory_CD4_T"), ])
B <- colSums(mydata[c("naive_B","memory_B"), ])
others <- mydata[!(rownames(mydata) %in% c("Naive_CD4_T","Memory_CD4_T","naive_B","memory_B")), ]
newdata <- rbind(others, CD4_T = CD4_T, B = B)

newdata
df_long <- melt(newdata, id.vars="celltype", 
                variable.name="Sample", 
                value.name="Count")
colnames(df_long) <- c("celltype", "Sample", "Count") 


myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

pdf('PBMC celltype ratio_2.pdf',width=10,height=4)

ggplot(df_long, aes(x=Count, y=Sample, fill=celltype)) +
  geom_bar(stat="identity",position="stack", color="black", width=0.7,size=0.25) +
  theme_bw() +
  labs(x="Cell Count", y="Sample", fill="Cell Type") + 
  scale_fill_manual(values = myUmapcolors)


dev.off()


###################### pie graph ######################
sample_result <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/MitoSort/Demultiplex_output_4/result_pvalue.txt')
sample_result <- as.data.frame(sample_result)
sample_result <- sample_result[sample_result$Demultiplex != 'Unassign',]

mydata = table(sample_result$Demultiplex)

pdf('PBMC sample ratio pie.pdf')

pie(
  mydata,
  labels = paste(names(mydata), mydata, sep = " : "),
  main = "Cell Type Distribution",
  col=c("#BFBFBF", "#499676", "#C26636", "#6F6BA6", "#C53A80")
)

dev.off()




















