cd /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/08_heatmap_mito_genome
conda activate new_r4_base
R

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(cowplot)


############################################ mitosort ############################################

specific_germline <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/MitoSort/Demultiplex_output_4/specific_germline.txt')
specific_germline <- as.data.frame(specific_germline)
colnames(specific_germline)[1] <- 'Specific_germline'
top10_by_sample <- specific_germline %>%
  group_by(Sample) %>%
  slice_head(n = 10) %>%
  ungroup()

sample_result <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/MitoSort/Demultiplex_output_4/result_pvalue.txt')
sample_result <- as.data.frame(sample_result)

SNP_matrix <- read.csv('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/MitoSort/SNP_matrix/alt.csv',row.names=1)
colnames(SNP_matrix) <- gsub('\\.','-',colnames(SNP_matrix))

filtered_cell <- sample_result$Barcode[!sample_result$Demultiplex %in% c('Doublet','Unassign')]

ref_matrix <- read.csv('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/MitoSort/SNP_matrix/ref.csv',row.names=1)
colnames(ref_matrix) <- gsub('\\.','-',colnames(ref_matrix))

all_matrix <- SNP_matrix + ref_matrix

specific_SNP_matrix <- SNP_matrix[top10_by_sample$Specific_germline,filtered_cell]
specific_all_matrix <- all_matrix[top10_by_sample$Specific_germline,filtered_cell]

table(rownames(specific_SNP_matrix) == rownames(specific_SNP_matrix))
table(colnames(specific_all_matrix) == colnames(specific_all_matrix))

data <- specific_SNP_matrix/specific_all_matrix
data[is.na(data)] <- 0
count <- data


#Column/Row Name - Sample Name Format Change
type <- factor(sample_result[match(colnames(count),sample_result$Barcode),'Demultiplex'])
annotation_col = data.frame(Sample_type = type)
rownames(annotation_col) = factor(colnames(count))

type <- factor(top10_by_sample[match(rownames(count),top10_by_sample$Specific_germline),'Sample'])
annotation_row = data.frame(Snv_type = type)
rownames(annotation_row) = factor(rownames(count))
annotation_row$Snv_type <- paste0(annotation_row$Snv_type,'-spcific')

ann_colors = list(
  Sample_type = c("Sample0"="#499676", "Sample1"="#C26636", "Sample2"="#6F6BA6", "Sample3"="#C53A80"),
  Snv_type =  c("Sample0-spcific"="#344E9B",
                "Sample1-spcific"="#481E66",
                "Sample2-spcific"="#22948F",
                "Sample3-spcific"="#F7DD77"))

sample_order <- c(sample_result$Barcode[sample_result$Demultiplex=='Sample0'],
           sample_result$Barcode[sample_result$Demultiplex=='Sample1'],
           sample_result$Barcode[sample_result$Demultiplex=='Sample2'],
           sample_result$Barcode[sample_result$Demultiplex=='Sample3'])

snv_order <- c(top10_by_sample$Specific_germline[top10_by_sample$Sample=='Sample0'],
           top10_by_sample$Specific_germline[top10_by_sample$Sample=='Sample1'],
           top10_by_sample$Specific_germline[top10_by_sample$Sample=='Sample2'],
           top10_by_sample$Specific_germline[top10_by_sample$Sample=='Sample3'])

count <- count[snv_order,sample_order]   
bk <- c(seq(0,1,by=0.01))

pdf(paste0('./',"mitosort_top10.pdf"))
p3<-pheatmap(count,cluster_cols = F,cluster_rows = F,
        color = c(colorRampPalette(colors = c("#FDF6F3","#5C191A"))(length(bk))),
        legend_breaks=seq(0,1,0.25),
        breaks=bk,
        annotation_col = annotation_col, 
        #annotation_row = annotation_row,
        annotation_colors = ann_colors,
        show_rownames=T,show_colnames=F)
print(p3)
dev.off()

