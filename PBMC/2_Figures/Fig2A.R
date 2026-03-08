cd /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA
# clone MitoSort repository
git clone https://github.com/tangzhj/MitoSort.git

# create a conda environment and install all the required python packages
#demultiplex sample pipeline of mitoSort
#resource : https://github.com/tangzhj/MitoSort

conda activate MitoSort
#data 
lib1 path /data/R03/zhangwx/project/human_PBMC/lib1
bam_path="/md01/nieyg/project/AmbientRNA_Benchmarking/01_data/PBMC/joint_FA/joint_FA/outs/atac_possorted_bam.bam"

#Realign mitochondrial reads using GATK
##python MitoSort_pipeline.py -b /path/to/possorted_bam.bam -f /path/to/reference.fasta --gatk_path /path/to/GenomeAnalysisTK_3.5-0.jar -o /path/to/output_dir
fasta_path="/md01/nieyg/ref/10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
gatk_path="/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/MitoSort/GenomeAnalysisTK_3.5-0.jar"
output_path="/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/"
python /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/MitoSort_pipeline_2.py mt-realign \
        -b ${bam_path} \
        -f ${fasta_path} \
        --gatk_path ${gatk_path} \
        -o ${output_path}


##Generate SNP matrices
##python MitoSort_pipeline.py generate-snp-matrix -b /path/to/possorted_chrM_realign.bam -f /path/to/reference.fasta -c /path/to/singlecell.csv -m /path/to/MitoSort/data/hg38_chrM.bed --varscan_path /path/to/VarScan.v2.3.7.jar -o /path/to/output_dir



realign_bam_file="/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/MitoSort/BAM/possorted_chrM_realign.bam" #file from the first step output
barcode_path="/md01/nieyg/project/AmbientRNA_Benchmarking/01_data/PBMC/joint_FA/joint_FA/outs/per_barcode_metrics.csv"
chrM_bed="/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/MitoSort/data/hg38_chrM.bed" #containing chrM region
varscan_path="/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/MitoSort/VarScan.v2.3.7.jar"
output_path="/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/"

python /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/MitoSort/MitoSort_pipeline.py generate-snp-matrix \
        -b ${realign_bam_file} \
        -f ${fasta_path} \
        -c ${barcode_path} \
        -m ${chrM_bed} \
        --cell_tag CB \
        --varscan_path ${varscan_path} \
        -o ${output_path}

# --cell_tag CB if your bamdile barcode tag is CB, when is OR change it to CR

##Doublet identification and sample demultiplexing
##python MitoSort_pipeline.py demultiplex -o /path/to/output_dir -k number_of_pooled_individuals
output_path="/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/"
python /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/MitoSort/MitoSort_pipeline.py demultiplex \
        -o ${output_path} \
        -k 3 \
        --p1_cutoff 0.8 \
        --p2_cutoff 0.2
# -k 'the sample number pooled'
# --p1_cutoff 'maximum cutoff of p1 for doublet identification.'
# --p2_cutoff 'minimun cutoff of p2 for doublet identification.'







############################################ mitosort ############################################

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



specific_germline <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/MitoSort/Demultiplex_output_4/specific_germline.txt')
specific_germline <- as.data.frame(specific_germline)
colnames(specific_germline)[1] <- 'Specific_germline'

sample_result <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/MitoSort/Demultiplex_output_4/result_pvalue.txt')
sample_result <- as.data.frame(sample_result)

SNP_matrix <- read.csv('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/MitoSort/SNP_matrix/alt.csv',row.names=1)
colnames(SNP_matrix) <- gsub('\\.','-',colnames(SNP_matrix))

filtered_cell <- sample_result$Barcode[!sample_result$Demultiplex %in% c('Doublet','Unassign')]

ref_matrix <- read.csv('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/Mitosort_output_new/MitoSort/SNP_matrix/ref.csv',row.names=1)
colnames(ref_matrix) <- gsub('\\.','-',colnames(ref_matrix))

all_matrix <- SNP_matrix + ref_matrix

specific_SNP_matrix <- SNP_matrix[specific_germline$Specific_germline,filtered_cell]
specific_all_matrix <- all_matrix[specific_germline$Specific_germline,filtered_cell]

table(rownames(specific_SNP_matrix) == rownames(specific_SNP_matrix))
table(colnames(specific_all_matrix) == colnames(specific_all_matrix))

data <- specific_SNP_matrix/specific_all_matrix
data[is.na(data)] <- 0
count <- data


type <- factor(sample_result[match(colnames(count),sample_result$Barcode),'Demultiplex'])
annotation_col = data.frame(Sample_type = type)
rownames(annotation_col) = factor(colnames(count))

type <- factor(specific_germline[match(rownames(count),specific_germline$Specific_germline),'Sample'])
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

snv_order <- c(specific_germline$Specific_germline[specific_germline$Sample=='Sample0'],
           specific_germline$Specific_germline[specific_germline$Sample=='Sample1'],
           specific_germline$Specific_germline[specific_germline$Sample=='Sample2'],
           specific_germline$Specific_germline[specific_germline$Sample=='Sample3'])

count <- count[snv_order,sample_order]   
bk <- c(seq(0,1,by=0.01))

pdf(paste0('./',"mitosort.pdf"))
p3<-pheatmap(count,cluster_cols = F,cluster_rows = F,
        color = c(colorRampPalette(colors = c("#FDF6F3","#5C191A"))(length(bk))),
        legend_breaks=seq(0,1,0.25),
        breaks=bk,
        annotation_col = annotation_col, 
        #annotation_row = annotation_row,
        annotation_colors = ann_colors,
        show_rownames=F,show_colnames=F)
print(p3)
dev.off()





