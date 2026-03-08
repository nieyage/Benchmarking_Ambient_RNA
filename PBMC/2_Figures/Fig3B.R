cd /data/R03/chenxy957/project//Miss.nie/contamination/04_human_mix_joint_FA/04_4_bam
conda activate gatk

# ps aux | awk '{print $1,$4}' | awk '{arr[$1]+=$2} END {for (i in arr) {print i,arr[i]"%"}}'
# specify the number of cores to use
cores=8

bam_1='/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/04_4_bam/rawdata/Sample0'
bam_1='/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/04_4_bam/rawdata/Sample1'
bam_1='/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/04_4_bam/rawdata/Sample2'
bam_1='/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/04_4_bam/rawdata/Sample3'

samplename=`basename ${bam_1}`
echo "Sample name is $samplename"      

# directory with the genome and transcriptome index files + name of the gene annotation file
genome=/data/R03/chenxy957/reference/STAR-index/hg38/hg38
gtf=/data/R03/chenxy957/reference/STAR-index/hg38/gencode.v48.annotation.gtf
ref_genome=/data/R03/chenxy957/reference/STAR-index/hg38/GRCh38.p14.genome.fa
dbSNP_vcf=/data/R03/chenxy957/reference/snp/hg38/00-All.vcf
new_ref_genome=/data/R03/chenxy957/project/Miss.nie/contamination/09_gatk_learning/new_fa/genome.nochr.fa
hg38db=/data/R03/chenxy957/project/Miss.nie/contamination/09_gatk_learning/9_Variants_Annotation/humandb

# make all of the output directories
output_dir=/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/04_4_bam/
mkdir -p ${output_dir}3_picard
mkdir -p ${output_dir}4_SplitNCigarReads
mkdir -p ${output_dir}5_BaseRecalibrator
mkdir -p ${output_dir}6_ApplyBQSR
mkdir -p ${output_dir}7_HaplotypeCaller
mkdir -p ${output_dir}8_gVCF
mkdir -p ${output_dir}9_Variants_Annotation

# set up output directories
picard_out=${output_dir}3_picard/
SplitNCigarReads_out=${output_dir}4_SplitNCigarReads/
BaseRecalibrator_out=${output_dir}5_BaseRecalibrator/
ApplyBQSR_out=${output_dir}6_ApplyBQSR/
HaplotypeCaller_out=${output_dir}7_HaplotypeCaller/
gVCF_out=${output_dir}8_gVCF/
Variants_Annotation_out=${output_dir}9_Variants_Annotation/

## No.3 data cleanup.
echo "Starting picard data cleanup for $samplename"
picard AddOrReplaceReadGroups I=${bam_1}_filtered.bam O=${picard_out}${samplename}_rg.tmp RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample
picard MarkDuplicates I=${picard_out}${samplename}_rg.tmp O=${picard_out}${samplename}_markdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${picard_out}${samplename}.metrics


# No.4 cigar n reads 

# cd /data/R03/chenxy957/project/Miss.nie/contamination/09_gatk_learning/5_BaseRecalibrator
# wget https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz
# bwa index -a bwtsw Homo_sapiens_assembly38.fasta
# samtools faidx Homo_sapiens_assembly38.fasta
# picard CreateSequenceDictionary \
#   R=Homo_sapiens_assembly38.fasta \
#   O=Homo_sapiens_assembly38.dict
echo "Starting SplitNCigarReads for $samplename"
gatk SplitNCigarReads \
  -R ${ref_genome} \
  -I ${picard_out}${samplename}_markdup.bam \
  -O ${SplitNCigarReads_out}${samplename}_split.bam


# No.5 Base Quality Recalibration
echo "Starting Base Quality Recalibration for $samplename"
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz.md5

# gatk IndexFeatureFile -I ${dbSNP_vcf} 
# gatk IndexFeatureFile -I $INDEL_vcf




# No.6 Base Quality Recalibration
samtools view -h ${SplitNCigarReads_out}${samplename}_split.bam \
  | awk '{if($0 ~ /^@/){if($0 ~ /^@SQ/ && $2 !~ /^SN:chr/) next; else print} else {if($3 ~ /^chr/) print}}' \
  | samtools view -b -o ${SplitNCigarReads_out}${samplename}_chr_split.bam

samtools view -H ${SplitNCigarReads_out}${samplename}_chr_split.bam > ${SplitNCigarReads_out}${samplename}_chr_split_header.sam

sed '/^@SQ/{s/SN:chr\([0-9XY]\)/SN:\1/g;s/SN:chrM/SN:MT/g}' ${SplitNCigarReads_out}${samplename}_chr_split_header.sam  > ${SplitNCigarReads_out}${samplename}_chr_split_new.header.sam

samtools reheader ${SplitNCigarReads_out}${samplename}_chr_split_new.header.sam ${SplitNCigarReads_out}${samplename}_chr_split.bam > ${SplitNCigarReads_out}${samplename}_chr_split_output.bam

samtools index ${SplitNCigarReads_out}${samplename}_chr_split_output.bam

samtools idxstats ${SplitNCigarReads_out}${samplename}_chr_split_output.bam | cut -f1 


cd /data/R03/chenxy957/project/Miss.nie/contamination/09_gatk_learning/new_fa

sed -e 's/^>chr/>/' -e 's/^>M/>MT/' ${ref_genome} > genome.nochr.fa

samtools faidx genome.nochr.fa

gatk CreateSequenceDictionary \
   -R /data/R03/chenxy957/project/Miss.nie/contamination/09_gatk_learning/new_fa/genome.nochr.fa \
   -O /data/R03/chenxy957/project/Miss.nie/contamination/09_gatk_learning/new_fa/genome.nochr.dict


gatk BaseRecalibrator \
  -R ${new_ref_genome} \
  -I ${SplitNCigarReads_out}${samplename}_chr_split_output.bam \
  -O ${BaseRecalibrator_out}${samplename}_recal_table \
  -known-sites ${dbSNP_vcf} 

gatk ApplyBQSR \
  -R ${new_ref_genome} \
  -I ${SplitNCigarReads_out}${samplename}_chr_split_output.bam \
  --bqsr-recal-file ${BaseRecalibrator_out}${samplename}_recal_table \
  -O ${ApplyBQSR_out}${samplename}_recal.bam


# No.7 Variant Calling
gatk HaplotypeCaller \
  -R ${new_ref_genome} \
  -I ${ApplyBQSR_out}${samplename}_recal.bam \
  -O ${HaplotypeCaller_out}${samplename}_variants.vcf \
  --dont-use-soft-clipped-bases \
  -stand-call-conf 20.0 \
  --dbsnp ${dbSNP_vcf}

# No.8 Variant filter
gatk VariantFiltration \
  -R ${new_ref_genome} \
  -V ${HaplotypeCaller_out}${samplename}_variants.vcf \
  -O ${HaplotypeCaller_out}${samplename}_filtered_variants.vcf \
  --window 35 \
  --cluster 3 \
  --filter-name "FS" \
  --filter "FS>30.0" \
  --filter-name "QD" \
  --filter "QD<2.0"


# No.9 Variants_Annotation
gatk SelectVariants \
   --select-type-to-include SNP \
   -R ${new_ref_genome} \
   -V ${HaplotypeCaller_out}${samplename}_filtered_variants.vcf \
   -O ${Variants_Annotation_out}${samplename}_fil_snp.vcf

/data/R03/chenxy957/project/Miss.nie/contamination/09_gatk_learning/9_Variants_Annotation/annovar/convert2annovar.pl  \
   --format vcf4 ${Variants_Annotation_out}${samplename}_fil_snp.vcf > ${Variants_Annotation_out}${samplename}_fil_snp.avinput

/data/R03/chenxy957/project/Miss.nie/contamination/09_gatk_learning/9_Variants_Annotation/annovar/table_annovar.pl \
  ${Variants_Annotation_out}${samplename}_fil_snp.avinput \
  ${hg38db} \
  -otherinfo \
  --build hg38 \
  -out ${Variants_Annotation_out}${samplename}_fil_snp.anno \
  -protocol refGene,cytoBand \
  -operation g,r \
  -remove \
  -nastring '.'


##################################################################################################################
cd /data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/09_contamination_gene_ratio
conda activate new_r4_base
R

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(cowplot)


S0anno <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/04_4_bam/9_Variants_Annotation/Sample0_fil_snp.anno.hg38_multianno.txt')
S0anno <- as.data.frame(S0anno)
S0anno$snv_ID <- paste0('chr',S0anno$'Chr',':',S0anno$Start,":",S0anno$Ref,'-',S0anno$Alt)
S0anno <- S0anno[,c('snv_ID','Gene.refGene')]

S1anno <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/04_4_bam/9_Variants_Annotation/Sample1_fil_snp.anno.hg38_multianno.txt')
S1anno <- as.data.frame(S1anno)
S1anno$snv_ID <- paste0('chr',S1anno$'Chr',':',S1anno$Start,":",S1anno$Ref,'-',S1anno$Alt)
S1anno <- S1anno[,c('snv_ID','Gene.refGene')]

S2anno <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/04_4_bam/9_Variants_Annotation/Sample2_fil_snp.anno.hg38_multianno.txt')
S2anno <- as.data.frame(S2anno)
S2anno$snv_ID <- paste0('chr',S2anno$'Chr',':',S2anno$Start,":",S2anno$Ref,'-',S2anno$Alt)
S2anno <- S2anno[,c('snv_ID','Gene.refGene')]

S3anno <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/04_4_bam/9_Variants_Annotation/Sample3_fil_snp.anno.hg38_multianno.txt')
S3anno <- as.data.frame(S3anno)
S3anno$snv_ID <- paste0('chr',S3anno$'Chr',':',S3anno$Start,":",S3anno$Ref,'-',S3anno$Alt)
S3anno <- S3anno[,c('snv_ID','Gene.refGene')]


all_anno <- rbind(S0anno, S1anno, S2anno, S3anno)
all_anno <- unique(all_anno)
all_anno_expanded <- all_anno %>%
  separate_rows(Gene.refGene, sep = ",")
all_anno_expanded <- as.data.frame(all_anno_expanded)
all_anno_expanded <- unique(all_anno_expanded)

##################################################################################################################

S0_vcf <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/04_4_bam/7_HaplotypeCaller/Sample0_filtered_variants.vcf')
S1_vcf <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/04_4_bam/7_HaplotypeCaller/Sample1_filtered_variants.vcf')
S2_vcf <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/04_4_bam/7_HaplotypeCaller/Sample2_filtered_variants.vcf')
S3_vcf <- data.table::fread('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/04_4_bam/7_HaplotypeCaller/Sample3_filtered_variants.vcf')

S0_vcf <- as.data.frame(S0_vcf)
S1_vcf <- as.data.frame(S1_vcf)
S2_vcf <- as.data.frame(S2_vcf)
S3_vcf <- as.data.frame(S3_vcf)

S0_vcf$snv_ID <- paste0('chr',S0_vcf$'#CHROM',':',S0_vcf$POS,":",S0_vcf$REF,'-',S0_vcf$ALT)
S1_vcf$snv_ID <- paste0('chr',S1_vcf$'#CHROM',':',S1_vcf$POS,":",S1_vcf$REF,'-',S1_vcf$ALT)
S2_vcf$snv_ID <- paste0('chr',S2_vcf$'#CHROM',':',S2_vcf$POS,":",S2_vcf$REF,'-',S2_vcf$ALT)
S3_vcf$snv_ID <- paste0('chr',S3_vcf$'#CHROM',':',S3_vcf$POS,":",S3_vcf$REF,'-',S3_vcf$ALT)

# mut read
S0_vcf$mut_reads_SO <- sapply(strsplit(as.character(S0_vcf$sample), ":"), function(x) {
  AD <- x[2]
  as.numeric(strsplit(AD, ",")[[1]][2])
})


S1_vcf$mut_reads_S1 <- sapply(strsplit(as.character(S1_vcf$sample), ":"), function(x) {
  AD <- x[2]
  as.numeric(strsplit(AD, ",")[[1]][2])
})


S2_vcf$mut_reads_S2 <- sapply(strsplit(as.character(S2_vcf$sample), ":"), function(x) {
  AD <- x[2]
  as.numeric(strsplit(AD, ",")[[1]][2])
})


S3_vcf$mut_reads_S3 <- sapply(strsplit(as.character(S3_vcf$sample), ":"), function(x) {
  AD <- x[2]
  as.numeric(strsplit(AD, ",")[[1]][2])
})


S0_vcf$samplt_ID <- 'S0'
S1_vcf$samplt_ID <- 'S1'
S2_vcf$samplt_ID <- 'S2'
S3_vcf$samplt_ID <- 'S3'

S0_vcf_brief <- S0_vcf[,c("snv_ID",'mut_reads_SO')]
S1_vcf_brief <- S1_vcf[,c("snv_ID",'mut_reads_S1')]
S2_vcf_brief <- S2_vcf[,c("snv_ID",'mut_reads_S2')]
S3_vcf_brief <- S3_vcf[,c("snv_ID",'mut_reads_S3')]


# combine
merged_df <- merge(S0_vcf_brief, S1_vcf_brief, by = "snv_ID", all = TRUE)
merged_df <- merge(merged_df, S2_vcf_brief, by = "snv_ID", all = TRUE)
merged_df <- merge(merged_df, S3_vcf_brief, by = "snv_ID", all = TRUE)
head(merged_df)
merged_df[is.na(merged_df)] <- 0
alt_df <- merged_df
rownames(alt_df) <- alt_df$snv_ID
alt_df$snv_ID <- NULL
head(alt_df)

mut_cols <- c("mut_reads_SO", "mut_reads_S1", "mut_reads_S2", "mut_reads_S3")

# Calculate the proportion by row
prop_matrix <- t(apply(merged_df[, mut_cols], 1, function(x) {
  total <- sum(x)
  if(total == 0) {
    rep(0, length(x))
  } else {
    x / total
  }
}))

# Convert to data frame
prop_df <- as.data.frame(prop_matrix)
colnames(prop_df) <- mut_cols
row.names(prop_df)<-merged_df$snv_ID

mut_cols <- c("mut_reads_SO", "mut_reads_S1", "mut_reads_S2", "mut_reads_S3")

# Add a new column named "max_reads" to store the maximum value for each row.
prop_df$max_proportion <- apply(prop_df[, mut_cols], 1, max)

selected_snv <- rownames(prop_df)[prop_df$max_proportion > 0.85 & prop_df$max_proportion < 1] # 4193

alt_df_selected <- alt_df[selected_snv,]
colnames(alt_df_selected) <- paste0('S',0:3)
alt_df_selected$max_col <- apply(
  alt_df_selected[, c("S0","S1","S2","S3")],
  1,
  function(x) names(x)[which.max(x)]
)
conta_alt_info <- alt_df_selected


############################################################################################################
S0_vcf <- S0_vcf %>%
  mutate(dp_S0 = as.numeric(sapply(strsplit(sample, ":"), function(x) x[3])))
S1_vcf <- S1_vcf %>%
  mutate(dp_S1 = as.numeric(sapply(strsplit(sample, ":"), function(x) x[3])))
S2_vcf <- S2_vcf %>%
  mutate(dp_S2 = as.numeric(sapply(strsplit(sample, ":"), function(x) x[3])))
S3_vcf <- S3_vcf %>%
  mutate(dp_S3 = as.numeric(sapply(strsplit(sample, ":"), function(x) x[3])))


S0_dp_brief <- S0_vcf[,c("snv_ID",'dp_S0')]
S1_dp_brief <- S1_vcf[,c("snv_ID",'dp_S1')]
S2_dp_brief <- S2_vcf[,c("snv_ID",'dp_S2')]
S3_dp_brief <- S3_vcf[,c("snv_ID",'dp_S3')]

# combine
merged_df <- merge(S0_dp_brief, S1_dp_brief, by = "snv_ID", all = TRUE)
merged_df <- merge(merged_df, S2_dp_brief, by = "snv_ID", all = TRUE)
merged_df <- merge(merged_df, S3_dp_brief, by = "snv_ID", all = TRUE)
head(merged_df)
merged_df[is.na(merged_df)] <- 0

dp_df <- merged_df
rownames(dp_df) <- dp_df$snv_ID
dp_df$snv_ID <- NULL
head(dp_df)



dp_df_selected <- dp_df[selected_snv,] 
alt_df_selected <- alt_df[selected_snv,]
colnames(alt_df_selected) <- paste0('S',0:3)
colnames(dp_df_selected) <- paste0('S',0:3)

conta_snv_ratio <- c()
conta_snv_alt <- c()
conta_snv_dp <- c()


for (i in 1:length(selected_snv)) {
    print(selected_snv[i])
    current_snv <- selected_snv[i]
    current_sample <- conta_alt_info[current_snv,]$max_col
    print(current_sample)
    
    needed_sample <- setdiff(paste0('S',0:3),current_sample)
    all_alt_contamination_reads <- sum(alt_df_selected[current_snv,needed_sample])
    all_dp_contamination_reads <- sum(dp_df_selected[current_snv,needed_sample])
    current_ratio <- all_alt_contamination_reads/all_dp_contamination_reads
    conta_snv_ratio[i] <- current_ratio
    conta_snv_alt[i] <- all_alt_contamination_reads
    conta_snv_dp[i] <- all_dp_contamination_reads
}

conta_snv_ratio_df <- data.frame(snv_ID = selected_snv,
                                 conta_alt_reads = conta_snv_alt,
                                 conta_dp_reads = conta_snv_dp,
                                 conta_ratio = conta_snv_ratio)
            
conta_snv_ratio_df_filter <- conta_snv_ratio_df[conta_snv_ratio_df$conta_dp_reads > 10,]
dim(conta_snv_ratio_df_filter) #513

conta_snv_ratio_df_filter_anno <- merge(conta_snv_ratio_df_filter,all_anno_expanded,by='snv_ID')

conta_snv_ratio_df_filter_anno_sorted <- conta_snv_ratio_df_filter_anno %>%
  arrange(desc(conta_ratio)) %>%
  distinct(Gene.refGene, .keep_all = TRUE)

head(conta_snv_ratio_df_filter_anno_sorted)


prop_df$snv_ID <- rownames(prop_df)

conta_snv_ratio_df_filter_anno_sorted <- merge(conta_snv_ratio_df_filter_anno_sorted,prop_df[selected_snv,c('snv_ID','max_proportion')],by='snv_ID')
conta_snv_ratio_df_filter_anno_sorted <- conta_snv_ratio_df_filter_anno_sorted[order(conta_snv_ratio_df_filter_anno_sorted$conta_ratio,decreasing=T),]
dim(conta_snv_ratio_df_filter_anno_sorted) #441
write.csv(conta_snv_ratio_df_filter_anno_sorted,'conta_snv_ratio_df_filter_anno_sorted.csv')



conta_alt_info$snv_ID <- rownames(conta_alt_info)
conta_alt_info_anno <- merge(conta_alt_info,all_anno_expanded,by='snv_ID')
write.csv(conta_alt_info_anno,'conta_alt_info_anno.csv')


conta_gene_ratio <- read.csv('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/09_contamination_gene_ratio/conta_snv_ratio_df_filter_anno_sorted.csv',row.names=1)
top_gene <- c(conta_gene_ratio$Gene.refGene[1:3],top_marker_gene)
conta_gene_ratio$rank <- 1:nrow(conta_gene_ratio)
conta_gene_ratio$gene <- conta_gene_ratio$Gene.refGene
conta_gene_ratio$contamination_ratio <- conta_gene_ratio$conta_ratio


# Load the necessary packages
library(ggplot2)
library(ggrepel)

pdf('rank_dotplot.pdf', width = 5, height = 5)

p <- ggplot(conta_gene_ratio, aes(x = rank, y = contamination_ratio)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_point(
    data = subset(conta_gene_ratio, gene %in% top_marker_gene),
    color = "red", alpha = 0.8
  ) +
  geom_text_repel(
    data = subset(conta_gene_ratio, gene %in% top_gene),
    aes(label = gene),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.5,
    max.overlaps = Inf
  ) +
  labs(
    title = "Contamination Gene Ratio by Rank",
    x = "Rank",
    y = "Ratio",
    size = "Ratio"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

print(p)

dev.off()