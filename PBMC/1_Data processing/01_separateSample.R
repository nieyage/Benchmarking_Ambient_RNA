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

# Changes made to MitoSort_pipeline_2.py, otherwise possorted_chrM_realign.bam cannot be generated
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


#lib5
#step1
python ./MitoSort_pipeline.py generate-snp-matrix -b /md01/Chenzh275/Project/Zhangwx/lib1/MitoSort/BAM/possorted_chrM_realign.bam -f /public/home/chenbzh5/DB/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa -c /data/R03/zhangwx/project/human_PBMC/lib1/outs/per_barcode_metrics.csv -m /md01/Chenzh275/Project/Zhangwx/lib1/MitoSort/data/hg38_chrM.bed --varscan_path /md01/Chenzh275/Project/Zhangwx/lib1/MitoSort/VarScan.v2.3.7.jar --cell_tag CB -o /md01/Chenzh275/Project/Zhangwx/lib1/
#step2
python ./MitoSort_pipeline.py generate-snp-matrix -b /md01/Chenzh275/Project/Zhangwx/lib5/MitoSort/BAM/possorted_chrM_realign.bam -f /public/home/chenbzh5/DB/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa -c /data/R03/zhangwx/project/human_PBMC/NhPBMC/NhPBMC_joint/outs/per_barcode_metrics.csv -m /md01/Chenzh275/Project/Zhangwx/lib5/MitoSort/data/hg38_chrM.bed --varscan_path /md01/Chenzh275/Project/Zhangwx/lib5/MitoSort/VarScan.v2.3.7.jar --cell_tag CB -o /md01/Chenzh275/Project/Zhangwx/lib5/
#step3
python ./MitoSort_pipeline.py demultiplex -o /md01/Chenzh275/Project/Zhangwx/lib5 -k 4 --p1_cutoff 0.8 --p2_cutoff 0.2

#lib1
#step1
python ./MitoSort_pipeline.py mt-realign -b /data/R03/zhangwx/project/human_PBMC/lib1/outs/atac_possorted_bam.bam -f /public/home/chenbzh5/DB/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa --gatk_path /md01/zhangwx/project/separeteLib1/MitoSort/GenomeAnalysisTK_3.5-0.jar -o /md01/zhangwx/project/separeteLib1
#step2
python ./MitoSort_pipeline.py generate-snp-matrix -b /md01/zhangwx/project/separeteLib1/MitoSort/BAM/possorted_chrM_realign.bam -f /public/home/chenbzh5/DB/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa -c /data/R03/zhangwx/project/human_PBMC/lib1/outs/per_barcode_metrics.csv -m /md01/Chenzh275/Project/Zhangwx/lib1/MitoSort/data/hg38_chrM.bed --varscan_path /md01/zhangwx/project/separeteLib1/MitoSort/VarScan.v2.3.7.jar --cell_tag CB -o /md01/zhangwx/project/separeteLib1/
#step3
#routine pipeline for demultiplex sample
#python ./MitoSort_pipeline.py demultiplex -o /md01/zhangwx/project/separeteLib1/ -k 2 --p1_cutoff 0.8 --p2_cutoff 0.2 

#because the minimal cell number of lib1 we should add the parmeter "--method 'direct'"
python ./MitoSort_pipeline.py demultiplex -o /md01/zhangwx/project/separeteLib1/ -k 2 --p1_cutoff 0.7 --p2_cutoff 0.3 --method 'direct'