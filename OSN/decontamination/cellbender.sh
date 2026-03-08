cellbender remove-background \
  --cuda \
  --input /data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_1_merged/GSM5608794/outs/raw_feature_bc_matrix.h5 \
  --output /data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/NC1/cellbender_output_file.h5

cellbender remove-background \
  --cuda \
  --input /data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_2_merged/GSM5608795/outs/raw_feature_bc_matrix.h5 \
  --output /data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/NC2/cellbender_output_file.h5

cellbender remove-background \
  --cuda \
  --input /data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_3_merged/GSM5608796/outs/raw_feature_bc_matrix.h5 \
  --output /data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/NC3/cellbender_output_file.h5

cellbender remove-background \
  --cuda \
  --input /data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/NC_2022/Project_OSN_scRNAseq/Library_4_merged/GSM5608797/outs/raw_feature_bc_matrix.h5 \
  --output /data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/NC4/cellbender_output_file.h5

cellbender remove-background \
  --cuda \
  --input /data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752985/cellranger/GSM4752985/outs/raw_feature_bc_matrix.h5 \
  --output /data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/SA3/cellbender_output_file.h5

cellbender remove-background \
  --cuda \
  --input /data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752986/cellranger/GSM4752986/outs/raw_feature_bc_matrix.h5 \
  --output /data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/SA4/cellbender_output_file.h5

cellbender remove-background \
  --cuda \
  --input /data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752989/cellranger/GSM4752989/outs/raw_feature_bc_matrix.h5 \
  --output /data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/SA5/cellbender_output_file.h5

cellbender remove-background \
  --cuda \
  --input /data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/01_rawdata/GSE157100/GSM4752990/cellranger/GSM4752990/outs/raw_feature_bc_matrix.h5 \
  --output /data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/SA6/cellbender_output_file.h5  这可以优化吗