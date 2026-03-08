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
cell_expression_ratio <- cell_expression_ratio[cell_expression_ratio$gene == 'VCAN',]


cell_expression_ratio$cell<-factor(cell_expression_ratio$cell,levels=c('Mono','NK',"DC","CD4_T","Naive_CD8_T","Memory_CD8_T","B_cell","plasma"))
cell_expression_ratio$tool <- factor(cell_expression_ratio$tool,levels=c('Rawdata','scAR'))

# 你的数据框 cell_expression_ratio 已经有两列: tool 和 expression_ratio

pdf('cell_expression_ratio.pdf', width = 5, height = 4)

p <- ggplot(cell_expression_ratio, aes(x = cell, y = expression_ratio, fill = tool)) +
  geom_col(position = "dodge", width = 0.7) +
  labs(title = "Cell Expression Ratio by Tool",
       x = "Tools",
       y = "Percentage",
       fill = "Cell Type") +
  scale_fill_manual(values = c("Rawdata" = "#BFBFBF", "scAR" = "#0052A1")) +
  scale_y_continuous(labels = percent_format(accuracy = 0.01)) +   # 转百分比
  theme_classic() +                                             # 去掉背景框和网格
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


print(p)
dev.off()




#################################################################################################################################################
cell_umi_ratio <-  read.csv('/data/R03/chenxy957/project/Miss.nie/contamination/04_human_mix_joint_FA/13_scAR/03_anno/all_umi_ratio.csv',header=T, row.names=1)
cell_umi_ratio <- cell_umi_ratio[cell_umi_ratio$gene == 'VCAN',]

cell_umi_ratio$cell<-factor(cell_umi_ratio$cell,levels=c('Mono','NK',"DC","CD4_T","Naive_CD8_T","Memory_CD8_T","B_cell","plasma"))
cell_umi_ratio$tool <- factor(cell_umi_ratio$tool,levels=c('Rawdata','scAR'))

# 你的数据框 cell_umi_ratio 已经有两列: tool 和 umi_ratio

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

print(p)

dev.off()


