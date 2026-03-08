library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(readxl)

integrated_data <- readRDS("/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/integrated_8samples_filtered.rds")

# Load OR gene information
mm10_genes <- read.table("/data/R02/tanj93/reference/mouse_ORs/gene_symbols.txt", header = FALSE)
OR_genes <- read_excel("/data/R02/tanj93/reference/mouse_ORs/PMC7055050_mouse_OR_genes.xlsx") %>% 
  select("Gene symbol", "Ensembl gene ID") %>% 
  distinct()

merged_OR_genes <- merge(mm10_genes, OR_genes, by.x = "V2", by.y = "Ensembl gene ID")
or_gene_list <- merged_OR_genes$V1

get_single_OR <- function(seurat_obj, or_genes, threshold = 4) {
  available_OR <- intersect(or_genes, rownames(seurat_obj))
  OR_expr <- LayerData(seurat_obj, assay = "RNA", layer = "counts")[available_OR, ]
  
  apply(OR_expr, 2, function(x) {
    expressed <- which(x >= threshold)
    if(length(expressed) == 1) {
      names(expressed)[1]
    } else {
      NA
    }
  })
}

# Extract mature OSNs
Idents(integrated_data) <- "celltype"
mature_OSNs <- subset(integrated_data, idents = "Mature OSNs")

# Set default assay and process data
DefaultAssay(mature_OSNs) <- "RNA"

# Identify single OR expressing cells
mature_OSNs$top_OR <- get_single_OR(mature_OSNs, or_gene_list, threshold = 4)
mature_OSNs_filtered <- subset(mature_OSNs, !is.na(top_OR))

# Normalization and scaling
mature_OSNs_filtered <- NormalizeData(mature_OSNs_filtered)
mature_OSNs_filtered <- FindVariableFeatures(mature_OSNs_filtered, nfeatures = 2000)
mature_OSNs_filtered <- ScaleData(mature_OSNs_filtered)
mature_OSNs_filtered <- RunPCA(mature_OSNs_filtered, npcs = 50)

focus_ORs <- c("Olfr728", "Olfr1226", "Olfr1440")
focus_colors <- c("#FF7AA3", "#9D1A42", "#F46C44")

mature_OSNs_filtered <- RunUMAP(mature_OSNs_filtered,
                               dims = 1:40,
                               n.neighbors = 50,
                               min.dist = 0.3)

plot_data <- as.data.frame(mature_OSNs_filtered[["umap"]]@cell.embeddings)
colnames(plot_data) <- c("UMAP_1", "UMAP_2")
plot_data$OR <- mature_OSNs_filtered$top_OR
plot_data$group <- ifelse(plot_data$OR %in% focus_ORs, plot_data$OR, "Other")

label_positions <- plot_data %>%
  filter(group != "Other") %>%
  group_by(group) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  ) %>%
  mutate(
    UMAP_1 = UMAP_1 + c(0, 0.5, -0.5)[1:n()],
    UMAP_2 = UMAP_2 + c(0.5, 0, -0.5)[1:n()]
  ) %>%
  mutate(
    UMAP_2 = case_when(
      group == "Olfr728" ~ UMAP_2 - 0.4,
      group == "Olfr1226" ~ UMAP_2 + 0.4,
      TRUE ~ UMAP_2
    )
  )

# Create UMAP plot
umap_plot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  
  geom_point(data = subset(plot_data, group == "Other"),
            color = "gray90", size = 1, alpha = 0.4) +
  
  geom_point(
    data = subset(plot_data, group != "Other"),
    aes(fill = group),
    size = 1.8,
    alpha = 0.8,
    stroke = 0.5,
    shape = 21,
    color = "gray40"
  ) +
  
  scale_fill_manual(values = setNames(focus_colors, focus_ORs),
                   name = "OR Genes",
                   breaks = focus_ORs) +
  
  geom_text(
    data = label_positions,
    aes(label = group),
    color = "black",
    size = 4,
    fontface = "plain"
  ) +
  
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "right",
    aspect.ratio = 1
  )

# Save UMAP plot
pdf(paste0(output_dir, "OR_UMAP.pdf"), width = 7, height = 6, useDingbats = FALSE)
print(umap_plot)
dev.off()

png(paste0(output_dir, "OR_UMAP.png"), width = 2100, height = 1800, res = 300)
print(umap_plot)
dev.off()

target_genes <- c("Olfr728", "Olfr1226", "Olfr1440")

output_dir_feature <- "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/OR/"
if(!dir.exists(output_dir_feature)) dir.create(output_dir_feature, recursive = TRUE)

# Feature plots
featureplots <- function(seurat_obj, gene) {
  FeaturePlot(seurat_obj,
             features = gene,
             reduction = "umap",
             pt.size = 0.8,
             order = TRUE,
             cols = c("lightgrey", "red"),
             raster = FALSE) +
    ggtitle(gene) +
    xlab("umap_1") +
    ylab("umap_2") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      legend.position = "right",
      legend.title = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )
}

for(gene in target_genes) {
  # Generate plot
  p <- featureplots(mature_OSNs_filtered, gene)
  
  # Save as PDF
  pdf(paste0(output_dir_feature, gene, "_featureplot.pdf"), width = 7, height = 6, useDingbats = FALSE)
  print(p)
  dev.off()
  
  # Save as PNG
  png(paste0(output_dir_feature, gene, "_featureplot.png"), width = 2100, height = 1800, res = 300)
  print(p)
  dev.off()
}
