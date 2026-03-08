### Fig5J
library(Seurat)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(reshape2)
library(proxy)
library(ggpubr)
library(dplyr)
library(writexl)
library(lsa)
library(readxl)

# Set color
light_blue <- "#A6CEE3"
light_red <- "#FB9A99"

main_output_dir <- "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/"
dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)

mm10_genes <- read.table("/data/R02/tanj93/reference/mouse_ORs/gene_symbols.txt", header = FALSE)
OR_genes <- read_xlsx("/data/R02/tanj93/reference/mouse_ORs/PMC7055050_mouse_OR_genes.xlsx") %>%
  select("Gene symbol", "Ensembl gene ID") %>%
  distinct()

merged_OR_genes <- merge(mm10_genes, OR_genes, 
                         by.x = "V2", by.y = "Ensembl gene ID",
                         all = FALSE)

all_or_genes <- unique(merged_OR_genes$`Gene symbol`)

target_genes <- c("Olfr1440", "Olfr166", "Olfr1384")
threshold_value <- 4

get_top_OR <- function(seurat_obj, or_genes, threshold = 4) {
  available_OR <- intersect(or_genes, rownames(seurat_obj))
  if(length(available_OR) == 0) {
    cat("Warning: No OR genes found in the dataset\n")
    return(NULL)
  }
  
  OR_expr <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[available_OR, ]
  
  apply(OR_expr, 2, function(x) {
    expressed <- which(x >= threshold)
    if(length(expressed) == 1) {
      names(x)[expressed]
    } else if(length(expressed) > 1) {
      names(x)[which.max(x[expressed])]
    } else {
      NA
    }
  })
}

combine_all_layers <- function(seurat_obj, dataset_name) {
  layers_to_combine <- paste0("data.", 1:8)
  available_layers <- layers_to_combine[layers_to_combine %in% Layers(seurat_obj, assay = "RNA")]
  if (length(available_layers) == 0) {
    cat("No data layers found, using default counts\n")
    if ("counts" %in% Layers(seurat_obj, assay = "RNA")) {
      combined_counts <- LayerData(seurat_obj, assay = "RNA", layer = "counts")
    } else {
      stop("No counts layer found in the Seurat object")
    }
    return(combined_counts)
  }

  all_cells <- c()
  all_genes <- c()
  
  for (layer in available_layers) {
    layer_data <- LayerData(seurat_obj, assay = "RNA", layer = layer)
    all_cells <- c(all_cells, colnames(layer_data))
    all_genes <- c(all_genes, rownames(layer_data))
  }
  
  all_cells <- unique(all_cells)
  all_genes <- unique(all_genes)

  combined_matrix <- matrix(0, 
                           nrow = length(all_genes), 
                           ncol = length(all_cells),
                           dimnames = list(all_genes, all_cells))
  
  for (layer in available_layers) {
    cat("Processing layer", layer, "...\n")
    layer_data <- LayerData(seurat_obj, assay = "RNA", layer = layer)
    layer_cells <- colnames(layer_data)
    layer_genes <- rownames(layer_data)
    combined_matrix[layer_genes, layer_cells] <- as.matrix(layer_data)
  }
  
  combined_sparse <- as(combined_matrix, "dgCMatrix")
  
  return(combined_sparse)
}

analyze_dataset <- function(dataset_path, dataset_name, output_dir) {
  cat("Analyzing", dataset_name, "dataset\n")
  cat("Loading dataset from:", dataset_path, "\n")
  seurat_data <- readRDS(dataset_path)
  combined_counts <- combine_all_layers(seurat_data, dataset_name)
  seurat_obj <- CreateSeuratObject(
    counts = combined_counts,
    project = dataset_name,
    min.cells = 3,
    min.features = 200
  )

  cat("Assigning main OR gene to each cell (threshold =", threshold_value, ")...\n")
  main_or <- get_top_OR(seurat_obj, all_or_genes, threshold = threshold_value)
  if(is.null(main_or)) {
    cat("No OR genes available for analysis\n")
    return(NULL)
  }
  
  barcode_label <- data.frame(
    barcode = names(main_or),
    label = main_or
  )
  
  barcode_label <- na.omit(barcode_label)
  if(nrow(barcode_label) == 0) {
    cat("No cells have assigned OR genes\n")
    return(NULL)
  }
  
  barcode_label <- barcode_label[barcode_label$label %in% target_genes, ]
  cat("Cells expressing target OR genes:", nrow(barcode_label), "\n")
  
  if (nrow(barcode_label) == 0) {
    cat("No cells express the target OR genes!\n")
    return(NULL)
  }
  
  cat("Distribution of target OR genes:\n")
  or_counts_table <- table(barcode_label$label)
  print(or_counts_table)
  
  if (any(or_counts_table < 3)) {
    warning(paste("Some OR genes have very few cells in", dataset_name, ". Results may not be reliable."))
    cat("Cells per OR gene (min 3 recommended):\n")
    print(or_counts_table)
  }
  
  barcode_label <- barcode_label[order(barcode_label$label), ]
  seurat_subset <- subset(seurat_obj, cells = barcode_label$barcode)
  if (nrow(barcode_label) < 10) {
    cat("Using scaled gene expression for distance calculation...\n")
    seurat_subset <- NormalizeData(seurat_subset, verbose = FALSE)
    seurat_subset <- ScaleData(seurat_subset, verbose = FALSE)
    scaled_data <- GetAssayData(seurat_subset, assay = "RNA", slot = "scale.data")
    if (is.null(scaled_data) || nrow(scaled_data) == 0) {
      scaled_data <- GetAssayData(seurat_subset, assay = "RNA", layer = "counts")
      scaled_data <- log1p(scaled_data)
    }

    cat("Calculating cosine distances from gene expression...\n")
    trans_dist <- 1 - cosine(t(as.matrix(scaled_data)))
  } else {
    seurat_subset <- NormalizeData(seurat_subset, verbose = FALSE)
    seurat_subset <- FindVariableFeatures(seurat_subset, 
                                         selection.method = "vst", 
                                         nfeatures = min(2000, nrow(seurat_subset)),
                                         verbose = FALSE)
    seurat_subset <- ScaleData(seurat_subset, verbose = FALSE)
    seurat_subset <- RunPCA(seurat_subset, 
                           features = VariableFeatures(seurat_subset),
                           verbose = FALSE,
                           npcs = min(20, ncol(seurat_subset) - 1))
    n_pcs <- min(10, ncol(Embeddings(seurat_subset, "pca")))
    cat("Using", n_pcs, "PCs for distance calculation\n")
    
    embeddings <- Embeddings(seurat_subset, reduction = "pca")[, 1:n_pcs]
    
    if (!all(barcode_label$barcode == rownames(embeddings))) {
      cat("Reordering embeddings to match barcode label order...\n")
      embeddings <- embeddings[barcode_label$barcode, ]
    }
    
    cat("Calculating cosine distances...\n")
    trans_dist <- 1 - cosine(t(embeddings))
  }
  
  cat("Distance matrix dimensions:", dim(trans_dist), "\n")
  
  # Order by OR gene
  gene_order <- order(barcode_label$label)
  trans_dist_ordered <- trans_dist[gene_order, gene_order]
  
  barcode_label_pheatmap <- data.frame(
    OR = barcode_label$label
  )
  rownames(barcode_label_pheatmap) <- barcode_label$barcode
  
  or_colors <- c("Olfr1440" = "#377EB8", 
                 "Olfr166" = "#5E891B", 
                 "Olfr1384" = "#ED8F13")

  present_or_genes <- unique(barcode_label$label)
  ann_colors <- list(OR = or_colors[names(or_colors) %in% present_or_genes])
  
  # Calculate group boundaries for heatmap gaps
  or_counts_table <- table(barcode_label$label)
  group_boundaries <- cumsum(or_counts_table)[-length(or_counts_table)]
  
  # Generate heatmap exactly as specified
  p1 <- pheatmap(trans_dist_ordered,
                 cluster_cols = FALSE,
                 cluster_rows = FALSE,
                 color = colorRampPalette(c("#4B0C1A", "#C25B3F", "#F1ECEC", "#3082BD", "#1C214F"))(100),
                 annotation_col = barcode_label_pheatmap,
                 annotation_row = barcode_label_pheatmap,
                 annotation_colors = ann_colors,
                 annotation_legend = TRUE,
                 show_rownames = FALSE,
                 show_colnames = FALSE,
                 main = "Transcriptome distance",
                 legend = TRUE,
                 border_color = NA,
                 annotation_legend_side = "top",
                 legend_breaks = c(0, 0.5, 1.0, 1.5, 2.0),
                 legend_labels = c("0.0", "0.5", "1.0", "1.5", "2.0"),
                 gaps_row = group_boundaries,
                 gaps_col = group_boundaries,
                 fontsize = 10,
                 fontsize_row = 8,
                 fontsize_col = 8
  )
  
  # Calculate within-group and between-group distances
  within_group <- c()
  between_group <- c()
  
  for (i in 1:nrow(trans_dist)) {
    for (j in 1:ncol(trans_dist)) {
      if (i != j) {
        if (barcode_label[rownames(trans_dist)[i], 2] == barcode_label[colnames(trans_dist)[j], 2]) {
          within_group <- c(within_group, trans_dist[i, j])
        } else {
          between_group <- c(between_group, trans_dist[i, j])
        }
      }
    }
  }
  
  cat("Within-group distances:", length(within_group), "\n")
  cat("Between-group distances:", length(between_group), "\n")
  
  # Statistical test
  if (length(within_group) > 5 && length(between_group) > 5) {
    wilcox_test <- wilcox.test(within_group, between_group)
    p_value <- wilcox_test$p.value
    cat("Wilcoxon test p-value:", p_value, "\n")
  } else {
    warning("Not enough data points for statistical test")
    p_value <- NA
  }
  
  # Prepare plotting data
  plot_data <- data.frame(
    type = c(rep("within-OR", length(within_group)), 
             rep("between-OR", length(between_group))),
    distance = c(within_group, between_group)
  )
  plot_data$type <- factor(plot_data$type, levels = c("within-OR", "between-OR"))
  
  # Generate boxplot
  p3 <- ggplot(plot_data, aes(x = type, y = distance, fill = type)) +
    geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA) +
    scale_fill_manual(values = c("within-OR" = light_blue, "between-OR" = light_red)) +
    labs(x = "", y = "Transcriptome distance") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, hjust = 0.5),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)) +
    scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5))
  
  if (!is.na(p_value)) {
    p3 <- p3 + 
      labs(title = paste0("Wilcoxon test, p = ", format.pval(p_value, digits = 3))) +
      stat_compare_means(method = "wilcox.test", 
                        label = "p.signif",
                        label.x = 1.5,
                        label.y = 1.9,
                        size = 4)
  }
  
  # Combine all plots
  combined_plot <- plot_grid(
    p1$gtable, 
    p3,
    ncol = 2, 
    rel_widths = c(2, 1)
  )
  
  # Add title
  title <- ggdraw() + 
    draw_label(paste("Transcriptome Distance Analysis:", paste(target_genes, collapse = ", "), "-", dataset_name),
               fontface = "bold", size = 16, x = 0.5, hjust = 0.5) +
    theme(plot.margin = margin(0, 0, 10, 0))
  
  # Create final plot
  final_plot <- plot_grid(
    title, 
    combined_plot,
    ncol = 1, 
    rel_heights = c(0.08, 1),
    labels = NULL
  )
  
  # Save combined plot
  pdf(file.path(output_dir, paste0("OR_transcriptome_distance_analysis_", dataset_name, ".pdf")), 
      width = 14, height = 9)
  print(final_plot)
  dev.off()
  
  # Save individual plots
  pdf(file.path(output_dir, paste0("heatmap_", dataset_name, ".pdf")), width = 9, height = 8)
  print(p1)
  dev.off()
  
  pdf(file.path(output_dir, paste0("boxplot_", dataset_name, ".pdf")), width = 3, height = 4)
  print(p3)
  dev.off()
  
  # Save statistical results
  stats_result <- data.frame(
    Dataset = dataset_name,
    Genes = paste(target_genes, collapse = ", "),
    Threshold = threshold_value,
    Within_Group_Median = ifelse(length(within_group) > 0, median(within_group), NA),
    Between_Group_Median = ifelse(length(between_group) > 0, median(between_group), NA),
    Log2FC = ifelse(length(within_group) > 0 && length(between_group) > 0, 
                    log2(median(between_group)) - log2(median(within_group)), NA),
    P_value = p_value,
    Within_Group_Count = length(within_group),
    Between_Group_Count = length(between_group),
    Total_Cells = nrow(barcode_label),
    Cells_per_OR = paste(sapply(split(barcode_label$label, barcode_label$label), length), collapse = "; ")
  )
  
  # Save detailed data
  write.csv(plot_data, 
            file.path(output_dir, paste0("distance_data_", dataset_name, ".csv")), 
            row.names = FALSE)
  
  write.csv(barcode_label, 
            file.path(output_dir, paste0("cell_annotations_", dataset_name, ".csv")), 
            row.names = FALSE)
  
  return(list(
    stats = stats_result,
    plot_data = plot_data,
    barcode_label = barcode_label,
    trans_dist = trans_dist
  ))
}

dataset_paths <- list(Original = "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/integrated_8samples_filtered.rds")

all_results <- list()
all_stats <- data.frame()

for (dataset_name in names(dataset_paths)) {
  dataset_path <- dataset_paths[[dataset_name]]
  
  dataset_output_dir <- file.path(main_output_dir, paste0(dataset_name, "_analysis_Olfr1440_166_1384"))
  dir.create(dataset_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  tryCatch({
    result <- analyze_dataset(dataset_path, dataset_name, dataset_output_dir)
    if (!is.null(result)) {
      all_results[[dataset_name]] <- result
      all_stats <- rbind(all_stats, result$stats)
    }
  }, error = function(e) {
    cat("Error analyzing", dataset_name, ":", e$message, "\n")
  })
}

if (nrow(all_stats) > 0) {
  stats_output_dir <- file.path(main_output_dir, "Original_results_Olfr1440_166_1384")
  dir.create(stats_output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(all_stats, file.path(stats_output_dir, "Original_statistics.csv"), row.names = FALSE)
} 

### Fig5I
library(ggplot2)
library(dplyr)

fc_data <- data.frame(
  Method = c("rawdata", "CellClear", "SoupX", "CellBender", 
             "FastCAR", "scCDC", "scAR", "DecontX"),
  FoldChange = c(1.27164190275821, 1.4627110482949, 1.33598671013914, 
                 1.21448624121961, 0.97475254066223, 1.03462253485715, 
                 0.93746083422599, 1.11351710033694)
)

rawdata_value <- fc_data$FoldChange[fc_data$Method == "rawdata"]

fc_data <- fc_data %>%
  mutate(DeltaFC = FoldChange - rawdata_value) %>%
  filter(Method != "rawdata")

fc_data <- fc_data %>%
  arrange(desc(DeltaFC))

fc_data$Method <- factor(fc_data$Method, levels = fc_data$Method)

method_colors <- c(
  "CellClear" = "#FDCDE5",
  "SoupX" = "#4293C6",
  "CellBender" = "#88419D",
  "FastCAR" = "#AE017E",
  "scCDC" = "#F768A1",
  "scAR" = "#0052A1",
  "DecontX" = "#9ECAE1"
)

method_colors <- method_colors[levels(fc_data$Method)]

y_min <- min(fc_data$DeltaFC) * 1.1
y_max <- max(fc_data$DeltaFC) * 1.1

fc_plot <- ggplot(fc_data, aes(x = Method, y = DeltaFC, fill = Method)) +
  geom_col(width = 0.7) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "solid") +
  
  scale_fill_manual(values = method_colors) +
  labs(x = "", y = "Δ Fold Change") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
                               color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.x = element_blank(),
    axis.title.y = element_text(color = "black", size = 12, margin = margin(r = 10)),
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  scale_y_continuous(limits = c(y_min, y_max),
                     breaks = seq(floor(y_min*10)/10, ceiling(y_max*10)/10, by = 0.1),
                     expand = expansion(mult = c(0.1, 0.1)))

print(fc_plot)

pdf("DeltaFC_barplot.pdf", width = 4, height = 6)
print(fc_plot)
dev.off()

png("DeltaFC_barplot.png", width = 1600, height = 2400, res = 300)
print(fc_plot)
dev.off()
