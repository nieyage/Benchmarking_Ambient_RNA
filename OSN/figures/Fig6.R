
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(patchwork)

df <- read_excel("/data/R04/zhangchao/joint_nuclei/figures_new/summary.xlsx")

# The smaller, the better 
reverse_cols <- c(
  "PBMC RMSLE","PBMC sMAPE","PBMC MedAE",
  "non-OSN mean RMSLE", "non-OSN RMSLE variance", "non-OSN Kendall's τ variance", "non-OSN decontamination efficacy variance",
  "non-OSN mean RMSLE (cell number)", "non-OSN mean RMSLE (depth)",
  "OSN mean RMSLE"
)

# index classification
pbmc_cols <- c(
  "PBMC RMSLE","PBMC sMAPE","PBMC MedAE","PBMC Kendall's τ","PBMC R2","PBMC CCC",
  "PBMC prediction score","PBMC average silhouette width","PBMC purity","PBMC k-NN overlap"
)

nonosn_cols <- c(
  "non-OSN mean RMSLE", "non-OSN RMSLE variance",
  "non-OSN mean Kendall's τ", "non-OSN Kendall's τ variance",
  "non-OSN mean decontamination efficacy", "non-OSN decontamination efficacy variance", "non-OSN decontamination efficacy AUC",
  
  "non-OSN mean RMSLE (cell number)", "non-OSN mean RMSLE AUC (cell number)",
  "non-OSN mean Kendall's τ (cell number)", "non-OSN mean Kendall's τ AUC (cell number)",
  "non-OSN decontamination efficacy (cell number)", "non-OSN decontamination efficacy AUC (cell number)",

  "non-OSN mean RMSLE (depth)", "non-OSN mean RMSLE AUC (depth)",
  "non-OSN mean Kendall's τ (depth)", "non-OSN mean Kendall's τ AUC (depth)",
  "non-OSN decontamination efficacy (depth)", "non-OSN decontamination efficacy AUC (depth)"
)

osn_cols <- c(
  "OSN mean RMSLE","OSN mean Kendall's τ","OSN mean decontamination efficacy (cell)",
  "OSN mean decontamination efficacy (reads)","OSN FC of transcriptome distance"
)

# Lengthen format & Normalization & ranking
df_long <- df %>%
  pivot_longer(-method, names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  mutate(
    scaled_value = ifelse(
      metric %in% reverse_cols,
      (max(value, na.rm=TRUE) - value) / (max(value, na.rm=TRUE) - min(value, na.rm=TRUE)),
      (value - min(value, na.rm=TRUE)) / (max(value, na.rm=TRUE) - min(value, na.rm=TRUE))
    ),
    rank_value = rank(-scaled_value, ties.method = "min")
  ) %>%
  ungroup() %>%
  mutate(category = case_when(
    metric %in% pbmc_cols ~ "Accuracy",
    metric %in% nonosn_cols ~ "Robustness",
    metric %in% osn_cols ~ "Sensitivity",
    TRUE ~ "Other"
  ))

# color, rank from high to low
pbmc_colors <- (c("#ffffff","#a8dde1","#75b5dc","#478ecc","#326db6","#2c4ca0","#313772"))
nonosn_colors <- (c("#ffffff","#cfeadf","#a4cbb7","#81b095","#669877","#4d7e54","#376439"))
osn_colors <- (c("#ffffff","#fee3ce","#eabaa1","#dc917b","#d16d5b","#c44438","#b7282e"))
rank_colors <- (c("#ffffff", "#f0f0f0", "#e0e0e0", "#d0d0d0", "#c0c0c0", "#b0b0b0", "#a0a0a0"))

#====================================================================================================
# overall rank is calculated separately for each category
#====================================================================================================

category_rank <- df_long %>%
  filter(category %in% c("Accuracy", "Robustness", "Sensitivity")) %>%
  group_by(category, method) %>%
  summarize(rank1_count = sum(rank_value == 1, na.rm = TRUE), .groups = "drop") %>%
  group_by(category) %>%
  arrange(category, desc(rank1_count)) %>%
  mutate(overall_rank = row_number()) %>%
  ungroup()

# Generate the Rank column data for each category itself
rank_data_cat <- category_rank %>%
  mutate(metric = "Overall Rank") %>%
  select(method, category, metric, overall_rank)

# Merge back into the main data frame
df_long_with_rank <- bind_rows(
  df_long,
  rank_data_cat %>% mutate(
    scaled_value = NA, 
    rank_value = overall_rank, 
    value = overall_rank
  )
)

df_pbmc   <- df_long_with_rank %>% filter(category == "Accuracy")
df_nonosn <- df_long_with_rank %>% filter(category == "Robustness")
df_osn    <- df_long_with_rank %>% filter(category == "Sensitivity")

# Sort each panel separately by method
# Internal sorting of Accuracy
pbmc_order <- df_pbmc %>%
  filter(metric == "Overall Rank") %>%
  arrange(overall_rank) %>%
  pull(method)

# Internal sorting of Robustness
nonosn_order <- df_nonosn %>%
  filter(metric == "Overall Rank") %>%
  arrange(overall_rank) %>%
  pull(method)

# Internal sorting of Sensitivity 
osn_order <- df_osn %>%
  filter(metric == "Overall Rank") %>%
  arrange(overall_rank) %>%
  pull(method)

df_pbmc$method   <- factor(df_pbmc$method, levels = rev(pbmc_order))
df_nonosn$method <- factor(df_nonosn$method, levels = rev(nonosn_order))
df_osn$method    <- factor(df_osn$method, levels = rev(osn_order))


plot_bubble_with_rank <- function(df_sub, cols, color_scale, title, outname, w=12) {
  # color of title background
  strip_bg_color <- ifelse(grepl("Accuracy", title, ignore.case = TRUE), tail(pbmc_colors, 1),  # If it is "Accuracy", use the last color of the pbmc
                        ifelse(grepl("Robustness", title, ignore.case = TRUE), tail(nonosn_colors, 1), 
                        ifelse(grepl("Sensitivity", title, ignore.case = TRUE), tail(osn_colors, 1), 
                        "#808080"))) # Default gray

  # calculate bar length and color
  rank_df <- df_sub %>% filter(metric=="Overall Rank") %>% 
    mutate(
      bar_len = max(overall_rank) - overall_rank + 1, # Inverted ranking
      col_index = length(color_scale) - round((bar_len-min(bar_len))/(max(bar_len)-min(bar_len))*(length(color_scale)-1)),
      bar_col = color_scale[col_index],
      OverallRankTitle = "Overall" # Add a new column for the left facet title
    )

  # Extract non-rank indicators for bubble charts
  bubble_df <- df_sub %>% filter(metric!="Overall Rank")
  bubble_df$metric <- factor(bubble_df$metric, levels = cols) # Convert the indicators to factors and ensure that the bubbles are arranged in the specified order

  # Create a data frame with alternating gray and white backgrounds for each row
  method_levels <- unique(rank_df$method)
  bg_df <- data.frame(
    method = method_levels,
    ymin = as.numeric(factor(method_levels)) - 0.5, # The lower boundary of the background rectangle
    ymax = as.numeric(factor(method_levels)) + 0.5, # The upper boundary of the background rectangle
    fill_bg = rep(c("grey95", "white"), length.out = length(method_levels)) # Alternating gray and white
  )

  p1 <- ggplot() +
    # Draw an alternating gray and white background
    geom_rect(data=bg_df, inherit.aes = FALSE,
              aes(xmin=-Inf, xmax=Inf, ymin=ymin, ymax=ymax),
              fill=bg_df$fill_bg) + # Alternate fill colors
    # Draw bar
    geom_col(data=rank_df,
             aes(x=bar_len, y=method, fill=bar_col),
             width=0.6) + # bar width
    scale_fill_identity() +  # Use the assigned colors
    facet_grid(. ~ OverallRankTitle, scales="free_x", space="free_x") + # Title bar on the left
    theme_minimal() +
    theme(
      panel.grid = element_blank(), # Remove the grid lines
      axis.text.x = element_blank(), # Hide the text on the X-axis
      axis.text.y = element_text(size=11), # Y-axis text size
      axis.title = element_blank(), # Remove the axis title
      strip.background = element_rect(fill=strip_bg_color, color=NA),
      strip.text = element_text(color="white", size=11, face="bold")
    )

  p2 <- ggplot() +
    # Draw a background rectangle with alternating gray and white (consistent with the left side)
    geom_rect(data=bg_df, inherit.aes = FALSE,
              aes(xmin=-Inf, xmax=Inf, ymin=ymin, ymax=ymax),
              fill=bg_df$fill_bg) +
    # draw bubble plot
    geom_point(data=bubble_df, 
               aes(x=metric, y=method, size=scaled_value, fill=rank_value),
               shape=21, color="grey50", stroke=0.2) +
    # Set the color gradient (bubble-filled color)
    scale_fill_gradientn(
      colors=color_scale, 
      name=paste0("Rank (", title, ")") # Legend name, including the title
    ) + 
    scale_size_continuous(range=c(3,10)) + # bubble size
    facet_grid(. ~ category, scales="free_x", space="free_x") + # Create a title bar (strip) and display the value of this category as the title in the title bar
    theme_minimal() +
    theme(
      panel.grid = element_blank(), # Remove the grid lines
      axis.text.x = element_text(angle=45, hjust=1),  # The text on the X-axis is tilted at 45 degrees
      axis.text.y = element_blank(), # Hide the Y-axis text (already shown in the left image)
      axis.title = element_blank(), # Remove the axis title
      strip.background = element_rect(
        fill=strip_bg_color, # The background color of the title bar
        color=NA
      ), 
      strip.text = element_text(color="white", size=11, face="bold") # The text style of the title bar
    )
  # widths=c(2,10) indicates that the left bar chart takes up 2 parts of width and the right bubble chart takes up 10 parts of width
  ggsave(outname, plot=p1 + p2 + plot_layout(widths=c(2,10)), width=w, height=6, dpi=300)
}


# plot
# Accuracy
plot_bubble_with_rank(df_sub = df_pbmc, cols = pbmc_cols, color_scale = pbmc_colors,
  title = "Accuracy", outname = "/data/R04/zhangchao/joint_nuclei/figures_new/accuracy_ranked.pdf",w = 8)
# Robustness
plot_bubble_with_rank(df_sub = df_nonosn, cols = nonosn_cols, color_scale = nonosn_colors,
  title = "Robustness", outname = "/data/R04/zhangchao/joint_nuclei/figures_new/robustness_ranked.pdf", w = 12)
# Sensitivity
plot_bubble_with_rank(df_sub = df_osn, cols = osn_cols, color_scale = osn_colors, 
  title = "Sensitivity",outname = "/data/R04/zhangchao/joint_nuclei/figures_new/sensitivity_ranked.pdf", w = 6)


