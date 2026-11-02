library(ggplot2)
library(tidyr)
library(patchwork)

# SETUP!

# Mappings
title_mapping <- list(
  "PA" = "Positive Affect",
  "NA" = "Negative Affect",
  "ST" = "Stress",
  "RP" = "Reappraisal",
  "SP" = "Suppression",
  "DS" = "Distraction",
  "SA" = "Selective Attention",
  "SS" = "Situation Selection"
)

color_mapping <- list(
  "PA" = "#87C55F",
  "NA" = "#F89C74",
  "ST" = "#F6CF71",
  "RP" = "#9EB9F3",
  "SP" = "#9EB9F3",
  "DS" = "#9EB9F3",
  "SA" = "#9EB9F3",
  "SS" = "#9EB9F3"
)

variables <- c("PA", "NA", "ST", "RP", "SP", "DS", "SA", "SS")

# Load data and compute average transition matrices
base_dir <- "Z:/Projects/EMA_Project/Scripts/Output/Multiverse_Scratch"
average_transition_matrices <- list()

for (variable in variables) {
  transition_matrices_path <- file.path(
    base_dir,
    paste0(variable, "_transition_matrices_No_Filtering_Transition.csv")
  )
  df <- read.csv(transition_matrices_path)
  
  # Average each column (25 columns = flattened 5x5 matrix, first row is columns 1-5, second row is columns 6-10, etc.)
  avg_vec <- colMeans(df, na.rm = TRUE)
  avg_matrix <- matrix(avg_vec, nrow = 5, ncol = 5, byrow = TRUE)
  
  average_transition_matrices[[variable]] <- avg_matrix
}

# Functions

# Reshape a 5x5 matrix into long format for ggplot (slightly unnecessary to go to 5x5 and then back out, but easier to check intermediates)
mat_to_long <- function(mat) {
  df <- as.data.frame(mat)
  colnames(df) <- as.character(1:5)
  # Reverse factor levels for "From" so state 1 plots at the top,
  # matching matplotlib's imshow (row 0 at top, origin = 'upper')
  df$From <- factor(1:5, levels = 5:1)
  
  df_long <- pivot_longer(df, cols = -From, names_to = "To", values_to = "value")
  df_long$To <- factor(df_long$To, levels = as.character(1:5))
  df_long
}

# Build one heatmap panel
plot_transition_heatmap <- function(mat, variable) {
  df_long <- mat_to_long(mat)
  df_long$text_color <- ifelse(df_long$value > 0.5, "white", "black")
  
  var_color <- color_mapping[[variable]]
  
  ggplot(df_long, aes(x = To, y = From, fill = value)) +
    geom_tile(color = "grey90") +
    geom_text(aes(label = sprintf("%.2f", value), color = text_color), size = 3, fontface = "bold") +
    scale_color_identity() +
    scale_fill_gradient(low = "white", high = var_color, limits = c(0, 1), name = NULL) +
    labs(title = title_mapping[[variable]], x = "To State", y = "From State") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(size = 9, face = "bold"),
      axis.text.y = element_text(size = 9, face = "bold"),
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "none",
      aspect.ratio = 1
    )
}

# Section-header labels ("Affect Variables" / "ER Strategies"),
make_section_title <- function(label, size = 16) {
  ggplot() +
    theme_void() +
    labs(title = label) +
    theme(plot.title = element_text(size = size, face = "bold", hjust = 0.5))
}
# Assemble descriptive figure
affect_variables <- c("PA", "NA", "ST")
er_strategies <- c("RP", "SP", "DS", "SA", "SS")
affect_plots <- lapply(affect_variables, function(v) {
  plot_transition_heatmap(average_transition_matrices[[v]], v)
})
er_plots <- lapply(er_strategies, function(v) {
  plot_transition_heatmap(average_transition_matrices[[v]], v)
})
affect_title <- make_section_title("Affect Variables")
er_title <- make_section_title("ER Strategies")

design <- "
AAA
BCD
EEE
FGH
IJ#
"

combined <- wrap_plots(
  A = affect_title,
  B = affect_plots[[1]], C = affect_plots[[2]], D = affect_plots[[3]],
  E = er_title,
  F = er_plots[[1]], G = er_plots[[2]], H = er_plots[[3]],
  I = er_plots[[4]], J = er_plots[[5]],
  design = design,
  heights = c(0.1, 1, 0.1, 1, 1)
) +
  plot_annotation(
    title = "Average Transition Matrices",
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  )

# Save
output_dir <- "Z:/Projects/EMA_Project/Scripts/Figures/EMA"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_path <- file.path(output_dir, "Average_Transition_Matrices.png")
ggsave(output_path, combined, width = 12, height = 13.8, dpi = 300, bg = "white")