library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(readr)

# =========================
# Function Definitions
# =========================

parse_specification <- function(spec) {
  list(
    More_Than_10_Responses = ifelse(grepl("More_Than_10_Responses|MT10R_AND_MT20S", spec), 1, 0),
    More_Than_20_Seconds = ifelse(grepl("More_Than_20_Seconds|MT10R_AND_MT20S", spec), 1, 0),
    Start_Date = ifelse(grepl("Start_Date|Date_and_Time", spec), 1, 0),
    Time_of_Day = ifelse(grepl("Time_of_Day|Date_and_Time", spec), 1, 0),
    Demographics = ifelse(grepl("Demographics", spec) & !grepl("No_Demographics", spec), 1, 0),
    Baseline_ER = ifelse(grepl("_Baseline_ER", spec) & !grepl("No_Baseline_ER", spec), 1, 0),
    Prospective = ifelse(grepl("Prospective", spec), 1, 0)
  )
}

create_coefficient_plot <- function(file_path, title) {
  data <- read.csv(file_path)
  data <- data[order(data$Coefficient), ]
  data$Color <- ifelse(data$CI.Lower > 0 | data$CI.Upper < 0, "green3", "tomato")
  data$Index <- seq_len(nrow(data))
  
  ggplot(data, aes(x = Index, y = Coefficient)) +
    geom_point(aes(color = Color), size = 4, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 2) +
    scale_color_identity() +
    ylim(-0.35, 0.5) +
    labs(
      title = title,
      x = "",
      y = "Beta Coefficient"
    ) +
    theme_minimal(base_size = 6) +
    theme(
      plot.title = element_text(size = 44, hjust = 0.2, face = "bold"),
      axis.text.x = element_text(size = 34, face = "bold"),
      axis.text.y = element_text(size = 34, face = "bold"),
      axis.title.y = element_text(size = 34, face = "bold"),
      axis.title.x = element_text(size = 1),
      legend.position = "none"
    )
}

create_specification_plot <- function(file_path) {
  data <- read.csv(file_path)
  data <- data[order(data$Coefficient), ]
  data$Covariates <- paste(data$Filter, data$Analysis, data$Covariates, sep = "_")
  parsed_data <- do.call(rbind, lapply(data$Covariates, parse_specification))
  parsed_data <- as.data.frame(parsed_data)
  parsed_data$Index <- seq_len(nrow(parsed_data))
  parsed_data <- parsed_data %>% mutate(across(everything(), as.numeric))
  variable_order <- c(
    "Prospective", ">10_Resp", ">20_Sec", "Start_Date", 
    "Time_of_Day", "Demographics", "Baseline_ER"
  )
  long_data <- parsed_data %>%
    pivot_longer(cols = -Index, names_to = "Variable", values_to = "Value") %>%
    mutate(
      Variable = gsub("More_Than_10_Responses", ">10_Resp", Variable),
      Variable = gsub("More_Than_20_Seconds", ">20_Sec", Variable),
      Variable = factor(Variable, levels = variable_order),
      Value = as.factor(Value)
    )
  color_mapping <- c(
    "Prospective" = "purple3", ">10_Resp" = "firebrick", ">20_Sec" = "firebrick",
    "Start_Date" = "sienna", "Time_of_Day" = "sienna", 
    "Demographics" = "sienna", "Baseline_ER" = "sienna"
  )
  long_data <- long_data %>%
    mutate(Fill = ifelse(Value == "1", as.character(Variable), NA))
  
  ggplot(long_data, aes(x = Index, y = Variable, fill = Fill)) +
    geom_tile(height = 0.8) +
    scale_fill_manual(
      values = color_mapping,
      na.value = "white",
      name = "Value"
    ) +
    labs(
      x = "Specification",
      y = "Analytic Dimension"
    ) +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.x = element_text(size = 30, face = "bold"),
      axis.text.y = element_text(size = 30, face = "bold"),
      axis.title.x = element_text(size = 34, face = "bold"),
      axis.title.y = element_text(size = 34, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

# =========================
# Patchwork Creation Function
# =========================

make_patchwork <- function(file_paths, titles) {
  plots <- list()
  for (i in seq_along(file_paths)) {
    coefficient_plot <- create_coefficient_plot(file_paths[i], titles[i])
    specification_plot <- create_specification_plot(file_paths[i])
    plots[[i]] <- coefficient_plot / specification_plot + plot_layout(heights = c(2, 1))
  }
  plots
}

# =========================
# File Paths and Titles
# =========================

titles_na <- c(
  "NA and Reappraisal",
  "NA and Suppression",
  "NA and Distraction",
  "NA and Selective Attention",
  "NA and Situation Selection"
)
titles_pa <- c(
  "PA and Reappraisal",
  "PA and Suppression",
  "PA and Distraction",
  "PA and Selective Attention",
  "PA and Situation Selection"
)
titles_st <- c(
  "ST and Reappraisal",
  "ST and Suppression",
  "ST and Distraction",
  "ST and Selective Attention",
  "ST and Situation Selection"
)

# Joint
paths_na_joint <- c(
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_RP_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_SP_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_DS_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_SA_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_SS_joint.csv"
)
paths_pa_joint <- c(
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_RP_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_SP_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_DS_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_SA_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_SS_joint.csv"
)
paths_st_joint <- c(
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_RP_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_SP_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_DS_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_SA_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_SS_joint.csv"
)

# Single
paths_na_single <- c(
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_RP_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_SP_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_DS_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_SA_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_SS_single.csv"
)
paths_pa_single <- c(
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_RP_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_SP_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_DS_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_SA_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_SS_single.csv"
)
paths_st_single <- c(
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_RP_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_SP_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_DS_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_SA_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_SS_single.csv"
)

# =========================
# Create Joint Patchwork
# =========================

plots_joint <- c(
  make_patchwork(paths_na_joint, titles_na),
  make_patchwork(paths_pa_joint, titles_pa),
  make_patchwork(paths_st_joint, titles_st)
)

final_joint_all <- wrap_plots(plots_joint, ncol = 5, nrow = 3, guides = "collect") +
  plot_annotation(
    title = "First-Order Multiverse Analysis",
    theme = theme(plot.title = element_text(size = 64, hjust = 0.5, face = "bold")),
    tag_levels = "a"
  )

ggsave(
  filename = "FO_joint_multiverse_analysis.png",
  plot = final_joint_all,
  width = 60,
  height = 60,
  dpi = 300,
  limitsize = FALSE
)

# =========================
# Create Single Patchwork
# =========================

plots_single <- c(
  make_patchwork(paths_na_single, titles_na),
  make_patchwork(paths_pa_single, titles_pa),
  make_patchwork(paths_st_single, titles_st)
)

final_single_all <- wrap_plots(plots_single, ncol = 5, nrow = 3, guides = "collect") +
  plot_annotation(
    title = "First-Order Multiverse Analysis (Independent Models)",
    theme = theme(plot.title = element_text(size = 64, hjust = 0.5, face = "bold")),
    tag_levels = "a"
  )

ggsave(
  filename = "FO_single_multiverse_analysis.png",
  plot = final_single_all,
  width = 60,
  height = 60,
  dpi = 300,
  limitsize = FALSE
)
# IN TEXT BELOW

# Function to read file, mark significance, and tag outcome/strategy/type
read_with_labels <- function(path) {
  df <- read.csv(path)
  df$Significant <- ifelse(df$CI.Lower > 0 | df$CI.Upper < 0, 1, 0)
  df$Specification <- paste(df$Filter, df$Analysis, df$Covariates, sep = "_")
  
  # Extract outcome, strategy, type from filename
  m <- regexec("FO_([a-zA-Z]+)_([A-Z]{2})_(joint|single)", basename(path))
  parts <- regmatches(basename(path), m)[[1]]
  
  df$Outcome <- parts[2]   # e.g., "negAff", "posAff", "stress"
  df$Strategy <- parts[3]  # e.g., "RP", "SP", ...
  df$Type <- parts[4]      # joint / single
  
  df
}

# Collect all files (multi/joint and single)
all_files <- c(
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_RP_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_SP_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_DS_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_SA_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_SS_joint.csv",
  
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_RP_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_SP_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_DS_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_SA_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_negAff_SS_single.csv",
  
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_RP_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_SP_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_DS_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_SA_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_SS_joint.csv",
  
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_RP_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_SP_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_DS_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_SA_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_posAff_SS_single.csv",
  
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_RP_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_SP_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_DS_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_SA_joint.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_SS_joint.csv",
  
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_RP_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_SP_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_DS_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_SA_single.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\FO_graph\\FO_stress_SS_single.csv"
)

# Read all files into one dataframe
all_data <- bind_rows(lapply(all_files, read_with_labels))

# ------------------------------
# (1) Outcome-level summary (joint/multi only)
# ------------------------------
outcome_summary <- all_data %>%
  filter(Type == "joint") %>%
  group_by(Outcome, Specification) %>%
  summarise(Num_Significant = sum(Significant), .groups = "drop") %>%
  count(Outcome, Num_Significant, name = "Num_Specifications") %>%
  group_by(Outcome) %>%
  mutate(Percent = round(100 * Num_Specifications / sum(Num_Specifications), 3))

# ------------------------------
# (2) Outcome × Strategy × Type summary
# ------------------------------
strategy_summary <- all_data %>%
  group_by(Outcome, Strategy, Type, Specification) %>%
  summarise(Sig_Spec = any(Significant == 1), .groups = "drop") %>%
  group_by(Outcome, Strategy, Type) %>%
  summarise(
    Percent_Significant = round(100 * mean(Sig_Spec), 3),
    Num_Specifications = n(),
    .groups = "drop"
  )

# Show results
view(outcome_summary)
view(strategy_summary)
