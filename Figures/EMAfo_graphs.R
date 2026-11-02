library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(readr)
library(purrr)

# =========================
# Paths & Labels
# =========================

results_dir <- "Z:/Projects/EMA_Project/Scripts/Output/FO_graph"

# Outcome label used in titles/tables = name; code used in file names = value
outcomes <- c(
  "NA" = "negAff",
  "PA" = "posAff",
  "ST" = "stress"
)

strategies <- c(
  "RP" = "Reappraisal",
  "SP" = "Suppression",
  "DS" = "Distraction",
  "SA" = "Selective Attention",
  "SS" = "Situation Selection"
)

types <- c("joint", "single")

# Y-axis range
ylims_joint  <- c(-0.35, 0.5)
ylims_single <- c(-0.35, 0.5)

# Add an "n points outside y-axis range" subtitle to any panel that drops points
annotate_out_of_range <- TRUE

# Read every CSV once

read_with_labels <- function(outcome, strategy, type) {
  fp <- file.path(results_dir, paste0("FO_", outcomes[[outcome]], "_", strategy, "_", type, ".csv"))
  df <- read.csv(fp, stringsAsFactors = FALSE)
  df$Specification <- paste(df$Filter, df$Analysis, df$Covariates, sep = "_")
  df$Sig      <- df$CI.Lower > 0 | df$CI.Upper < 0
  df$Outcome  <- outcome
  df$Strategy <- strategy
  df$Type     <- type
  df
}

file_grid <- expand.grid(
  outcome  = names(outcomes),
  strategy = names(strategies),
  type     = types,
  stringsAsFactors = FALSE
)

all_data <- pmap_dfr(file_grid, read_with_labels)

# Specification parsing (vectorised)

parse_specification <- function(spec) {
  tibble(
    More_Than_10_Responses = as.integer(grepl("More_Than_10_Responses|MT10R_AND_MT20S", spec)),
    More_Than_20_Seconds   = as.integer(grepl("More_Than_20_Seconds|MT10R_AND_MT20S", spec)),
    Start_Date             = as.integer(grepl("Start_Date|Date_and_Time", spec)),
    Time_of_Day            = as.integer(grepl("Time_of_Day|Date_and_Time", spec)),
    Demographics           = as.integer(grepl("Demographics", spec) & !grepl("No_Demographics", spec)),
    Baseline_ER            = as.integer(grepl("_Baseline_ER", spec) & !grepl("No_Baseline_ER", spec)),
    Prospective            = as.integer(grepl("Prospective", spec))
  )
}

spec_levels <- c(
  "Prospective", ">10_Resp", ">20_Sec", "Start_Date",
  "Time_of_Day", "Demographics", "Baseline_ER"
)

spec_colors <- c(
  "Prospective" = "purple3", ">10_Resp" = "firebrick", ">20_Sec" = "firebrick",
  "Start_Date" = "sienna", "Time_of_Day" = "sienna",
  "Demographics" = "sienna", "Baseline_ER" = "sienna"
)

# =========================
# Plot builders
# =========================

create_coefficient_plot <- function(data, title, ylims) {
  data$Color <- ifelse(data$Sig, "green3", "tomato")
  n <- nrow(data)
  
  in_range <- !is.na(data$Coefficient) &
    data$Coefficient >= ylims[1] & data$Coefficient <= ylims[2]
  n_out <- sum(!is.na(data$Coefficient) & !in_range)
  
  subtitle <- if (annotate_out_of_range && n_out > 0) {
    paste0(n_out, " of ", n, " point", ifelse(n_out == 1, "", "s"),
           " outside y-axis range (not shown)")
  } else NULL
  
  ggplot(data[in_range, ], aes(x = Index, y = Coefficient)) +
    geom_point(aes(color = Color), size = 4, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 2) +
    scale_color_identity() +
    scale_y_continuous(limits = ylims) +
    scale_x_continuous(limits = c(0.5, n + 0.5)) +   # keeps alignment with the spec heatmap
    labs(title = title, subtitle = subtitle, x = "", y = "Beta Coefficient") +
    theme_minimal(base_size = 6) +
    theme(
      plot.title    = element_text(size = 44, hjust = 0.2, face = "bold"),
      plot.subtitle = element_text(size = 28, hjust = 0.2, face = "italic"),
      axis.text.x   = element_text(size = 34, face = "bold"),
      axis.text.y   = element_text(size = 34, face = "bold"),
      axis.title.y  = element_text(size = 34, face = "bold"),
      axis.title.x  = element_text(size = 1),
      legend.position = "none"
    )
}

create_specification_plot <- function(data) {
  n <- nrow(data)
  
  long_data <- parse_specification(data$Specification) %>%
    mutate(Index = data$Index) %>%
    pivot_longer(-Index, names_to = "Variable", values_to = "Value") %>%
    mutate(
      Variable = recode(Variable,
                        More_Than_10_Responses = ">10_Resp",
                        More_Than_20_Seconds   = ">20_Sec"),
      Variable = factor(Variable, levels = spec_levels),
      Fill     = ifelse(Value == 1L, as.character(Variable), NA_character_)
    )
  
  ggplot(long_data, aes(x = Index, y = Variable, fill = Fill)) +
    geom_tile(height = 0.8) +
    scale_fill_manual(values = spec_colors, na.value = "white", guide = "none") +
    scale_x_continuous(limits = c(0.5, n + 0.5)) +
    labs(x = "Specification", y = "Analytic Dimension") +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.x  = element_text(size = 30, face = "bold"),
      axis.text.y  = element_text(size = 30, face = "bold"),
      axis.title.x = element_text(size = 34, face = "bold"),
      axis.title.y = element_text(size = 34, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "none"
    )
}

# One panel = coefficient plot over spec heatmap, sorted by coefficient
make_panel <- function(data, title, ylims) {
  data <- data %>%
    arrange(Coefficient) %>%
    mutate(Index = row_number())
  create_coefficient_plot(data, title, ylims) /
    create_specification_plot(data) +
    plot_layout(heights = c(2, 1))
}

# All 15 panels for one model type: rows = NA / PA / ST, columns = strategies
make_patchwork <- function(type, ylims) {
  panels <- list()
  for (o in names(outcomes)) {
    for (s in names(strategies)) {
      d <- all_data %>% filter(Outcome == o, Strategy == s, Type == type)
      panels[[length(panels) + 1]] <- make_panel(d, paste0(o, " and ", strategies[[s]]), ylims)
    }
  }
  panels
}

# Patchwork
plots_joint <- make_patchwork("joint", ylims_joint)

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

# Patchwork
plots_single <- make_patchwork("single", ylims_single)

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

results_summary <- all_data %>%
  group_by(Outcome, Type, Strategy) %>%
  summarise(
    Total = n(),
    Sig   = sum(Sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Block    = ifelse(Type == "joint", "Joint", "Single"),
    Outcome  = factor(Outcome,  levels = names(outcomes)),
    Strategy = factor(Strategy, levels = names(strategies)),
    Percent  = round(100 * Sig / Total, 3)
  ) %>%
  select(Outcome, Block, Strategy, Total, Sig, Percent) %>%
  arrange(Block, Outcome, Strategy)

print(results_summary, n = Inf)

get_significance_summary <- function(outcome) {
  wide <- all_data %>%
    filter(Outcome == outcome, Type == "joint") %>%
    distinct(Specification, Strategy, .keep_all = TRUE) %>%
    select(Specification, Strategy, Sig) %>%
    pivot_wider(id_cols = Specification, names_from = Strategy, values_from = Sig,
                values_fill = FALSE)
  
  present <- intersect(names(strategies), names(wide))
  wide <- wide %>%
    mutate(
      Outcome = outcome,
      Num_Significant_Strategies = rowSums(across(all_of(present)), na.rm = TRUE)
    )
  
  total_specs <- nrow(wide)
  
  distribution <- wide %>%
    count(Outcome, Num_Significant_Strategies, name = "Num_Specifications") %>%
    complete(
      Outcome = outcome,
      Num_Significant_Strategies = 0:5,
      fill = list(Num_Specifications = 0)
    ) %>%
    mutate(Percent = round(100 * Num_Specifications / total_specs, 3)) %>%
    arrange(Outcome, Num_Significant_Strategies)
  
  list(
    by_spec = wide %>% arrange(Specification),
    distribution = distribution
  )
}

na_res <- get_significance_summary("NA")
pa_res <- get_significance_summary("PA")
st_res <- get_significance_summary("ST")

# Combined distribution table
signif_strategy_distribution <- bind_rows(
  na_res$distribution, pa_res$distribution, st_res$distribution
)

print(signif_strategy_distribution, n = Inf)