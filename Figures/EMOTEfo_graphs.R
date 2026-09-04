# ---- Packages ----
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(readr)
library(stringr)
library(tools)

# ---- Inputs: all result files you listed ----
all_files <- c(
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\concurrent_CLUSTERED_PA_SP.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\concurrent_CLUSTERED_NA_all.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\concurrent_CLUSTERED_NA_DS.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\concurrent_CLUSTERED_NA_RP.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\concurrent_CLUSTERED_NA_SP.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\concurrent_CLUSTERED_PA_all.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\concurrent_CLUSTERED_PA_DS.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\concurrent_CLUSTERED_PA_RP.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\prospective_CLUSTERED_PA_all.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\prospective_CLUSTERED_PA_DS.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\prospective_CLUSTERED_PA_RP.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\prospective_CLUSTERED_PA_SP.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\prospective_CLUSTERED_NA_all.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\prospective_CLUSTERED_NA_DS.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\prospective_CLUSTERED_NA_RP.csv",
  "Z:\\Projects\\EMA_Project\\Scripts\\Output\\EMOTE_Results\\prospective_CLUSTERED_NA_SP.csv"
)

# ---- Adaptive CI significance ----
compute_ci_sig <- function(beta, se, obs, alpha_base) {
  ok <- is.finite(beta) & is.finite(se) & is.finite(obs) & se > 0 & obs > 0
  alphai <- ifelse(ok, alpha_base / sqrt(obs / 100), NA_real_)
  z      <- ifelse(ok, qnorm(1 - alphai / 2), NA_real_)
  lower  <- beta - z * se
  upper  <- beta + z * se
  sig    <- ok & ((lower > 0) | (upper < 0))
  tibble(sig = sig, lower = lower, upper = upper, alphai = alphai, z = z)
}

# ---- Load, harmonize, and annotate ----
read_one <- function(fp) {
  df <- read_csv(fp, show_col_types = FALSE)
  df$Source_File <- fp
  df$Outcome[is.na(df$Outcome)] <- "NA"
  df <- df %>%
    mutate(
      Analysis       = trimws(Analysis),
      Predictor_Set  = trimws(Predictor_Set),
      Outcome        = trimws(Outcome),
      Predictor      = trimws(Predictor)
    ) %>%
    mutate(
      Spec_JointSingle = if_else(Predictor_Set == "All", "Joint", "Single")
    ) %>%
    filter(!tolower(Predictor) %in% c("intercept"))
  ci <- compute_ci_sig(
    beta = df$Coefficient,
    se   = df$Std_Err,
    obs  = df$Observations,
    alpha_base = 0.05
  )
  bind_cols(df, ci)
}

dat <- bind_rows(lapply(all_files, read_one))
dat <- dat %>%
  filter(Predictor %in% c("RP","SP","DS"),
         Outcome   %in% c("NA","PA")) %>%
  mutate(Color = if_else(sig, "green3", "tomato")) %>%
  group_by(Outcome, Predictor) %>%
  arrange(Coefficient, .by_group = TRUE) %>%
  mutate(Index = row_number()) %>%
  ungroup()

# ---- Specification long-form for tile plot ----
spec_long <- dat %>%
  transmute(
    Outcome, Predictor, Index,
    Concurrent = as.integer(Analysis == "Concurrent"),
    Joint      = as.integer(Spec_JointSingle == "Joint")
  ) %>%
  pivot_longer(cols = c(Concurrent, Joint), names_to = "Variable", values_to = "Value") %>%
  mutate(
    Variable = factor(Variable, levels = c("Concurrent","Joint")),
    Fill = if_else(Value == 1L, as.character(Variable), NA_character_)
  )

color_mapping <- c("Concurrent" = "purple3", "Joint" = "orange2")

# ---- Plotting helpers ----
create_coefficient_plot <- function(df_sub, title) {
  ggplot(df_sub, aes(x = Index, y = Coefficient)) +
    geom_point(aes(color = Color), size = 6, alpha = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
    scale_color_identity() +
    labs(title = title, x = "", y = "Beta Coefficient") +
    theme_minimal(base_size = 6) +
    coord_cartesian(ylim = c(-0.2, 0.4)) +
    theme(
      plot.title   = element_text(size = 19, hjust = 0.2),
      axis.text.x  = element_text(size = 17),
      axis.text.y  = element_text(size = 17),
      axis.title.y = element_text(size = 17),
      axis.title.x = element_text(size = 1),
      legend.position = "none"
    )
}

create_specification_plot <- function(spec_sub) {
  ggplot(spec_sub, aes(x = Index, y = Variable, fill = Fill)) +
    geom_tile(height = 0.5) +
    scale_fill_manual(values = color_mapping, na.value = "white") +
    labs(x = "Specification", y = "Analytic Dimension") +
    theme_minimal(base_size = 2) +
    theme(
      axis.text.x  = element_text(size = 15),
      axis.text.y  = element_text(size = 15),
      axis.title.x = element_text(size = 17),
      axis.title.y = element_text(size = 17),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

# ---- Build per-Outcome panels (RP, SP, DS), each stacking coef + spec ----
preds <- c("RP","SP","DS")
titles_na <- paste0("NA and ", c("Reappraisal","Suppression","Distraction"))
titles_pa <- paste0("PA and ", c("Reappraisal","Suppression","Distraction"))

plots_na <- lapply(seq_along(preds), function(i) {
  dfp   <- dat %>% filter(Outcome == "NA", Predictor == preds[i])
  specp <- spec_long %>% filter(Outcome == "NA", Predictor == preds[i])
  coef_plot <- create_coefficient_plot(dfp, titles_na[i])
  spec_plot <- create_specification_plot(specp)
  coef_plot / spec_plot + plot_layout(heights = c(2, 1))
})

plots_pa <- lapply(seq_along(preds), function(i) {
  dfp   <- dat %>% filter(Outcome == "PA", Predictor == preds[i])
  specp <- spec_long %>% filter(Outcome == "PA", Predictor == preds[i])
  coef_plot <- create_coefficient_plot(dfp, titles_pa[i])
  spec_plot <- create_specification_plot(specp)
  coef_plot / spec_plot + plot_layout(heights = c(2, 1))
})

# ---- Wrap each as a panel ----
final_plot_na <- wrap_plots(plots_na, ncol = 3, guides = "collect") +
  plot_annotation(
    title = "Negative Affect (NA) First-Order Analysis",
    theme = theme(plot.title = element_text(size = 25, hjust = 0.5))
  )

final_plot_pa <- wrap_plots(plots_pa, ncol = 3, guides = "collect") +
  plot_annotation(
    title = "Positive Affect (PA) First-Order Analysis",
    theme = theme(plot.title = element_text(size = 25, hjust = 0.5))
  )

# ---- Combine final NA and PA panels vertically ----
combined_fo_plot <- final_plot_na / final_plot_pa +
  plot_annotation(
    title = "First-Order Analysis",
    theme = theme(plot.title = element_text(size = 30, hjust = 0.5))
  )

# ---- Save the combined plot ----
ggsave(
  filename = "EMOTE_FirstOrder_Combined.png",
  plot = combined_fo_plot,
  width = 16, height = 14, dpi = 300, limitsize = FALSE
)
