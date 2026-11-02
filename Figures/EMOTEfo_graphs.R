library(tidyr)
library(dplyr)
library(readr)
library(purrr)

# Paths & Labels

results_dir <- "Z:/Projects/EMA_Project/Scripts/Output/EMOTE_Results"

outcomes <- c(
  "NA" = "Negative Affect",
  "PA" = "Positive Affect"
)

strategies <- c(
  "RP" = "Reappraisal",
  "SP" = "Suppression",
  "DS" = "Distraction"
)

timings        <- c("concurrent", "prospective")
predictor_sets <- c("all", names(strategies))   # "all" = joint model, RP/SP/DS = single models

# Read every CSV once

# Standard (unadjusted) CI significance
compute_ci_sig <- function(beta, se, alpha = 0.05) {
  ok    <- is.finite(beta) & is.finite(se) & se > 0
  z     <- qnorm(1 - alpha / 2)
  lower <- ifelse(ok, beta - z * se, NA_real_)
  upper <- ifelse(ok, beta + z * se, NA_real_)
  sig   <- ok & ((lower > 0) | (upper < 0))
  tibble(sig = sig, lower = lower, upper = upper)
}

read_with_labels <- function(timing, outcome, set) {
  fp <- file.path(results_dir, paste0(timing, "_CLUSTERED_", outcome, "_", set, ".csv"))
  df <- read_csv(fp, show_col_types = FALSE)
  df$Source_File <- fp
  df$Outcome[is.na(df$Outcome)] <- "NA"   # read_csv reads the string "NA" as missing
  df <- df %>%
    mutate(
      Analysis      = trimws(Analysis),
      Predictor_Set = trimws(Predictor_Set),
      Outcome       = trimws(Outcome),
      Predictor     = trimws(Predictor)
    ) %>%
    filter(!tolower(Predictor) %in% c("intercept"))
  bind_cols(df, compute_ci_sig(df$Coefficient, df$Std_Err, alpha = 0.05))
}

file_grid <- expand.grid(
  timing  = timings,
  outcome = names(outcomes),
  set     = predictor_sets,
  stringsAsFactors = FALSE
)

all_data <- pmap_dfr(file_grid, read_with_labels) %>%
  filter(Predictor %in% names(strategies), Outcome %in% names(outcomes)) %>%
  mutate(
    Strategy = Predictor,
    Type     = if_else(Predictor_Set == "All", "joint", "single"),
    Sig      = sig
  )

# TABLES: one per model type (Joint / Independent)
# Rows = Outcome x Concurrent; columns = N (obs.) + one beta per strategy
# * = 95% CI excludes 0

table_long <- all_data %>%
  mutate(
    Model         = factor(if_else(Type == "joint", "Joint", "Independent"),
                           levels = c("Joint", "Independent")),
    Outcome_Label = factor(outcomes[Outcome], levels = outcomes),
    Concurrent    = factor(if_else(Analysis == "Concurrent", "\u2714", "\u2718"),
                           levels = c("\u2714", "\u2718")),
    IV            = factor(strategies[Strategy], levels = strategies),
    Beta          = paste0(sprintf("%.3f", Coefficient), ifelse(Sig, "*", ""))
  )

# Sanity check: each cell should come from exactly one estimate
dupes <- table_long %>% count(Model, Outcome_Label, Concurrent, IV) %>% filter(n > 1)
if (nrow(dupes) > 0) {
  warning("Some cells have more than one estimate:")
  print(dupes)
}

fmt_n <- function(n) format(n, big.mark = ",", trim = TRUE)

for (m in levels(table_long$Model)) {
  d <- table_long %>% filter(Model == m)
  
  # N (obs.) per row. If the strategies in a row were fit on different Ns
  # (possible for independent models), all distinct Ns are listed.
  n_col <- d %>%
    group_by(Outcome_Label, Concurrent) %>%
    summarise(`N (obs.)` = paste(fmt_n(sort(unique(Observations))), collapse = " / "),
              .groups = "drop")
  
  wide <- d %>%
    select(Outcome_Label, Concurrent, IV, Beta) %>%
    pivot_wider(names_from = IV, values_from = Beta) %>%
    left_join(n_col, by = c("Outcome_Label", "Concurrent")) %>%
    arrange(Outcome_Label, Concurrent) %>%
    select(`Outcome Variable` = Outcome_Label, Concurrent, `N (obs.)`,
           any_of(unname(strategies)))
  
  cat("\n=====================================================\n")
  cat("EMOTE First-Order Analysis:", m, "Model\n")
  cat("(* = 95% CI excludes 0)\n")
  cat("=====================================================\n")
  print(as.data.frame(wide), row.names = FALSE, right = FALSE)
}