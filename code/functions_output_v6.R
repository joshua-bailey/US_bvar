###############################################################################
## functions_output_v6.R – Visualisation & CSV-export helpers                ##
###############################################################################
## This file assembles model outputs into tidy data frames, applies
## last-mile transformations (YoY growth, level back-transforms, etc.),
## and produces ggplot2 graphics suitable for publication.
##
## Dependencies: zoo, ggplot2, dplyr, tidyr, readr, purrr, scales
###############################################################################

library(zoo)         # explicit load so downstream scripts need not remember
library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

# ---------------------------------------------------------------------------
# 1.  calc_yoy()  – year-over-year % change                                  #
# ---------------------------------------------------------------------------
#   x      : numeric vector of levels
#   lag_n  : spacing in periods (defaults to 4 quarters)
#
#   Returns a vector of the same length as x, with NA for the first lag_n
#   observations.
calc_yoy <- function(x, lag_n = 4) {
  out <- rep(NA_real_, length(x))
  for (i in seq(lag_n + 1, length(x)))
    out[i] <- (x[i] / x[i - lag_n] - 1) * 100
  out
}

# ---------------------------------------------------------------------------
# 2.  gather_scenarios_for_plot()                                            #
# ---------------------------------------------------------------------------
#   scenario_list : named list of data frames produced by
#                   combine_fcst_with_history_raw().
#   df_original   : transformed historical df (used only for last actual date)
#   yoy_vars      : character vector to be converted to YoY growth
#   raw_vars      : character vector plotted in levels
#
#   Output: long tidy tibble ready for facet-wrap plotting.
gather_scenarios_for_plot <- function(
    scenario_list,
    df_original,
    yoy_vars = c("GDP", "PCE.Price.Index", "Core.PCE.Price.Index"),
    raw_vars = c("Unemployment.Rate", "Fed.Funds.Rate", "TR10y.Rate")) {
  
  df_original$yq  <- as.yearqtr(df_original$yq)
  last_actual     <- max(df_original$yq, na.rm = TRUE)
  
  window_start <- last_actual - 4   # 4 years back
  window_end   <- last_actual + 3   # 3 years forward
  
  key_vars <- c(yoy_vars, raw_vars)
  out_all  <- NULL
  
  for (s_name in names(scenario_list)) {
    
    df_s <- scenario_list[[s_name]] %>%
      mutate(yq = as.yearqtr(yq)) %>%
      select(any_of(c("yq", key_vars))) %>%
      arrange(yq)
    
    ## YoY conversion – done in raw units so plotting remains intuitive ------
    for (vv in yoy_vars)
      if (vv %in% names(df_s))
        df_s[[vv]] <- calc_yoy(as.numeric(df_s[[vv]]), 4)
    
    ## Ensure raw vars are numeric ------------------------------------------
    for (vv in raw_vars)
      if (vv %in% names(df_s))
        df_s[[vv]] <- as.numeric(df_s[[vv]])
    
    df_s <- df_s %>%
      mutate(scenario = s_name,
             phase    = if_else(yq <= last_actual, "History", "Forecast")) %>%
      filter(yq >= window_start, yq <= window_end)
    
    ## Pad any missing columns (keeps bind_rows tidy) -----------------------
    col_order     <- c("yq", key_vars, "scenario", "phase")
    missing_cols  <- setdiff(col_order, names(df_s))
    for (mc in missing_cols)
      df_s[[mc]] <- if (mc == "yq") as.yearqtr(NA) else NA
    
    df_s <- df_s[, col_order]
    out_all <- bind_rows(out_all, df_s)
  }
  
  ## Long format, factor ordering ------------------------------------------
  df_long <- out_all %>%
    pivot_longer(all_of(key_vars),
                 names_to  = "variable",
                 values_to = "value") %>%
    filter(!is.na(value))
  
  df_long$variable <- factor(df_long$variable,
                             levels = c("GDP", "PCE.Price.Index", "Core.PCE.Price.Index",
                                        "Unemployment.Rate", "Fed.Funds.Rate", "TR10y.Rate"))
  
  scenario_levels <- names(scenario_list)
  if ("Unconditional" %in% scenario_levels)
    scenario_levels <- c("Unconditional",
                         setdiff(scenario_levels, "Unconditional"))
  df_long$scenario <- factor(df_long$scenario, levels = scenario_levels)
  
  df_long
}

# ---------------------------------------------------------------------------
# 3.  plot_key_vars_multi_scenarios()                                        #
# ---------------------------------------------------------------------------
#   Produces a facet-wrapped line chart comparing scenarios for a handful
#   of headline variables.  Historical observations for non-baseline
#   scenarios are suppressed to avoid over-plotting.
plot_key_vars_multi_scenarios <- function(
    df_long,
    main_title = "US BVAR scenarios") {
  
  if ("Unconditional" %in% unique(df_long$scenario))
    df_long <- df_long %>%
      filter(scenario == "Unconditional" | phase == "Forecast")
  
  ## Colour palette – reserve a fixed blue for the baseline  ----------------
  scen <- unique(df_long$scenario)
  pal  <- setNames(rep(NA, length(scen)), scen)
  if ("Unconditional" %in% scen)
    pal["Unconditional"] <- "#1f78b4"
  others <- setdiff(scen, "Unconditional")
  if (length(others) > 0)
    pal[others] <- hue_pal()(length(others))
  
  pretty_labels <- c(
    "GDP"                  = "real gdp yoy",
    "PCE.Price.Index"      = "pce yoy",
    "Core.PCE.Price.Index" = "core pce yoy",
    "Unemployment.Rate"    = "unemployment rate",
    "Fed.Funds.Rate"       = "fed funds rate",
    "TR10y.Rate"           = "10-year treasury yield")
  
  ggplot(df_long,
         aes(yq, value,
             colour   = scenario,
             linetype = phase,
             group    = interaction(scenario, phase))) +
    geom_line(size = 0.75) +
    facet_wrap(~ variable, scales = "free_y",
               labeller = labeller(variable = pretty_labels)) +
    scale_color_manual(values = pal) +
    labs(x = NULL, y = NULL,
         colour = "scenario", linetype = "phase",
         title  = main_title) +
    theme_bw() +
    theme(strip.text  = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---------------------------------------------------------------------------
# 4.  make_fanchart_one_var_export()                                         #
# ---------------------------------------------------------------------------
#   Draws a fan chart (median + 66 / 90 % intervals) for a selected series,
#   back-transforms to real-world units, and optionally writes a CSV.
make_fanchart_one_var_export <- function(
    fcst_array,
    df,
    var_of_interest,
    var_names,
    var_names0,
    trans,
    start_plot = 2015,
    center     = c("median", "mean"),
    export_csv_name = NULL) {
  
  center <- match.arg(center)
  
  idx <- match(var_of_interest, var_names0)
  if (is.na(idx))
    stop("Variable not found in var_names0: ", var_of_interest)
  
  h     <- dim(fcst_array)[1]
  n_sim <- dim(fcst_array)[3]
  mat   <- matrix(fcst_array[, idx, ], h, n_sim)
  
  ## Undo model transforms ---------------------------------------------------
  if (trans[idx] == "log100") mat <- exp(mat / 100)
  if (trans[idx] == "raw")    mat <- mat / 100
  if (var_names0[idx] == "Potential GDP") mat <- mat / 10
  
  probs  <- c(0.05, 0.17, 0.50, 0.83, 0.95)
  fq     <- t(apply(mat, 1, quantile, probs = probs, na.rm = TRUE))
  colnames(fq) <- c("lo90","lo66","med","up66","up90")
  center_line <- if (center == "mean") rowMeans(mat) else fq[, "med"]
  
  ## Historical series (also back-transformed) ------------------------------
  df_hist <- df[, c("yq", make.names(var_of_interest))]
  names(df_hist) <- c("yq", "hist")
  if (trans[idx] == "log100") df_hist$hist <- exp(df_hist$hist / 100)
  if (trans[idx] == "raw")    df_hist$hist <- df_hist$hist / 100
  if (var_names0[idx] == "Potential GDP") df_hist$hist <- df_hist$hist / 10
  df_hist <- df_hist[df_hist$yq >= start_plot, ]
  
  ## Combine with forecasts --------------------------------------------------
  fut_dates <- seq(tail(df_hist$yq, 1) + 0.25,
                   by = 0.25, length.out = h)
  df_fcst <- data.frame(
    yq   = fut_dates,
    lo90 = fq[, "lo90"],
    lo66 = fq[, "lo66"],
    med  = center_line,
    up66 = fq[, "up66"],
    up90 = fq[, "up90"])
  
  ## Optional CSV export -----------------------------------------------------
  if (!is.null(export_csv_name)) {
    out_csv <- bind_rows(
      mutate(df_hist,
             lo90 = NA, lo66 = NA, med = NA, up66 = NA, up90 = NA),
      mutate(df_fcst, hist = NA))
    write_csv(out_csv, export_csv_name)
  }
  
  ## Plot --------------------------------------------------------------------
  ggplot() +
    geom_line(data = df_hist,
              aes(yq, hist, colour = "history"), size = 0.8) +
    geom_ribbon(data = df_fcst,
                aes(yq, ymin = lo90, ymax = up90, fill = "90% interval"),
                alpha = 0.20) +
    geom_ribbon(data = df_fcst,
                aes(yq, ymin = lo66, ymax = up66, fill = "66% interval"),
                alpha = 0.40) +
    geom_line(data = df_fcst,
              aes(yq, med, colour = "forecast"), size = 1) +
    scale_fill_manual(name = "",
                      values = c("90% interval" = "#2c7fb8",
                                 "66% interval" = "#41b6c4")) +
    scale_colour_manual(name = "",
                        values = c("history"  = "black",
                                   "forecast" = "blue")) +
    theme_bw() +
    labs(x = NULL, y = NULL,
         title = paste("fan chart –", var_of_interest)) +
    theme(legend.position = "bottom")
}
