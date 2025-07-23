########################
## functions_output_v6.R
########################

#--------------------------------------------
# 1) helper: calc_yoy (year-over-year % change)
#--------------------------------------------
calc_yoy <- function(x, lag_n = 4) {
  # this function calculates yoy by looking at x[t] / x[t-lag_n] - 1, in percent
  # if x is, say, real gdp, then yoy is ((gdp now)/(gdp last year)-1)*100
  # default lag is 4 for quarterly yoy
  out <- rep(NA_real_, length(x))
  for (i in seq(lag_n + 1, length(x))) {
    out[i] <- (x[i] / x[i - lag_n] - 1) * 100
  }
  return(out)
}

#--------------------------------------------
# 2) gather_scenarios_for_plot
#--------------------------------------------
gather_scenarios_for_plot <- function(
    scenario_list,
    df_original,     # the historical df so we know the last actual date
    yoy_vars = c("GDP", "PCE.Price.Index", "Core.PCE.Price.Index"),
    raw_vars = c("Unemployment.Rate", "Fed.Funds.Rate", "TR10y.Rate")
) {
  # this function merges multiple scenario data frames (like unconditional, fedOnly, etc.)
  # applies yoy transformations to yoy_vars, keeps raw_vars in their raw form,
  # then combines them in a long format for easy multi-scenario plotting.
  #
  # note that df_original has the earliest date -> we find the last actual date from it
  
  df_original$yq <- as.yearqtr(df_original$yq)
  last_actual_date <- max(df_original$yq, na.rm = TRUE)
  
  # define a date window
  start_date <- last_actual_date - 4  # show ~6 years back
  end_date   <- last_actual_date + 3  # show ~3 years forward
  
  key_vars <- c(yoy_vars, raw_vars)
  out_all  <- NULL
  
  for (s_name in names(scenario_list)) {
    df_scenario <- scenario_list[[s_name]]
    df_scenario$yq <- as.yearqtr(df_scenario$yq)
    
    # pick only yq + the key vars
    df_temp <- df_scenario %>%
      select(any_of(c("yq", key_vars))) %>%
      arrange(yq)
    
    # yoy transformations
    for (vv in yoy_vars) {
      if (vv %in% names(df_temp)) {
        df_temp[[vv]] <- as.numeric(df_temp[[vv]])
        df_temp[[vv]] <- calc_yoy(df_temp[[vv]], 4)
      }
    }
    # ensure raw vars are numeric
    for (vv in raw_vars) {
      if (vv %in% names(df_temp)) {
        df_temp[[vv]] <- as.numeric(df_temp[[vv]])
      }
    }
    
    # label scenario and phase
    df_temp$scenario <- s_name
    df_temp <- df_temp %>%
      mutate(phase = if_else(yq <= last_actual_date, "History", "Forecast"))
    
    # filter date window
    df_temp <- df_temp %>%
      filter(yq >= start_date & yq <= end_date)
    
    # make sure all expected cols exist
    col_order <- c("yq", key_vars, "scenario", "phase")
    missing_cols <- setdiff(col_order, names(df_temp))
    for (mc in missing_cols) {
      if (mc == "yq") {
        df_temp[[mc]] <- as.yearqtr(NA)
      } else if (mc %in% c("scenario", "phase")) {
        df_temp[[mc]] <- as.character(NA)
      } else {
        df_temp[[mc]] <- NA_real_
      }
    }
    df_temp <- df_temp[, col_order, drop = FALSE]
    
    # append
    out_all <- bind_rows(out_all, df_temp)
  }
  
  # pivot long
  df_long <- out_all %>%
    pivot_longer(
      cols = all_of(key_vars),
      names_to = "variable",
      values_to = "value"
    ) %>%
    filter(!is.na(value))
  
  # set factor levels
  df_long$variable <- factor(
    df_long$variable,
    levels = c("GDP", "PCE.Price.Index", "Core.PCE.Price.Index",
               "Unemployment.Rate", "Fed.Funds.Rate", "TR10y.Rate")
  )
  
  # reorder scenario factor so 'Unconditional' if present is first
  scenario_order <- names(scenario_list)
  if ("Unconditional" %in% scenario_order) {
    scenario_order <- c("Unconditional", scenario_order[scenario_order != "Unconditional"])
  }
  df_long$scenario <- factor(df_long$scenario, levels = scenario_order)
  
  return(df_long)
}

#--------------------------------------------
# 3) plot_key_vars_multi_scenarios
#--------------------------------------------
plot_key_vars_multi_scenarios <- function(df_long, main_title = "us bvar scenarios") {
  # this function creates a facet wrap of variables, each with lines for each scenario,
  # but only 'unconditional' scenario shows history while others show forecast.
  #
  # a color is assigned to each scenario, with unconditional possibly singled out.
  
  # drop history rows for non-unconditional scenarios
  if ("Unconditional" %in% unique(df_long$scenario)) {
    df_long <- df_long %>%
      filter(scenario == "Unconditional" | phase == "Forecast")
  }
  
  # define color mapping so unconditional is a fixed color (e.g. blue)
  unique_scenarios <- unique(df_long$scenario)
  scenario_order   <- levels(df_long$scenario)
  color_values     <- setNames(rep(NA, length(unique_scenarios)), unique_scenarios)
  if ("Unconditional" %in% unique_scenarios) {
    color_values["Unconditional"] <- "#1f78b4"  # fixed blue
  }
  other_scenarios <- setdiff(unique_scenarios, "Unconditional")
  if (length(other_scenarios) > 0) {
    other_colors <- hue_pal()(length(other_scenarios))
    names(other_colors) <- other_scenarios
    color_values[other_scenarios] <- other_colors
  }
  
  # rename variables for plotting
  var_labels <- c(
    "GDP"                  = "real gdp yoy",
    "PCE.Price.Index"      = "pce yoy",
    "Core.PCE.Price.Index" = "core pce yoy",
    "Unemployment.Rate"    = "unemployment rate",
    "Fed.Funds.Rate"       = "fed funds rate",
    "TR10y.Rate"           = "10-year treasury yield"
  )
  
  p <- ggplot(df_long,
              aes(x = yq, y = value,
                  colour = scenario,
                  linetype = phase,
                  group = interaction(scenario, phase))) +
    geom_line(size = 0.75) +
    facet_wrap(~ variable, scales = "free_y",
               labeller = labeller(variable = var_labels)) +
    scale_color_manual(values = color_values) +
    labs(x = NULL, y = NULL, colour = "scenario", linetype = "phase",
         title = main_title) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

#--------------------------------------------
# 4) make_fanchart_one_var_export
#--------------------------------------------
make_fanchart_one_var_export <- function(
    fcst_array,
    df,
    var_of_interest,
    var_names,
    var_names0,
    trans,
    start_plot     = 2015,
    center         = "median",
    impose_zero    = FALSE,    # not used by new main script, can remove if truly unneeded
    export_csv_name= NULL
) {
  # this function builds a single-variable fan chart from the forecast array [h x n_vars x n_sim]
  # then merges historical data from df (in transformed form),
  # reverts the transform, and plots ribbons for 5%,17%,83%,95% intervals.
  #
  # if export_csv_name is given, it also writes out a csv with the historical and forecast intervals
  
  idx_var <- which(var_names0 == var_of_interest)
  if (length(idx_var) == 0) {
    stop(paste("variable", var_of_interest, "not found in var_names0"))
  }
  
  h      <- dim(fcst_array)[1]
  n_sim  <- dim(fcst_array)[3]
  fcst_sub <- fcst_array[, idx_var, , drop = FALSE]
  
  # flatten to h x n_sim
  fcst_2d <- matrix(fcst_sub, nrow = h, ncol = n_sim)
  
  # revert transforms
  if (trans[idx_var] == "log100") {
    fcst_2d <- exp(fcst_2d / 100)
  } else if (trans[idx_var] == "raw") {
    fcst_2d <- fcst_2d / 100
  }
  # potential gdp scaling if needed
  if (var_names0[idx_var] == "Potential GDP") {
    fcst_2d <- fcst_2d / 10
  }
  
  quants <- c(0.05, 0.17, 0.50, 0.83, 0.95)
  fcst_quant <- apply(fcst_2d, 1, quantile, probs = quants, na.rm = TRUE)
  fcst_quant <- t(fcst_quant)
  colnames(fcst_quant) <- c("q05","q17","q50","q83","q95")
  
  if (center == "mean") {
    mean_line <- apply(fcst_2d, 1, mean, na.rm = TRUE)
    fcst_quant <- cbind(fcst_quant, center = mean_line)
  } else {
    fcst_quant <- cbind(fcst_quant, center = fcst_quant[,"q50"])
  }
  
  # if impose_zero is never used, we can remove it. but let's keep it in case needed.
  if (impose_zero) {
    fcst_quant <- pmax(fcst_quant, 0)
  }
  
  # revert historical
  df_hist <- df[, c("yq", make.names(var_of_interest))]
  colnames(df_hist) <- c("yq","hist_value")
  
  if (trans[idx_var] == "log100") {
    df_hist$hist_value <- exp(df_hist$hist_value / 100)
  } else if (trans[idx_var] == "raw") {
    df_hist$hist_value <- df_hist$hist_value / 100
  }
  if (var_names0[idx_var] == "Potential GDP") {
    df_hist$hist_value <- df_hist$hist_value / 10
  }
  
  df_hist <- df_hist[df_hist$yq >= start_plot, ]
  
  # future dates
  last_hist  <- tail(df_hist$yq, 1)
  future_dates <- seq(as.yearqtr(last_hist) + 0.25, by = 0.25, length.out = h)
  
  df_fcst <- data.frame(
    yq       = future_dates,
    lo90     = fcst_quant[,"q05"],
    lo66     = fcst_quant[,"q17"],
    median   = fcst_quant[,"center"],
    up66     = fcst_quant[,"q83"],
    up90     = fcst_quant[,"q95"]
  )
  
  # combine historical + forecast for optional csv
  n_hist <- nrow(df_hist)
  df_hist_out <- data.frame(
    yq         = df_hist$yq,
    hist_value = df_hist$hist_value,
    lo90       = rep(NA, n_hist),
    lo66       = rep(NA, n_hist),
    median     = rep(NA, n_hist),
    up66       = rep(NA, n_hist),
    up90       = rep(NA, n_hist)
  )
  n_fcst <- nrow(df_fcst)
  df_fcst_out <- data.frame(
    yq         = df_fcst$yq,
    hist_value = rep(NA, n_fcst),
    lo90       = df_fcst$lo90,
    lo66       = df_fcst$lo66,
    median     = df_fcst$median,
    up66       = df_fcst$up66,
    up90       = df_fcst$up90
  )
  df_plot <- rbind(df_hist_out, df_fcst_out)
  
  if (!is.null(export_csv_name)) {
    write.csv(df_plot, file = export_csv_name, row.names = FALSE)
  }
  
  # build ggplot
  p <- ggplot() +
    geom_line(
      data = df_hist_out,
      aes(x = yq, y = hist_value, colour = "history"),
      size = 0.8
    ) +
    geom_ribbon(
      data = df_fcst_out,
      aes(x = yq, ymin = lo90, ymax = up90, fill = "90% interval"),
      alpha = 0.2
    ) +
    geom_ribbon(
      data = df_fcst_out,
      aes(x = yq, ymin = lo66, ymax = up66, fill = "66% interval"),
      alpha = 0.4
    ) +
    geom_line(
      data = df_fcst_out,
      aes(x = yq, y = median, colour = "forecast"),
      size = 1
    ) +
    scale_fill_manual(name = "", values = c("90% interval" = "#2c7fb8", "66% interval" = "#41b6c4")) +
    scale_colour_manual(name = "", values = c("history" = "black", "forecast" = "blue")) +
    theme_bw() +
    labs(x = NULL, y = NULL, title = paste("fan chart for", var_of_interest)) +
    theme(legend.position = "bottom")
  
  return(p)
}