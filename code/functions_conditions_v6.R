###############################################################################
## functions_conditions_v6.R – Conditional & shock scenarios for the BVAR   ##
###############################################################################

## This file provides three building blocks used by the main workflow:
##
##  (1)  build_level_path_from_yoy()        – convert a user-supplied YoY %
##                                           path into a quarterly level path.
##  (2)  apply_conditions_flexible()        – map user conditions (rates, levels,
##                                           YoY growth) onto the Y-matrix that
##                                           feeds the Kalman smoother.
##  (3)  conditional_flexible()             – run a fully conditional forecast.
##  (4)  shock_after_baseline_partial_pin() – add ex-post “shocks” on top of any
##                                           baseline (unconditional or
##                                           conditional) and re-estimate.
##
## The code assumes **transformed** model space (log*100 or raw*100) internally.
## All necessary back-transformations are handled locally.

# ---------------------------------------------------------------------------
# 1. build_level_path_from_yoy()
# ---------------------------------------------------------------------------
# Convert a YoY growth profile (in **plain percent**, not log differences) into
# an explicit quarterly level path.  The anchor for each observation is either
# (a) the corresponding quarter one year earlier (for t ≤ 4) or
# (b) the level computed four steps back (for t > 4).
#
# Args
#   last_4_quarters_levels : numeric(4) – most-recent four observed levels.
#   yoy_vector             : numeric(h) – user-specified YoY % growth path.
#   var_label              : character  – used only for error messages.
#
# Returns
#   numeric(h) – level path suitable for insertion into the model state.
# ---------------------------------------------------------------------------

build_level_path_from_yoy <- function(last_4_quarters_levels,
                                      yoy_vector,
                                      var_label = "??") {
  
  if (length(last_4_quarters_levels) < 4)
    stop("Need ≥ 4 historical data points for YoY calculation")
  
  h <- length(yoy_vector)
  levels_out <- rep(NA_real_, h)
  
  for (t in seq_len(h)) {
    g <- yoy_vector[t]                       # plain % growth (e.g. 2.5 means 2.5 %)
    if (is.na(g)) next
    
    anchor <- if (t <= 4) last_4_quarters_levels[t] else levels_out[t - 4]
    if (is.na(anchor)) next
    
    levels_out[t] <- anchor * (1 + g / 100)
  }
  
  levels_out
}

# ---------------------------------------------------------------------------
# 2. apply_conditions_flexible()
# ---------------------------------------------------------------------------
# Overwrite the future rows of the Y-matrix with user-specified paths.
#
# Supported modes in `condition_specs`:
#   * "rate_raw"  – interest-rate style series supplied in percent. Stored as
#                   value*100 to match the model’s raw*100 convention.
#   * "level_log" – series supplied in natural units. Stored as log(x)*100.
#   * "yoy_log"   – YoY % growth path. Converted to levels by
#                   build_level_path_from_yoy(), then stored as log(level)*100.
#
# `condition_specs` structure:
#   list(
#     "Variable name in var_names0" = list(
#         mode = <character>,
#         path = numeric(h)
#     ),
#     ...
#   )
#
# Args
#   Y, end_ind  : internal state matrix and index of last historical row.
#   var_names0  : character vector of original variable names.
#   df_history_raw : data.frame of raw-level history (for YoY anchoring).
#   horizon     : forecast horizon (sanity check).
#
# Returns
#   Updated Y with selected future entries replaced.
# ---------------------------------------------------------------------------

apply_conditions_flexible <- function(Y, end_ind, var_names0, condition_specs,
                                      df_history_raw, horizon) {
  
  for (nm in names(condition_specs)) {
    idx <- which(var_names0 == nm)
    if (length(idx) == 0) {
      message("Variable '", nm, "' not found in var_names0; skipping.")
      next
    }
    
    spec       <- condition_specs[[nm]]
    user_path  <- spec$path
    mode       <- spec$mode
    
    if (length(user_path) > horizon)
      stop("Path for '", nm, "' exceeds forecast horizon")
    
    col_fixed <- make.names(nm)
    if (!col_fixed %in% names(df_history_raw))
      stop("Column '", col_fixed, "' not found in df_history_raw")
    
    if (mode == "rate_raw") {
      # Path supplied in percent; store as x*100.
      for (t in seq_along(user_path))
        if (!is.na(user_path[t]))
          Y[end_ind + t, idx] <- user_path[t] * 100
      
    } else if (mode == "level_log") {
      # Path supplied in levels; store as log(level)*100.
      for (t in seq_along(user_path))
        if (!is.na(user_path[t]))
          Y[end_ind + t, idx] <- log(user_path[t]) * 100
      
    } else if (mode == "yoy_log") {
      # YoY % path –> level path –> log(level)*100.
      hist_col   <- df_history_raw[[col_fixed]]
      last_4     <- tail(hist_col, 4)
      lvl_path   <- build_level_path_from_yoy(last_4, user_path, nm)
      
      for (t in seq_along(lvl_path))
        if (!is.na(lvl_path[t]))
          Y[end_ind + t, idx] <- log(lvl_path[t]) * 100
      
    } else {
      message("Unknown mode '", mode, "' for variable '", nm, "'")
    }
  }
  
  Y
}

# ---------------------------------------------------------------------------
# 3. conditional_flexible()
# ---------------------------------------------------------------------------
# Produce a fully conditional forecast: first pin the future according to
# `condition_specs`, then run draw-by-draw Kalman smoothing to obtain the
# predictive distribution.
#
# Returns
#   3-D array:  [h × n_vars × n_sim]
# ---------------------------------------------------------------------------

conditional_flexible <- function(data, df, df_history_raw, res, p, h, n_sim,
                                 var_names0,
                                 scenario_name  = "cond_flexible",
                                 condition_specs = list(),
                                 verbose = TRUE) {
  
  Y       <- rbind(as.matrix(data), matrix(NA, h, ncol(data)))
  end_ind <- nrow(data)
  
  # Apply user pins
  Y_cond  <- apply_conditions_flexible(
    Y, end_ind, var_names0, condition_specs,
    df_history_raw, horizon = h)
  
  # Draw-by-draw forecast
  forecasts <- array(0, c(h, ncol(Y), n_sim))
  if (verbose) pb <- txtProgressBar(min = 1, max = n_sim, style = 3)
  
  for (i in seq_len(n_sim)) {
    beta  <- res$beta[i, , ]
    sigma <- res$sigma[i, , ]
    
    forecasts[,, i] <- generate_forecasts_one_batch(
      Y_cond, df, beta, sigma, p, h,
      var_names0, end_ind, c_init = beta[1, ])
    if (verbose) setTxtProgressBar(pb, i)
  }
  if (verbose) cat("\nconditional_flexible completed.\n")
  
  forecasts
}

# ---------------------------------------------------------------------------
# 4. shock_after_baseline_partial_pin()
# ---------------------------------------------------------------------------
# Construct a scenario that layers deterministic “shocks” on top of a chosen
# baseline and re-estimates the system, partially pinning only the shocked
# variables.
#
#   • If `baseline_specs` is NULL → start from the unconditional baseline.
#   • Otherwise                         → use a conditional baseline defined by
#                                         `baseline_specs`.
#
# `shock_list` structure:
#   list(
#     "Variable name in var_names0" = list(
#         quarter_start = <int>,
#         quarter_end   = <int>,
#         shift         = <numeric>   # same units as transformed Y
#     ),
#     ...
#   )
#
# Returns
#   3-D array: [h × n_vars × n_sim]
# ---------------------------------------------------------------------------

shock_after_baseline_partial_pin <- function(data, df, df_history_raw,
                                             res, p, h, n_sim,
                                             var_names0,
                                             scenario_name  = "shock_after_baseline",
                                             baseline_specs = NULL,
                                             shock_list     = NULL,
                                             verbose        = TRUE) {
  
  # -- Step 1 : obtain the baseline forecast ---------------------------------
  baseline_fcst <- if (!is.null(baseline_specs)) {
    conditional_flexible(data, df, df_history_raw, res, p, h, n_sim,
                         var_names0,
                         scenario_name  = paste0(scenario_name, "_baseline"),
                         condition_specs= baseline_specs,
                         verbose        = verbose)
  } else {
    generate_unconditional(data, df, res, p, h, n_sim, var_names0)
  }
  
  baseline_median <- apply(baseline_fcst, c(1, 2), median)
  
  # -- Step 2 : inject deterministic shocks ----------------------------------
  if (!is.null(shock_list)) {
    for (nm in names(shock_list)) {
      idx <- which(var_names0 == nm)
      if (length(idx) == 0) next
      
      st  <- shock_list[[nm]]$quarter_start
      ed  <- shock_list[[nm]]$quarter_end
      dv  <- shock_list[[nm]]$shift
      
      for (t_i in seq(st, ed))
        baseline_median[t_i, idx] <- baseline_median[t_i, idx] + dv
    }
  }
  
  # -- Step 3 : partially pin only shocked variables -------------------------
  Y       <- rbind(as.matrix(data), matrix(NA, h, ncol(data)))
  end_ind <- nrow(data)
  
  shock_vars_idx <- if (!is.null(shock_list))
    match(names(shock_list), var_names0)
  else integer(0)
  
  for (t in seq_len(h))
    for (v in shock_vars_idx)
      Y[end_ind + t, v] <- baseline_median[t, v]
  
  # -- Step 4 : re-estimate with partial pins --------------------------------
  scenario_fcst <- array(0, c(h, ncol(Y), n_sim))
  if (verbose) pb <- txtProgressBar(min = 1, max = n_sim, style = 3)
  
  for (i in seq_len(n_sim)) {
    beta  <- res$beta[i, , ]
    sigma <- res$sigma[i, , ]
    
    scenario_fcst[,, i] <- generate_forecasts_one_batch(
      Y, df, beta, sigma, p, h,
      var_names0, end_ind, c_init = beta[1, ])
    if (verbose) setTxtProgressBar(pb, i)
  }
  if (verbose) cat("\nshock_after_baseline_partial_pin completed.\n")
  
  scenario_fcst
}
