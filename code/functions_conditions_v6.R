########################
## functions_conditions_v6.R
########################


#-------------------------------------------
# 1) build_level_path_from_yoy
#-------------------------------------------
build_level_path_from_yoy <- function(last_4_quarters_levels, yoy_vector, var_label = "??") {
  # this function helps if we want yoy constraints
  # we have the last 4 quarters of real levels, plus a yoy growth path
  # e.g. yoy_vector[t] is the yoy growth in percent
  # we then build the level path by compounding from the relevant quarter
  #
  # returns a vector of levels with length = length(yoy_vector)
  
  if (length(last_4_quarters_levels) < 4) {
    stop("need >=4 historical data points for yoy calculation")
  }
  h <- length(yoy_vector)
  levels_out <- rep(NA_real_, h)
  
  for (t in seq_len(h)) {
    yoy_rate <- yoy_vector[t]
    if (is.na(yoy_rate)) {
      levels_out[t] <- NA
    } else {
      if (t <= 4) {
        # for t=1..4, we anchor on the last_4_quarters_levels
        lvl_4_prior <- last_4_quarters_levels[t]
      } else {
        # beyond that, we anchor on the newly computed levels
        lvl_4_prior <- levels_out[t - 4]
      }
      if (is.na(lvl_4_prior)) {
        levels_out[t] <- NA
      } else {
        # yoy_rate in percent => yoy_rate / 100
        levels_out[t] <- lvl_4_prior * exp(yoy_rate / 100)
      }
    }
  }
  
  return(levels_out)
}

#-------------------------------------------
# 2) apply_conditions_flexible
#-------------------------------------------
apply_conditions_flexible <- function(Y, end_ind, var_names0, condition_specs,
                                      df_history_raw, horizon) {
  # this function pins certain future values in y based on the user's condition specs
  # e.g. yoy constraints or raw interest rates
  #
  # condition_specs is a list: variable_name -> (mode, path)
  #   mode can be "rate_raw", "level_log", or "yoy_log"
  #   path is a numeric vector for each forecast step
  #
  # we directly overwrite Y[end_ind+t, var_index] with the transformed path
  
  for (nm in names(condition_specs)) {
    idx <- which(var_names0 == nm)
    if (length(idx) == 0) {
      message("warning: variable '", nm, "' not found in var_names0")
      next
    }
    spec      <- condition_specs[[nm]]
    user_path <- spec$path
    mode      <- spec$mode
    
    # check horizon
    if (length(user_path) > horizon) {
      stop("user-specified path is longer than forecast horizon for ", nm)
    }
    
    # figure out the column in df_history_raw that we might use for yoy logic
    col_fixed <- make.names(nm)
    if (! col_fixed %in% names(df_history_raw)) {
      stop("column '", col_fixed, "' not found in df_history_raw. check naming.")
    }
    
    if (mode == "rate_raw") {
      # user_path[t] = e.g. 3.50 means 3.50 percent
      # we store it in Y as 350 if the raw transform is x*100
      for (t in seq_along(user_path)) {
        if (!is.na(user_path[t])) {
          Y[end_ind + t, idx] <- user_path[t] * 100
        }
      }
    } else if (mode == "level_log") {
      # user_path[t] is a real level, so we store log(x)*100
      for (t in seq_along(user_path)) {
        if (!is.na(user_path[t])) {
          Y[end_ind + t, idx] <- log(user_path[t]) * 100
        }
      }
    } else if (mode == "yoy_log") {
      # yoy -> build levels from yoy, then store them in log form
      hist_col <- df_history_raw[[col_fixed]]
      last_4   <- tail(hist_col, 4)
      yoy_levels <- build_level_path_from_yoy(last_4, user_path, var_label = nm)
      for (t in seq_along(yoy_levels)) {
        if (!is.na(yoy_levels[t])) {
          Y[end_ind + t, idx] <- log(yoy_levels[t]) * 100
        }
      }
    } else {
      message("unknown mode: ", mode, " for variable: ", nm)
    }
  }
  return(Y)
}

#-------------------------------------------
# 3) conditional_flexible
#-------------------------------------------
conditional_flexible <- function(data, df, df_history_raw, res, p, h, n_sim,
                                 var_names0, scenario_name = "cond_flexible",
                                 condition_specs = list(),
                                 verbose = TRUE) {
  # this function applies conditions to the forecast horizon,
  # then re-runs the bvar draws to incorporate those constraints.
  #
  # data is the full historical matrix
  # df is the historical data frame (transformed)
  # df_history_raw is the historical data frame reverted to real-world levels
  # res is the bvar result
  # condition_specs is a named list describing how each variable is pinned
  #
  # returns a 3d array [h x n_vars x n_sim]
  
  Y <- rbind(as.matrix(data), matrix(NA, h, ncol(data)))
  end_ind <- nrow(data)
  
  # pin the future values as specified
  Y_cond <- apply_conditions_flexible(Y, end_ind, var_names0, condition_specs,
                                      df_history_raw, horizon = h)
  
  # we produce draws
  forecasts <- array(0, c(h, ncol(Y), n_sim))
  if (verbose) {
    pb <- txtProgressBar(min = 1, max = n_sim, style = 3)
  }
  
  for (i in seq_len(n_sim)) {
    beta  <- res$beta[i, , ]
    sigma <- res$sigma[i, , ]
    # run the existing forecast routine
    fcst_i <- generate_forecasts_one_batch(
      Y_cond, df, beta, sigma, p, h,
      var_names0, end_ind, c_init = beta[1,],
      base0 = 0  # no special offset
    )
    forecasts[,, i] <- fcst_i
    
    if (verbose) setTxtProgressBar(pb, i)
  }
  if (verbose) {
    cat("\nconditional_flexible forecast completed.\n")
  }
  
  return(forecasts)
}

#-------------------------------------------
# 4) shock_after_baseline_partial_pin
#-------------------------------------------
shock_after_baseline_partial_pin <- function(data, df, df_history_raw,
                                             res, p, h, n_sim,
                                             var_names0,
                                             scenario_name  = "shock_after_baseline",
                                             baseline_specs = NULL,
                                             shock_list     = NULL,
                                             verbose        = TRUE) {
  # this function:
  # 1) generates a baseline forecast (unconditional or conditional)
  # 2) picks the median path of that baseline
  # 3) adds the shock to the pinned variable(s) for the chosen horizon
  # 4) partially pins those shocked variables, leaving others as NA
  # 5) re-runs the bvar forecast
  #
  # returns a 3d array [h x n_vars x n_sim]
  
  # step 1: baseline scenario
  if (!is.null(baseline_specs)) {
    # conditional scenario
    baseline_fcst <- conditional_flexible(
      data           = data,
      df             = df,
      df_history_raw = df_history_raw,
      res            = res,
      p              = p,
      h              = h,
      n_sim          = n_sim,
      var_names0     = var_names0,
      scenario_name  = paste0(scenario_name, "_baseline"),
      condition_specs= baseline_specs,
      verbose        = verbose
    )
  } else {
    # unconditional
    baseline_fcst <- generate_unconditional(
      data        = data,
      df          = df,
      res         = res,
      p           = p,
      h           = h,
      n_sim       = n_sim,
      var_names0  = var_names0
    )
  }
  
  # get the median path
  baseline_median <- apply(baseline_fcst, c(1, 2), median)
  
  # step 2: add shocks
  if (!is.null(shock_list)) {
    for (nm in names(shock_list)) {
      idx <- which(var_names0 == nm)
      if (length(idx) == 0) next
      
      st        <- shock_list[[nm]]$quarter_start
      ed        <- shock_list[[nm]]$quarter_end
      shift_val <- shock_list[[nm]]$shift
      
      # apply shift
      for (t_i in seq(st, ed)) {
        baseline_median[t_i, idx] <- baseline_median[t_i, idx] + shift_val
      }
    }
  }
  
  # step 3: partial pin
  Y <- rbind(as.matrix(data), matrix(NA, h, ncol(data)))
  end_ind <- nrow(data)
  
  # figure out which variables got shocked
  shock_vars_idx <- integer(0)
  if (!is.null(shock_list)) {
    for (nm in names(shock_list)) {
      idx <- which(var_names0 == nm)
      shock_vars_idx <- c(shock_vars_idx, idx)
    }
  }
  
  # fill in the pinned variable(s) with the new baseline+shift
  for (t in seq_len(h)) {
    for (v in shock_vars_idx) {
      Y[end_ind + t, v] <- baseline_median[t, v]
    }
  }
  
  # step 4: re-run the bvar scenario
  scenario_fcst <- array(0, c(h, ncol(Y), n_sim))
  if (verbose) {
    pb <- txtProgressBar(min = 1, max = n_sim, style = 3)
  }
  for (i in seq_len(n_sim)) {
    beta  <- res$beta[i, , ]
    sigma <- res$sigma[i, , ]
    fcst_i<- generate_forecasts_one_batch(
      Y, df, beta, sigma, p, h,
      var_names0, end_ind, c_init = beta[1,],
      base0 = 0
    )
    scenario_fcst[,, i] <- fcst_i
    if (verbose) setTxtProgressBar(pb, i)
  }
  if (verbose) {
    cat("\nshock_after_baseline_partial_pin scenario completed.\n")
  }
  
  return(scenario_fcst)
}
