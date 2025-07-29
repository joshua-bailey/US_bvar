###############################################################################
## functions_data_v6.R – Data ingest & pre-processing for the BVAR pipeline  ##
###############################################################################
## Tasks performed in this file
##   • Download monthly / quarterly series from FRED
##   • Aggregate to quarterly frequency (simple averages)
##   • Apply CBO-style transformations
##   • Provide helper utilities for inverse transforms
##   • Supply a convenience wrapper for BVAR estimation
###############################################################################

library(tidyverse)
library(zoo)     # yearqtr class
library(fredr)

# ---------------------------------------------------------------------------
# 0.  Setup
# ---------------------------------------------------------------------------
# Users are encouraged to store their API key in an environment variable
# (e.g. FRED_API_KEY).  The hard-coded key below is public and rate-limited.
fredr_set_key("3f66fc645f49a16dba8369e2539076f8")
set.seed(1234)   # deterministic behaviour in small helper routines

# ---------------------------------------------------------------------------
# 1.  Utility functions
# ---------------------------------------------------------------------------

#' Convert a decimal “YYYY.q” year to a zoo::yearqtr
#'
#' Example: 2024.25 → "2024 Q2"
decimal_to_yearqtr <- function(x) {
  yr   <- floor(x)
  frac <- x - yr
  q    <- case_when(
    abs(frac - 0.00) < 1e-9 ~ 1,
    abs(frac - 0.25) < 1e-9 ~ 2,
    abs(frac - 0.50) < 1e-9 ~ 3,
    abs(frac - 0.75) < 1e-9 ~ 4,
    TRUE                    ~ 4)
  as.yearqtr(sprintf("%d Q%d", yr, q))
}

#' Download a *monthly* FRED series and convert to quarterly averages
fetch_monthly_to_quarterly <- function(series_id, start_date) {
  fredr(series_id,
        observation_start = as.Date(sprintf("%d-01-01", floor(start_date)))) %>%
    transmute(yq = as.yearqtr(date), val = value) %>%
    group_by(yq) %>%
    summarise(value = mean(val, na.rm = TRUE), .groups = "drop")
}

#' Download a *quarterly* FRED series (average of intra-quarter observations)
fetch_quarterly <- function(series_id, start_date) {
  fredr(series_id,
        observation_start = as.Date(sprintf("%d-01-01", floor(start_date))),
        frequency          = "q",
        aggregation_method = "avg") %>%
    distinct(date, .keep_all = TRUE) %>%
    transmute(yq = as.yearqtr(date), value)
}

# ---------------------------------------------------------------------------
# 2.  Master download + transform
# ---------------------------------------------------------------------------

#' Fetch a set of FRED series and apply model-ready transformations
#'
#' @param var_config  Tibble with columns: var_names0, fred_id, freq, transform
#'                    • var_names0 : descriptive names (used in plots / tables)
#'                    • fred_id    : FRED series ID
#'                    • freq       : "monthly" or "quarterly"
#'                    • transform  : "log100" (log(x)*100) or "raw" (x*100)
#' @param start_date  First observation in decimal-year notation (e.g. 1986)
#' @param end_date    Last observation in decimal-year notation (e.g. 2025.75)
#'
#' @return list(
#'           df         = wide data.frame (quarterly, transformed),
#'           data       = numeric matrix for BVAR,
#'           trans      = character vector of transform codes,
#'           var_names  = safe column names (make.names),
#'           var_names0 = original descriptive names
#'         )
fetch_data_from_fred <- function(var_config,
                                 start_date = 1959,
                                 end_date   = 2022.50) {
  
  stopifnot(all(c("var_names0","fred_id","freq","transform") %in%
                  names(var_config)))
  
  start_q <- decimal_to_yearqtr(start_date)
  end_q   <- decimal_to_yearqtr(end_date)
  
  # --- 2.1  Download each series -------------------------------------------
  lst <- vector("list", nrow(var_config))
  for (i in seq_len(nrow(var_config))) {
    tmp <- if (var_config$freq[i] == "monthly")
      fetch_monthly_to_quarterly(var_config$fred_id[i], start_date)
    else
      fetch_quarterly(var_config$fred_id[i],    start_date)
    lst[[i]] <- rename(tmp, !!var_config$var_names0[i] := value)
  }
  
  # --- 2.2  Merge & clip to requested sample -------------------------------
  df <- reduce(lst, full_join, by = "yq") %>%
    arrange(yq) %>%
    filter(yq >= start_q, yq <= end_q)
  
  # Use syntactically safe column names inside the model
  safe_names <- make.names(var_config$var_names0)
  names(df)[-1] <- safe_names
  
  # --- 2.3  Apply CBO transforms -------------------------------------------
  for (j in seq_along(safe_names)) {
    v <- safe_names[j]
    if (var_config$transform[j] == "log100") {
      df[[v]] <- log(pmax(df[[v]], .Machine$double.eps)) * 100
    } else {
      df[[v]] <- df[[v]] * 100
    }
  }
  
  list(
    df        = df,
    data      = as.matrix(df[, safe_names, drop = FALSE]),
    trans     = var_config$transform,
    var_names = safe_names,
    var_names0= var_config$var_names0
  )
}

#' Reverse the model transforms back to real-world levels
#'
#' Useful for human-readable outputs and YoY anchoring.
rebuild_raw_levels <- function(df_trans, trans, var_names0) {
  df_raw <- df_trans
  for (j in seq_along(trans)) {
    v <- names(df_raw)[j + 1]      # skip yq column
    if (trans[j] == "log100") df_raw[[v]] <- exp(df_raw[[v]] / 100)
    if (trans[j] == "raw")    df_raw[[v]] <- df_raw[[v]] / 100
    if (var_names0[j] == "Potential GDP")
      df_raw[[v]] <- df_raw[[v]] / 10      # CBO scaling quirk
  }
  df_raw
}

# ---------------------------------------------------------------------------
# 3.  BVAR estimation helpers (unchanged logic, documented here)
# ---------------------------------------------------------------------------

#' Compute variable-specific Psi hyper-parameters for the Minnesota prior
calculate_psi <- function(data, p) {
  sapply(seq_len(ncol(data)), function(i) {
    y   <- data[(p + 1):nrow(data), i]
    lag <- embed(data[, i], p + 1)[, -1]
    sd(y - lag %*% solve(crossprod(lag)) %*% crossprod(lag, y))
  })
}

#' Wrapper around bvar::bvar() with Minnesota (+ optional SOC) prior
#'
#' @param data     Numeric matrix (T × n)
#' @param p        Lag order
#' @param n_draw   Total MCMC draws
#' @param n_burn   Burn-in draws
#' @param use_soc  Logical – include stochastic volatility in coefficients?
estimate_bvar <- function(data, p, n_draw, n_burn,
                          verbose = FALSE, use_soc = FALSE) {
  
  mn <- bv_mn(bv_lambda(mode = 0.1),           # overall tightness
              bv_alpha (mode = 2),             # cross-eqn tightness
              bv_psi   (mode = calculate_psi(data, p)))  # variable-specific
  
  pri <- if (use_soc)
    bv_priors(hyper = c("lambda","soc"),
              mn    = mn,
              soc   = bv_soc(mode = 0.20, sd = 1,
                             min = 1e-4, max = 1))
  else
    bv_priors(hyper = "lambda", mn = mn)
  
  bvar(data,
       lags    = p,
       n_draw  = n_draw,
       n_burn  = n_burn,
       priors  = pri,
       mh      = bv_metropolis(scale_hess = 0.01),
       verbose = verbose)
}
