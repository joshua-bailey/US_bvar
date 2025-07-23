########################
## functions_data_v6.R
########################


# set your FRED API key
fredr_set_key("3f66fc645f49a16dba8369e2539076f8")
set.seed(1234)

#-----------------------------------
# helper: convert monthly data to quarterly by averaging
#-----------------------------------
fetch_monthly_to_quarterly <- function(series_id, start_date) {
  # this grabs monthly data from fred and then groups it by quarter
  df_m <- fredr(
    series_id         = series_id,
    observation_start = as.Date(paste0(floor(start_date), "-01-01"))
  ) %>%
    mutate(
      yq  = as.yearqtr(date),
      val = value
    ) %>%
    select(yq, val)
  
  df_q <- df_m %>%
    group_by(yq) %>%
    summarise(value = mean(val, na.rm = TRUE)) %>%
    ungroup()
  
  return(df_q)
}

#-----------------------------------
# helper: fetch direct quarterly data
#-----------------------------------
fetch_quarterly <- function(series_id, start_date) {
  # this grabs quarterly data directly from fred
  df_q <- fredr(
    series_id         = series_id,
    observation_start = as.Date(paste0(floor(start_date), "-01-01")),
    frequency         = "q",
    aggregation_method = "avg"
  ) %>%
    mutate(
      yq = as.yearqtr(date),
      value = value
    ) %>%
    select(yq, value) %>%
    distinct(yq, .keep_all = TRUE)
  
  return(df_q)
}

#-----------------------------------
# main function to fetch data from fred and apply transformations
#-----------------------------------
fetch_data_from_fred <- function(var_config, start_date = 1959, end_date = 2022.50) {
  # var_config should have columns:
  #   var_names0 (descriptive name)
  #   fred_id    (the fred series id)
  #   freq       (either 'monthly' or 'quarterly')
  #   transform  (either 'log100' or 'raw')
  #
  # we pull each series from fred, rename columns, merge by yq, and filter by date range
  # then we apply transformations: log100 => log(x)*100, raw => x*100
  # we return a list with:
  #   df        -> the data frame with columns yq + the variables
  #   data      -> a numeric matrix for the bvar
  #   trans     -> the transform vector
  #   var_names -> safe column names in r
  #   var_names0-> original descriptive variable names
  
  decimal_to_yearqtr <- function(x) {
    # converts something like 2022.50 to a yearqtr => 2022 q3
    year_part <- floor(x)
    frac <- x - year_part
    quarter_index <- case_when(
      abs(frac - 0.00) < 1e-9 ~ 1,
      abs(frac - 0.25) < 1e-9 ~ 2,
      abs(frac - 0.50) < 1e-9 ~ 3,
      abs(frac - 0.75) < 1e-9 ~ 4,
      TRUE ~ 4
    )
    as.yearqtr(paste(year_part, quarter_index, sep = " q"))
  }
  
  start_yq <- decimal_to_yearqtr(start_date)
  end_yq   <- decimal_to_yearqtr(end_date)
  
  series_dfs <- list()
  for (i in seq_len(nrow(var_config))) {
    freq_i  <- var_config$freq[i]
    fred_id <- var_config$fred_id[i]
    
    # fetch monthly or quarterly as specified
    if (freq_i == "monthly") {
      df_temp <- fetch_monthly_to_quarterly(fred_id, start_date)
    } else {
      df_temp <- fetch_quarterly(fred_id, start_date)
    }
    
    col_new <- var_config$var_names0[i]
    df_temp <- df_temp %>% rename(!!col_new := value)
    series_dfs[[i]] <- df_temp
  }
  
  # merge all variables by yq
  df_merged <- reduce(series_dfs, full_join, by = "yq") %>%
    arrange(yq)
  
  # restrict to start_yq:end_yq
  df_filtered <- df_merged %>%
    filter(yq >= start_yq & yq <= end_yq)
  
  var_names0 <- var_config$var_names0
  var_names  <- make.names(var_names0)
  
  # rename columns from descriptive to safe r names
  for (i in seq_along(var_names0)) {
    oldn <- var_names0[i]
    newn <- var_names[i]
    if (oldn %in% names(df_filtered)) {
      names(df_filtered)[names(df_filtered) == oldn] <- newn
    }
  }
  
  # apply transformations
  trans <- var_config$transform
  for (i in seq_along(var_names)) {
    col_i <- var_names[i]
    if (trans[i] == "log100") {
      df_filtered[[col_i]] <- log(df_filtered[[col_i]]) * 100
    } else if (trans[i] == "raw") {
      df_filtered[[col_i]] <- df_filtered[[col_i]] * 100
    }
  }
  
  data_matrix <- as.matrix(df_filtered[, var_names, drop = FALSE])
  
  return(list(
    df        = df_filtered,
    data      = data_matrix,
    trans     = trans,
    var_names = var_names,
    var_names0= var_names0
  ))
}

# -------------------------------------------------------------------
#  estimate_bvar  – Minnesota‑prior BVAR with optional SOC prior
#                   (fully aligned with CBO implementation)
# -------------------------------------------------------------------
estimate_bvar <- function(
    data,
    p,
    n_draw,
    n_burn,
    verbose = FALSE,
    use_soc = FALSE          # <- toggle for sum‑of‑coefficients dummy
){
  
  ## 1. prepare Minnesota prior (λ = 0.06, α = 2, ψ = OLS residual s.d.)
  mode_psi <- calculate_psi(data, p)
  mn0 <- bv_mn(
    bv_lambda(mode = 0.06),
    bv_alpha(mode = 2),
    bv_psi(mode = mode_psi)   # fixes prior scale to sample residuals
  )
  
  ## 2. optional SOC prior
  if (use_soc) {
    soc0   <- bv_soc(mode = 0.20, sd = 1, min = 1e-4, max = 50)
    priors <- bv_priors(hyper = c("lambda", "soc"), mn = mn0, soc = soc0)
    mh0    <- bv_metropolis(scale_hess = c(0.01, 0.01))  # TWO params
  } else {
    priors <- bv_priors(hyper = c("lambda"), mn = mn0)
    mh0    <- bv_metropolis(scale_hess = c(0.01))        # ONE param
  }
  
  ## 3. run MCMC
  res <- bvar(
    data    = data,
    lags    = p,
    n_draw  = n_draw,
    n_burn  = n_burn,
    n_thin  = 1L,
    priors  = priors,
    mh      = mh0,
    fcast   = NULL,
    irf     = NULL,
    verbose = verbose
  )
  
  return(res)
}

#-----------------------------------
# calculate_psi: used by estimate_bvar for the initial scale
#-----------------------------------
calculate_psi <- function(data, p) {
  # we do a quick ols for each variable and get the residual stdev
  # that stdev is used as an initial guess for the prior
  psi_vals <- c()
  for (i in seq_len(ncol(data))) {
    y <- data[(p + 1):nrow(data), i]
    x <- matrix(1, nrow = length(y), ncol = 1) # intercept
    for (j in 1:p) {
      x <- cbind(x, data[(p + 1 - j):(nrow(data) - j), i])
    }
    beta_ols <- solve(t(x) %*% x) %*% t(x) %*% y
    e        <- y - x %*% beta_ols
    psi_vals <- c(psi_vals, sd(e))
  }
  return(psi_vals)
}