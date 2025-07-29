###############################################################################
## main.R: download data ▸ estimate BVAR ▸ build scenarios
###############################################################################
#
## main.R – Script :                                       ##
##          • download FRED data                                             ##
##          • estimate a Minnesota-prior BVAR (optionally SOC-augmented)     ##
##          • produce unconditional and conditional forecasts                ##
##          • overlay ex-post “shocks”                                       ##
##          • generate comparison plots & fan charts                         ##
#
#  Author:  Original CBO code (Hark Yoo) – ported / extended for FRED by
#           Joshua Bailey. 
###############################################################################

## Required helper files (sourced below):
##   ─ functions_data_v6.R          – data download / transforms
##   ─ functions_forecasting_v6.R   – state-space wrapper + Kalman smoother
##   ─ functions_conditions_v6.R    – conditional & shock utilities
##   ─ functions_output_v6.R        – visualisation / export helpers
###############################################################################

## ── 0.  Libraries ──────────────────────────────────────────────────────────
library(tidyverse)   # data wrangling
library(readr)       # CSV I/O
library(zoo)         # yearqtr helper
library(BVAR)        # posterior sampler
library(fredr)       # FRED API
library(ggplot2)     # plotting
library(scales)      # pretty axes

## ── 1.  Global options ─────────────────────────────────────────────────────
##   Reproducible RNG seed: honour environment variable if present.
seed_env <- Sys.getenv("BVAR_SEED")
set.seed(if (nzchar(seed_env)) as.integer(seed_env) else 1234)

## ── 2.  Source helper scripts ──────────────────────────────────────────────
source("code/functions_data_v6.R")
source("code/functions_forecasting_v6.R")
source("code/functions_conditions_v6.R")
source("code/functions_output_v6.R")

## ── 3.  Variable configuration (FRED ticker-map) ───────────────────────────
##   Each tibble row: var_names0 – human-readable label (used for pinning)
##                    fred_id    – FRED series code
##                    freq       – "monthly"/"quarterly"
##                    transform  – "log100" (log-level ×100) or "raw" (×100)
nipa <- tribble(
  ~var_names0,                                             ~fred_id,  ~freq,       ~transform,
  "GDP",                                                   "GDPC1",   "quarterly", "log100",
  "Consumption",                                           "PCECC96", "quarterly", "log100",
  "Investment (nonres)",                                   "GPDIC1",  "quarterly", "log100",
  "Exports",                                               "EXPGSC1", "quarterly", "log100",
  "Imports",                                               "IMPGSC1", "quarterly", "log100",
  "GDP Deflator",                                          "GDPDEF",  "quarterly", "log100",
  "PCE Price Index",                                       "PCEPI",   "quarterly", "log100",
  "Core PCE Price Index",                                  "PCEPILFE","quarterly", "log100",
  "CPI-U",                                                 "CPIAUCSL","quarterly", "log100",
  "Potential GDP",                                         "GDPPOT",  "quarterly", "log100",
  "Real Govt. Consumption Expenditures & Gross Investment","GCEC1",   "quarterly", "log100"
)

labour <- tribble(
  ~var_names0,         ~fred_id,   ~freq,    ~transform,
  "Payroll Employment","PAYEMS",   "monthly","log100",
  "Labor Force",       "CLF16OV",  "monthly","log100",
  "Wages",             "A4102C1Q027SBEA","quarterly","log100",
  "Unemployment Rate", "UNRATE",   "monthly","raw"
)

financial <- tribble(
  ~var_names0,        ~fred_id, ~freq,    ~transform,
  "Fed Funds Rate",   "FEDFUNDS","monthly","raw",
  "TR3m Rate",        "TB3MS",  "monthly","raw",
  "TR10y Rate",       "GS10",   "monthly","raw",
  "TR1y Rate",        "GS1",    "monthly","raw",
  "Moody's Aaa",      "AAA",    "monthly","raw",
  "Industrial Production","INDPRO","monthly","log100"
)

additional <- tribble(
  ~var_names0,                               ~fred_id,    ~freq,    ~transform,
  "Cleveland Fed 10-yr Inflation Expectations","EXPINF10YR","monthly","raw",
  "Economic Policy Uncertainty Index",        "USEPUINDXM","monthly","log100",
  "WTI Spot Oil Price",                       "MCOILWTICO","monthly","log100"
)

var_config <- bind_rows(nipa, labour, financial, additional)

## ── 4.  Data download & transform ──────────────────────────────────────────
start_date <- 1986                 # inclusive
end_date   <- 2025.00              # decimal year (2025Q4)

dfs <- fetch_data_from_fred(var_config, start_date, end_date)
df          <- dfs$df          # transformed levels (tibble)
data_mat    <- dfs$data        # numeric matrix (T × N)
trans       <- dfs$trans
var_names   <- dfs$var_names   # safe column names
var_names0  <- dfs$var_names0  # original labels

## Convenience: raw-level historical df (for YoY constraints)
df_raw <- rebuild_raw_levels(df, trans, var_names0)

## ── 5.  BVAR estimation ────────────────────────────────────────────────────
p       <- 4        # lags
n_draw  <- 2000L    # posterior draws
n_burn  <- 1000L    # burn-in
n_sim   <- n_draw - n_burn

res <- estimate_bvar(
  data_mat, p,
  n_draw = n_draw, n_burn = n_burn,
  verbose = TRUE,
  use_soc = TRUE)        # SOC hyper-prior enabled

## ── 6.  Unconditional forecast ─────────────────────────────────────────────
h <- 41                           # quarters ahead

forecasts          <- list()
forecasts$base     <- generate_unconditional(
  data_mat, df, res, p, h, n_sim, var_names0)
uncond_fcst        <- forecasts$base

df_combined_raw <- combine_fcst_with_history_raw(
  uncond_fcst, df, h,
  trans, var_names, var_names0,
  center = "median")
write_csv(df_combined_raw, "unconditional_forecast.csv")

## ── 7.  Conditional forecasts – illustrative examples ----------------------
h_cond <- 12   # shorter horizon for narrative scenarios

### Example A – GDP & CPI YoY path ------------------------------------------
cond_gdp_cpi <- list(
  "GDP"   = list(mode = "yoy_log",
                 path = c(-0.16, -1.20, -1.56, -0.92,
                          0.99,  1.77,  2.05, rep(NA, 5))),
  "CPI-U" = list(mode = "yoy_log",
                 path = c(NA, 4.0, 3.5, rep(NA, 9)))
)

forecasts$cond_gdp_cpi <- conditional_flexible(
  data           = data_mat,
  df             = df,
  df_history_raw = df_raw,
  res            = res,
  p              = p,
  h              = h_cond,
  n_sim          = n_sim,
  var_names0     = var_names0,
  scenario_name  = "cond_gdp_cpi",
  condition_specs= cond_gdp_cpi,
  verbose        = TRUE)

df_cond_gdp_cpi <- combine_fcst_with_history_raw(
  forecasts$cond_gdp_cpi, df, h_cond,
  trans, var_names, var_names0,
  center = "median")
write_csv(df_cond_gdp_cpi, "conditional_flex_gdp_cpi.csv")

### Example B – Fed funds path + CPI-U YoY -----------------------------------
cond_fed_cpi <- list(
  "Fed Funds Rate" = list(mode = "rate_raw",
                          path = c(4.33, 3.83, 3.58, 3.33,
                                   rep(NA, 8))),
  "CPI-U"          = list(mode = "yoy_log",
                          path = c(NA, 4.0, 3.5, rep(NA, 9)))
)

forecasts$cond_fed_cpi <- conditional_flexible(
  data           = data_mat,
  df             = df,
  df_history_raw = df_raw,
  res            = res,
  p              = p,
  h              = h_cond,
  n_sim          = n_sim,
  var_names0     = var_names0,
  scenario_name  = "cond_fed_cpi",
  condition_specs= cond_fed_cpi,
  verbose        = TRUE)

df_cond_fed_cpi <- combine_fcst_with_history_raw(
  forecasts$cond_fed_cpi, df, h_cond,
  trans, var_names, var_names0,
  center = "median")
write_csv(df_cond_fed_cpi, "conditional_flex_fed_cpi.csv")

## ── 8.  Ex-post shock on baseline ------------------------------------------
shock_list_corepce <- list(
  "Core PCE Price Index" = list(quarter_start = 2,
                                quarter_end   = 3,
                                shift         = 1)   # +1 (log*100 units)
)

forecasts$shock_corepce <- shock_after_baseline_partial_pin(
  data           = data_mat,
  df             = df,
  df_history_raw = df_raw,
  res            = res,
  p              = p,
  h              = 12,
  n_sim          = n_sim,
  var_names0     = var_names0,
  scenario_name  = "shock_corepce_on_uncond",
  baseline_specs = NULL,         # start from unconditional
  shock_list     = shock_list_corepce,
  verbose        = TRUE)

df_shock_corepce_uncond_raw <- combine_fcst_with_history_raw(
  forecasts$shock_corepce, df, 12,
  trans, var_names, var_names0,
  center = "median")
write_csv(df_shock_corepce_uncond_raw, "shock_corepce_on_uncond.csv")

## ── 9.  Quick visual outputs -----------------------------------------------
scenario_data_list <- list(
  "Unconditional"       = df_combined_raw,
  "Cond GDP & CPI YoY"  = df_cond_gdp_cpi,
  "Cond FedFunds & CPI" = df_cond_fed_cpi,
  "Shock Core PCE"      = df_shock_corepce_uncond_raw
)

my_yoy_vars <- c("GDP", "PCE.Price.Index", "Core.PCE.Price.Index")
my_raw_vars <- c("Unemployment.Rate", "TR10y.Rate", "Fed.Funds.Rate")

df_long_all <- gather_scenarios_for_plot(
  scenario_data_list, df,
  yoy_vars = my_yoy_vars,
  raw_vars = my_raw_vars)

p_all <- plot_key_vars_multi_scenarios(
  df_long_all,
  main_title = "US BVAR – baseline, conditions & shocks")
ggsave("comparison_scenarios.png", p_all, width = 15, height = 6)

## Fan chart for real GDP ----------------------------------------------------
p_one_var <- make_fanchart_one_var_export(
  uncond_fcst, df,
  var_of_interest = "GDP",
  var_names       = var_names,
  var_names0      = var_names0,
  trans           = trans,
  start_plot      = 2015,
  center          = "median",
  export_csv_name = "gdp_unconditional_fanchart.csv")
ggsave("gdp_fcst_cond_fed_gdp.pdf", p_one_var, width = 8, height = 5)

## End of script #############################################################

