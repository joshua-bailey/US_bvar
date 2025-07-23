########################
## main.R
########################

# this has soc forecasting 
# need to do the IRFs inc the conditioning problem

# nb. the shock forecasts seem to work better with some conditioned paths. Maybe more to look into there

# ****************************************************************
# US macro forecasting and scenarios using a BVAR model
#
# this script shows how to use the new bvar pipeline to:
#   1) load and transform data from fred
#   2) estimate a bvar
#   3) generate unconditional forecasts
#   4) generate conditional/scenario forecasts (e.g. yoy constraints or fixed interest rate)
#   5) layer shocks on top of an existing baseline
#   6) produce example multi-scenario plots and single-variable fan charts
#
#
# Author: Original BVAR work by Hark Yoo, CBO, updated by Joshua Bailey for FRED data + new conditional forecasting function with shock scenarios
# ****************************************************************
library(tidyverse)
library(readxl)
library(openxlsx)
library(zoo)
library(data.table)
library(BVAR)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(fredr)
library(dplyr)
library(purrr)
library(scales)

# set seed 
set.seed(1234)

# source the function files
source("code/functions_data_v6.R")
source("code/functions_forecasting_v6.R")
source("code/functions_conditions_v6.R")
source("code/functions_output_v6.R")
source("code/functions_irf_v6.R")

#-----------------------------------
# 1. load/transform data from fred
#-----------------------------------
# define a var_config data frame that describes each variable:
#   - var_names0 is the descriptive name
#   - fred_id is the series id in fred
#   - freq is 'monthly' or 'quarterly'
#   - transform is 'log100' or 'raw'
#
# NIPA Data (Macro variables)
nipa <- data.frame(
  var_names0 = c(
    "GDP", 
    "Consumption", 
    "Investment (nonres)", 
    "Exports", 
    "Imports", 
    "GDP Deflator", 
    "PCE Price Index", 
    "Core PCE Price Index", 
    "CPI-U", 
    "Potential GDP",
    "Real Govt. Consumption Expenditures & Gross Investment"
  ),
  fred_id = c(
    "GDPC1", 
    "PCECC96", 
    "GPDIC1", 
    "EXPGSC1", 
    "IMPGSC1", 
    "GDPDEF", 
    "PCEPI", 
    "PCEPILFE", 
    "CPIAUCSL", 
    "GDPPOT",
    "GCEC1"
  ),
  freq = rep("quarterly", 11),
  transform = rep("log100", 11),
  stringsAsFactors = FALSE
)

# Labour Market Variables
labour <- data.frame(
  var_names0 = c(
    "Payroll Employment",
    "Labor Force",
    "Wages",
    "Unemployment Rate"
  ),
  fred_id = c(
    "PAYEMS",
    "CLF16OV",
    "A4102C1Q027SBEA",
    "UNRATE"
  ),
  freq = c("monthly", "monthly", "quarterly", "monthly"),
  transform = c("log100", "log100", "log100", "raw"),
  stringsAsFactors = FALSE
)

# Financial Market Variables
financial <- data.frame(
  var_names0 = c(
    "Fed Funds Rate",
    "TR3m Rate",
    "TR10y Rate",
    "TR1y Rate",
    "Moody's Aaa",
    "Industrial Production"
  ),
  fred_id = c(
    "FEDFUNDS",
    "TB3MS",
    "GS10",
    "GS1",
    "AAA",
    "INDPRO"
  ),
  freq = c("monthly", "monthly", "monthly", "monthly", "monthly", "monthly"),
  transform = c("raw", "raw", "raw", "raw", "raw", "log100"),
  stringsAsFactors = FALSE
)

# Additional Variables
additional <- data.frame(
  var_names0 = c(
    "Cleveland Fed 10-yr Inflation Expectations",
    "Economic Policy Uncertainty Index",
    "WTI Spot Oil Price"
  ),
  fred_id = c(
    "EXPINF10YR",
    "USEPUINDXM",
    "MCOILWTICO"
  ),
  freq = rep("monthly", 3),
  transform = c("raw", "log100", "log100"),
  stringsAsFactors = FALSE
)

# Combine all groups into one configuration
var_config <- rbind(nipa, labour, financial, additional)



# specify start/end in decimal year format
start_date = 1986
end_date   = 2025.00  # e.g. 2024q4

# fetch data
dfs <- fetch_data_from_fred(var_config, start_date, end_date)
df        = dfs$df       # transformed data frame
data      = dfs$data     # numeric matrix (for bvar)
trans     = dfs$trans
var_names = dfs$var_names
var_names0= dfs$var_names0

# create df_raw by reversing transforms for yoy constraints or just for clarity
df_raw <- df
for(j in seq_along(var_names)){
  col_j <- var_names[j]
  if(trans[j] == "log100"){
    df_raw[[col_j]] <- exp(df_raw[[col_j]] / 100)
  } else if(trans[j] == "raw"){
    df_raw[[col_j]] <- df_raw[[col_j]] / 100
  }
  # special scaling for potential gdp
  if(var_names0[j] == "Potential GDP"){
    df_raw[[col_j]] <- df_raw[[col_j]] / 10
  }
}

#-----------------------------------
# 2. bvar estimation
#-----------------------------------
p       = 4
n_draw  = 2000L
n_burn  = 1000L
n_sim   = n_draw - n_burn  # we do not hard-code 1000 again

# we estimate up to row 184 => example for historical window
data_est = data[1:nrow(data),]

res = estimate_bvar(
  data     = data_est,
  p        = 4,
  n_draw   = 2000L,
  n_burn   = 1000L,
  verbose  = TRUE,
  use_soc  = TRUE
)

#-----------------------------------
# 3. unconditional forecast
#-----------------------------------
# generate an unconditional forecast with horizon h
h = 41
forecasts = list()

forecasts[["base"]] = generate_unconditional(data, df, res, p, h, n_sim, var_names0)

# save that forecast to csv
uncond_fcst = forecasts[["base"]]
df_combined_raw = combine_fcst_with_history_raw(
  fcst_array     = uncond_fcst,
  df_transformed = df,
  h              = h,
  trans          = trans,
  var_names      = var_names,
  var_names0     = var_names0,
  center         = "median"
)
write.csv(df_combined_raw, "unconditional_forecast.csv")
# 
# #-----------------------------------
# # 4. flexible conditional forecast (example usage)
# #-----------------------------------

h_cond <- 12  # horizon in quarters for the scenario
# 
# # (c) condition: fed funds + yoy gdp
# cond_fed_gdp <- list(
#   "Fed Funds Rate" = list(
#     mode = "rate_raw",
#     path = c(4.33, 4.08, 3.83, 3.58, 3.33, NA, NA, NA, NA, NA, NA, NA)
#   ),
#   "GDP" = list(
#     mode = "yoy_log",
#     path = c(2.22, 1.66, 1.1, 0.76, 1.01, 1.25, 1.6, 2, 2.06, 2.09, 2.01, 1.85)
#   )
# )
# fcst_cond_fed_gdp <- conditional_flexible(
#   data           = data,
#   df             = df,
#   df_history_raw = df_raw,
#   res            = res,
#   p              = 4,
#   h              = h_cond,
#   n_sim          = n_sim,
#   var_names0     = var_names0,
#   scenario_name  = "cond_fed_gdp",
#   condition_specs= cond_fed_gdp,
#   verbose        = TRUE
# )
# df_cond_fed_gdp <- combine_fcst_with_history_raw(
#   fcst_array     = fcst_cond_fed_gdp,
#   df_transformed = df,
#   h              = h_cond,
#   trans          = trans,
#   var_names      = var_names,
#   var_names0     = var_names0,
#   center         = "median"
# )
# write.csv(df_cond_fed_gdp, "conditional_flex_fed_gdp.csv")
# forecasts[["cond_fed_gdp"]] <- fcst_cond_fed_gdp

# (d) march forecast 
cond_gdp_cpi <- list(
  "GDP" = list(
    mode = "yoy_log",
    path = c(-0.16, -1.2, -1.56, -0.92, 0.99, 1.77, 2.05, NA, NA, NA, NA, NA)
  ),
  "CPI-U" = list(
    mode = "yoy_log",
    path = c(NA, 4, 3.5,	NA, NA,	NA,	NA,	NA, NA, NA, NA, NA)
  )
)



fcst_cond_gdp_cpi <- conditional_flexible(
  data           = data,
  df             = df,
  df_history_raw = df_raw,
  res            = res,
  p              = 4,
  h              = h_cond,
  n_sim          = n_sim,
  var_names0     = var_names0,
  scenario_name  = "cond_gdp_cpi",
  condition_specs= cond_gdp_cpi,
  verbose        = TRUE
)
df_cond_gdp_cpi <- combine_fcst_with_history_raw(
  fcst_array     = fcst_cond_gdp_cpi,
  df_transformed = df,
  h              = h_cond,
  trans          = trans,
  var_names      = var_names,
  var_names0     = var_names0,
  center         = "median"
)
write.csv(df_cond_gdp_cpi, "conditional_flex_gdp_cpi.csv")
forecasts[["cond_gdp_cpi"]] <- fcst_cond_gdp_cpi


# (d) march forecast 
cond_fed_cpi <- list(
  "Fed Funds Rate" = list(
    mode = "rate_raw",
    path = c(4.33, 3.83, 3.58, 3.33, NA, NA, NA, NA, NA, NA, NA, NA)
  ),
  "CPI-U" = list(
    mode = "yoy_log",
    path = c(NA,	4, 3.5,	NA, NA,	NA,	NA,	NA, NA, NA, NA, NA)
  )
)


fcst_cond_fed_cpi <- conditional_flexible(
  data           = data,
  df             = df,
  df_history_raw = df_raw,
  res            = res,
  p              = 4,
  h              = h_cond,
  n_sim          = n_sim,
  var_names0     = var_names0,
  scenario_name  = "cond_fed_cpi",
  condition_specs= cond_fed_cpi,
  verbose        = TRUE
)
df_cond_fed_cpi <- combine_fcst_with_history_raw(
  fcst_array     = fcst_cond_fed_cpi,
  df_transformed = df,
  h              = h_cond,
  trans          = trans,
  var_names      = var_names,
  var_names0     = var_names0,
  center         = "median"
)
write.csv(df_cond_fed_cpi, "conditional_flex_fed_cpi.csv")
forecasts[["cond_fed_cpi"]] <- fcst_cond_fed_cpi

#-----------------------------------
# 5. adding shocks
#-----------------------------------
# we can shock certain variables after building a baseline
# example: (a) shock core pce on unconditional, (b) shock 10y rate on fed funds scenario

# (i) shock core pce 
shock_list_corepce <- list(
  "Core PCE Price Index" = list(
    quarter_start = 2,
    quarter_end   = 3,
    shift         = 1  
  )
)
fcst_shock_corepce_uncond <- shock_after_baseline_partial_pin(
  data           = data,
  df             = df,
  df_history_raw = df_raw,
  res            = res,
  p              = 4,
  h              = 12,
  n_sim          = n_sim,
  var_names0     = var_names0,
  scenario_name  = "shock_corepce_on_uncond",
  baseline_specs = NULL,  
  shock_list     = shock_list_corepce,
  verbose        = TRUE
)
df_shock_corepce_uncond_raw <- combine_fcst_with_history_raw(
  fcst_array     = fcst_shock_corepce_uncond,
  df_transformed = df,
  h              = 12,
  trans          = trans,
  var_names      = var_names,
  var_names0     = var_names0,
  center         = "median"
)
write.csv(df_shock_corepce_uncond_raw, "shock_corepce_on_uncond.csv")
forecasts[["shock_corepce_on_uncond"]] <- fcst_shock_corepce_uncond


# (ii) shock core pce, EPU increase
shock_list_corepce_infexp <- list(
  "PCE Price Index" = list(
    quarter_start = 2,
    quarter_end   = 3,
    shift         = 2  # +0.5 in log(*100) => about +0.5% level
  ),
  "Economic Policy Uncertainty Index" = list(
    quarter_start = 1,
    quarter_end   = 3,
    shift         = 50  # +0.5 in log(*100) => about +0.5% level
  )
)
fcst_shock_corepce_infexp_uncond <- shock_after_baseline_partial_pin(
  data           = data,
  df             = df,
  df_history_raw = df_raw,
  res            = res,
  p              = 4,
  h              = 12,
  n_sim          = n_sim,
  var_names0     = var_names0,
  scenario_name  = "shock_corepce_infexp_on_uncond",
  baseline_specs = NULL,  
  shock_list     = shock_list_corepce_infexp,
  verbose        = TRUE
)
df_shock_corepce_infexp_uncond_raw <- combine_fcst_with_history_raw(
  fcst_array     = fcst_shock_corepce_infexp_uncond,
  df_transformed = df,
  h              = 12,
  trans          = trans,
  var_names      = var_names,
  var_names0     = var_names0,
  center         = "median"
)
write.csv(df_shock_corepce_infexp_uncond_raw, "shock_corepce_infexp_on_uncond.csv")
forecasts[["shock_corepce_infexp_on_uncond"]] <- fcst_shock_corepce_infexp_uncond






#-----------------------------------
# 6. multi-scenario tile plot examples
#-----------------------------------
# here we compare yoy-based transformations for gdp or inflation,
# and raw rates for e.g. unemployment, fed funds, 10y yield.


# scenario_data_list <- list(
#   "Unconditional" = df_combined_raw,
#   "FedGDP" = df_cond_fed_gdp,
#   "FedGDPCpi" = df_cond_fed_gdp_cpi
# )

scenario_data_list <- list(
  "Unconditional" = df_combined_raw
)

my_yoy_vars <- c("GDP", "PCE.Price.Index", "Core.PCE.Price.Index")
my_raw_vars <- c("Unemployment.Rate", "TR10y.Rate", "Fed.Funds.Rate")

# gather scenarios into long format
df_long_all <- gather_scenarios_for_plot(
  scenario_list = scenario_data_list,
  df_original   = df,
  yoy_vars      = my_yoy_vars,
  raw_vars      = my_raw_vars
)

# plot the multi-scenario tile
p_all <- plot_key_vars_multi_scenarios(
  df_long_all,
  main_title = "US unconditional Q2 forecast"
)
print(p_all)
ggsave("comparison_scenarios.png", p_all, width = 15, height = 6)

#-----------------------------------
# 7. single-variable fan chart example
#-----------------------------------
# we create a fan chart for unemployment rate from a particular scenario

p_one_var <- make_fanchart_one_var_export(
  fcst_array      = uncond_fcst,
  df              = df,
  var_of_interest = "GDP",
  var_names       = var_names,
  var_names0      = var_names0,
  trans           = trans,
  start_plot      = 2015,          # earliest year to show
  center          = "median",
  impose_zero     = FALSE,
  export_csv_name = "gdp_fcst_cond_fed_gdp.csv"
)
print(p_one_var)
ggsave("gdp_fcst_cond_fed_gdp.pdf", p_one_var, width = 8, height = 5)


#-----------------------------------
# 8. structural impulse response example
#-----------------------------------
# compute an IRF to a Fed Funds Rate shock and plot responses

irf_ffr <- generate_irf(
  res       = res,
  p         = p,
  horizon   = 12,
  shock_var = "Fed Funds Rate"
)

p_irf_ffr <- plot_irf_core(
  irf_array  = irf_ffr,
  var_names0 = var_names0,
  shock_name = "Fed Funds Rate"
)
print(p_irf_ffr)
ggsave("irf_fedfunds.png", p_irf_ffr, width = 10, height = 6)



