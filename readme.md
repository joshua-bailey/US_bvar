# US BVAR

This repository contains utilities for estimating a Bayesian VAR and producing scenario forecasts with U.S. macro data.

## Impulse response functions

Structural impulse responses can be generated using helper functions in `code/functions_irf_v6.R`.
The `generate_irf` function takes draws from an estimated BVAR, applies a Cholesky identification and
returns an array of responses. `plot_irf_core` quickly plots the responses for each variable.

### Example usage

```r
source("code/functions_data_v6.R")
source("code/functions_forecasting_v6.R")
source("code/functions_irf_v6.R")

# estimate model (simplified)
# data_obj <- fetch_data_from_fred(var_config)
# res <- estimate_bvar(data_obj$data, p = 4, n_draw = 2000, n_burn = 1000)

# generate IRF to a federal funds rate shock
# irf_ff <- generate_irf(res, p = 4, horizon = 12, shock_var = "Fed Funds Rate")
# plot_irf_core(irf_ff, var_names0 = var_config$var_names0, shock_name = "Fed Funds Rate")
```

The plot function uses median responses across draws and facets by variable.
