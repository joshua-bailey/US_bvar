# US‑BVAR Toolkit

A self‑contained workflow for constructing, estimating, and visualising a Bayesian VAR of the United States economy using publicly available FRED data.

---

## 1 Overview

| Step          | What happens                                                                                                    |
| ------------- | --------------------------------------------------------------------------------------------------------------- |
| **Data**      | Download monthly & quarterly series from the FRED API and convert them to seasonally‑adjusted quarterly levels. |
| **Transform** | Apply CBO‑style transformations (`log*100` for real levels / prices, `×100` for interest‑rate & % variables).   |
| **Estimate**  | Fit a Minnesota‑prior BVAR (with optional Slope‑Of‑Coefficients hyper‑prior) via the `BVAR` R package.          |
| **Forecast**  | Produce unconditional paths; add conditional scenarios or ex‑post shocks via a state‑space / Kalman smoother.   |
| **Output**    | Write CSVs and plots (scenario tiles, fan charts).                                       			  |


---

## 2 Dependencies

* R 4.2+
* CRAN packages: `tidyverse`, `zoo`, `fredr`, `BVAR`, `ggplot2`, `scales`
* FRED API key (see below)

Install any missing package with:

```r
install.packages(c("tidyverse","zoo","fredr","BVAR","ggplot2","scales"))
```

---

## 3 FRED access

Export your key once per session (or set in `~/.Renviron`):

```r
Sys.setenv(FRED_API_KEY = "3f66fc645f49a16dba8369e2539076f8")
```

The key in the repository is a public/demo token; feel free to replace it.

---

## 4 Running the model

1. Clone the repo 
2. Open `main.R` in R Studio.
3. Adjust `start_date`, `end_date`, or `var_config` if required.
4. `Source` the script – it will:

   * download data
   * estimate the BVAR
   * write forecasts to CSV
   * save plots in PNG/PDF

Run‑time with the default 20‑variable system (2 000 draws) is **≈ 2–3 min** on a modern laptop.

---

## 6 Conditional scenarios

`functions_conditions_v6.R` allows three pinning modes

| Mode        | Interpretation                | How to specify in `condition_specs`                                      |
| ----------- | ----------------------------- | ------------------------------------------------------------------------ |
| `rate_raw`  | Interest rate (level, % p.a.) | path supplied in **percent** (e.g. `3.50`)                               |
| `level_log` | Direct level (e.g. CPI index) | path supplied in **index units**                                         |
| `yoy_log`   | Year‑over‑year % change       | path supplied in **percent**; internally converted to levels then logged |

Missing values (`NA`) leave the draw unconstrained.

---

## 7 Shocks on top of a baseline

`shock_after_baseline_partial_pin()` shifts one or more series for selected quarters
and re‑filters the system while allowing other variables to respond endogenously.

```r
shock_list <- list(
  "Core PCE Price Index" = list(quarter_start = 2,
                                 quarter_end   = 3,
                                 shift         = 1)   # +1 (log*100)
)
```

---

## 8 Outputs

| File                                    | Description                                          |
| --------------------------------------- | ---------------------------------------------------- |
| `unconditional_forecast.csv`            | Median path + history (raw units).                   |
| `comparison_scenarios.png`              | Six‑panel tile plot comparing baseline vs scenarios. |
| `gdp_unconditional_fanchart.csv / .pdf` | 66/90 % fan chart for real GDP.                      |


---


