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
| **Output**    | Write tidy CSVs and publication‑grade plots (scenario tiles, fan charts).                                       |

The design philosophy is **transparent** (no hidden steps), **reproducible** (single `main.R` script) and **flexible** (easy to add/replace series or scenarios).

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

## 4 Variable mapping

All price & quantity indices are **rebased to 2017 = 100** (FRED’s current convention).

| Human label                   | Safe column            | FRED ID  | Freq. | Transform |
| ----------------------------- | ---------------------- | -------- | ----- | --------- |
| Real GDP                      | `GDP`                  | GDPC1    | Q     | log100    |
| Real PCE                      | `Consumption`          | PCECC96  | Q     | log100    |
| Real Private Fixed Investment | `Investment.nonres.`   | GPDIC1   | Q     | log100    |
| Real Exports                  | `Exports`              | EXPGSC1  | Q     | log100    |
| Real Imports                  | `Imports`              | IMPGSC1  | Q     | log100    |
| GDP Deflator                  | `GDP.Deflator`         | GDPDEF   | Q     | log100    |
| PCE Price Index               | `PCE.Price.Index`      | PCEPI    | Q     | log100    |
| Core PCE Price Index          | `Core.PCE.Price.Index` | PCEPILFE | Q     | log100    |
| CPI‑U                         | `CPI.U`                | CPIAUCSL | Q     | log100    |
| Potential GDP                 | `Potential.GDP`        | GDPPOT   | Q     | log100    |
| Federal Funds Rate            | `Fed.Funds.Rate`       | FEDFUNDS | M → Q | ×100      |
| …                             | …                      | …        | …     | …         |

> **Transform codes**  `log100` = `log(x) * 100`;   `×100` = levels or rates multiplied by 100 so that “1 percent” is stored as `100`.

The full mapping is in `main.R`; edit or extend it there.

---

## 5 Running the model

1. Clone the repo and `cd` into it.
2. Open `main.R` in R Studio (or your favourite IDE).
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

All CSVs follow a tidy, wide format with a `yq` column (`zoo::yearqtr`).

---

## 9 Licence

The code is released under the MIT Licence; FRED data are public domain.
