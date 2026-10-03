[![MIT License](https://custom-icon-badges.herokuapp.com/badge/license-MIT-8BB80A.svg?logo=law&logoColor=white)](LICENSE)
[![R](https://custom-icon-badges.herokuapp.com/badge/R-198CE7.svg?logo=R&logoColor=white)]()
[![COVIVIS](https://img.shields.io/badge/COVIVIS-v2.0-CCCCCC?link=https%3A%2F%2Fcovivis.soken.ac.jp%2F)](https://covivis.soken.ac.jp/)

# COVIVIS: Functions for Epidemiological and Wastewater Data Analysis

<p align="center">
  <img src="COVIVIS_logo.png" alt="COVIVIS logo" width="640">
</p>

`EpiPredFunctions20260819.R` is a collection of R functions for estimating disease onset, relating epidemiological data to wastewater virus concentrations, and evaluating prediction accuracy.

It includes onset back-projection, linear regression, and shedding profile models.

## Requirements

Install and load the following R packages.

```r
install.packages(c("dplyr", "magrittr", "nleqslv", "nloptr", "nnls"))

library(dplyr)
library(magrittr)
library(nleqslv)
library(nloptr)
library(nnls)
```

Load the functions in an R console or script.

```r
source("EpiPredFunctions20260819.R")
```

## Data Format

Most functions expect a `data.frame` with dates in the first column and numeric values in the second. Some functions rename columns internally, so place the date column first and the value column second.

| Use | Expected columns |
| --- | --- |
| Reported, onset, or infected cases | `date`, count |
| Sentinel surveillance data | `date`, value per sentinel |
| Wastewater data | `date`, virus concentration |
| Regression data | `date`, positive value |

Inputs used for log transformation must be positive. Missing values are removed or handled depending on the function.

## Functions

### 1. Back-projecting Disease Onset from Report Dates

#### `find_me_from_sd(xmean, xsd)`

Numerically finds the shape and scale of a Weibull distribution from the mean (`xmean`) and standard deviation (`xsd`) of reporting delay.

- **Returns**: `c(shape, scale)`
- **Used by**: `generate_onset()` and `generate_onset_sentinel()`

<br>

#### `generate_onset(reporteddata, xmu, xsd)`

Uses Monte Carlo back-projection with a Weibull reporting-delay distribution to estimate daily onset counts from daily reported cases.

- `reporteddata`: two columns containing date and reported cases
- `xmu`, `xsd`: mean and standard deviation of reporting delay
- **Returns**: a data frame with `date` and estimated `onset`

<br>

#### `generate_onset_sentinel(sentineldata, xmu, xsd)`

Temporarily scales sentinel values by 100 or 1,000, applies the same method as `generate_onset()`, and restores the original scale.

- `sentineldata`: two columns containing date and sentinel value
- **Returns**: a data frame with `date` and estimated `onset`

<br>

### 2. Time-series Processing and Interval Estimation

#### `lerp_data(timeseries)`

Converts weekly sentinel data to a daily series. It divides the second column by seven and linearly interpolates values between dates.

- `timeseries`: two columns containing date and weekly value
- **Returns**: a data frame with daily `date` and `value`

<br>

#### `calc_ts_CI(obs, pred, n_param = 3, conf_level = 0.95)`

Calculates residual variance, effective sample size or effective sample size based on lag-1 autocorrelation, and confidence-interval (CI) and prediction-interval (PI) widths. Inputs are assumed to be on a log scale.

- `obs`, `pred`: observed and predicted values of equal length
- `n_param`: number of estimated parameters (default: 3)
- `conf_level`: confidence level (default: 0.95)
- **Returns**: a list including `ve`, `n_eff`, `rho1`, `ci_width`, and `pi_width`

<br>

#### `plot_result_cipi(dataframe, ci, pi)`

Adds CI and PI limits to predictions on a log scale.

- `dataframe`: two columns containing date and predicted value
- **Returns**: `date`, `predicted.log`, `CI_lower`, `CI_upper`, `PI_lower`, and `PI_upper`

<br>

#### `plot_result_cipi_nonlog(dataframe, ci, pi)`

Converts base-10 log predictions and interval widths to the original scale for plotting.

- `dataframe`: two columns containing date and base-10 log prediction
- **Returns**: `date`, `predicted`, `CI_lower`, `CI_upper`, `PI_lower`, and `PI_upper`

<br>

#### `weekly_prediction(weeklyobserved, dailypredicted)`

Sums predicted daily cases over the seven days ending on each weekly observation date.

- `weeklyobserved`: weekly observations with a `date` column
- `dailypredicted`: daily predictions with `date` and `estimated_cases` columns
- **Returns**: a data frame with `date` and `pred_sum`

<br>

### 3. Regression and R-squared

#### `Rsquared_in_lm(xdata, ydata, pa, pb)`

Uses supplied regression coefficients `pa` (intercept) and `pb` (slope) to calculate SSE, SSR, SST, and two R-squared values on the base-10 log scale after matching data by date.

- **Returns**: `c(SSE, SSR, SST, R2_1, R2_2)`
- `R2_1 = 1 - SSE / SST`
- `R2_2 = SSR / SST`

<br>

#### `Rsquared(observeddata, predicteddata)`

Calculates SSE, SSR, SST, and two R-squared values from observed and predicted values. It returns an error when vector lengths differ or SST is zero.

- **Returns**: `c(SSE, SSR, SST, R2_1, R2_2)`

<br>

### 4. Linear Regression Model

#### `param_estim_by_LM(xdata, ydata)`

Fits a simple linear model, `log10(y) ~ log10(x)`, to two series matched by date. It also returns statistics needed to calculate prediction intervals.

- **Returns**: `c(intercept, slope, t_value, n, x_mean, Sxx, residual_variance, SSE, SSR, SST, R2_1, R2_2)`

<br>

#### `epi.prediction_by_LM(xdata, pa, pb, tval, num.d, xmean, sxx, uv, rsq)`

Uses the regression parameters and statistics returned by `param_estim_by_LM()` to calculate base-10 log predictions, 95% CIs, and 95% PIs for `xdata`.

- **Returns**: a data frame with `date`, `xdata`, `prediction.y`, `ci.up`, `ci.lw`, `pi.up`, and `pi.lw`
- Note: `rsq` is not used by the current implementation.

<br>

### 5. Shedding Profile Model

#### `param_estim_by_SPM(epidemicdata, sewagedata)`

Fits a forward model from epidemiological data (I) to wastewater virus concentration (V) using a Gaussian shedding profile. Weekly epidemiological data are interpolated to daily data when needed.

- `epidemicdata`: two columns containing date and infected cases; weekly or daily data are accepted
- `sewagedata`: two columns containing date and virus concentration
- **Returns**: a list with two elements
  - `params`: `v`, `m`, `sigma`, `ci_width`, `pi_width`, and `R2`
  - `predicted`: daily predictions for the fitting period with `date` and `pred.virus`

<br>

#### `virus_prediction_by_SPM(epidemicdata, v, m, sigma)`

Predicts daily wastewater virus concentrations from infected cases with specified shedding-profile parameters. Predictions are returned on the base-10 log scale.

- `v`: magnitude of the shedding profile
- `m`: time shift of the shedding peak
- `sigma`: width of the shedding profile
- **Returns**: a data frame with `date` and `pred.virus`

<br>

#### `infected_prediction_by_SPM(sewagedata, est_v, est_m, est_sigma, lambda = 10^2)`

Estimates infected cases (I) from wastewater virus concentrations (V). It minimizes log-scale reconstruction error and a second-difference penalty on the log case series. Estimated case counts are always positive and are multiplied by a symptomatic fraction of `2/3` before return.

- `sewagedata`: wastewater data with `date` and `virus` columns
- `est_v`, `est_m`, `est_sigma`: shedding-profile parameters estimated by `param_estim_by_SPM()`
- `lambda`: smoothing penalty for the log case series (default: `10^2`). Larger values produce smoother estimates.
- **Returns**: a list with two elements
  - `resI`: a data frame with `date` and `estimated_cases`
  - `resV`: a numeric vector of wastewater virus concentrations reconstructed from the estimated cases

## Basic Example

The following example predicts wastewater virus concentration from epidemiological data, then estimates infected cases from wastewater data.

```r
source("EpiPredFunctions20260819.R")

# epidemicdata: date, infected
# sewagedata:   date, virus
fit <- param_estim_by_SPM(epidemicdata, sewagedata)

# Forward prediction of virus concentration
virus_pred <- virus_prediction_by_SPM(
  epidemicdata,
  v = fit$params$v,
  m = fit$params$m,
  sigma = fit$params$sigma
)

# Backward estimation of infected cases
infected_fit <- infected_prediction_by_SPM(
  sewagedata,
  est_v = fit$params$v,
  est_m = fit$params$m,
  est_sigma = fit$params$sigma
)

# Estimated infected cases
infected_pred <- infected_fit$resI
```

## Notes

- Results from `generate_onset()` and `generate_onset_sentinel()` vary because they use Monte Carlo simulation. Use `set.seed()` before running them when reproducible results are needed.
- Some functions rename columns internally. Check column order and date types before use.
- Interpret CIs, PIs, and R-squared values on the scale used by each function.
- This code is intended for research and analysis. Check data quality, model assumptions, and external validity before using results for epidemiological interpretation or operations.

## License

This source code is released under the [MIT License](LICENSE). See `LICENSE` for details.
