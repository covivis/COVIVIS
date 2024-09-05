# COVIVIS functions on R

<p align="center">
<img  alt="COVIVIS-logo" src="https://covivis.soken.ac.jp/images/logo.svg" width="60%">
</p>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://custom-icon-badges.herokuapp.com/badge/R-198CE7.svg?logo=R&logoColor=white)]()
[![Static Badge](https://img.shields.io/badge/COVIVIS-v1.1.0-CCCCCC?link=https%3A%2F%2Fcovivis.soken.ac.jp%2F)](https://covivis.soken.ac.jp/)


## Overview
[COVIVIS](https://covivis.soken.ac.jp/) is a web tool that predicts the number of infected people from time-series data of the virus concentrations in waste water. This R script provides functions estimating some epidemiological parameters and predicting the number of disease cases with them. The simulation results may not perfectly match the output from COVIVIS due to the use of random numbers and the specification of internal functions, but they are qualitatively identical. In addition, the script does not include the codes for data formatting or processing.

## Functions
### find_me_from_sd(m,sd)
- Calculation shape and scale parameters of Weibull distribution from the mean and the standard deviation of data.

### generate_onset(reporteddata,xmu,xsd)
- Estimate the number of disease onsets from the number of reported cases with Monte Carlo method. The distribution of time delay from  onset of disease to reporting follows 　Weibull distribution with mean $m$ and standard deviation $SD$.

### generate_onset_sentinel(sentineldata,xmu,xsd)
- Estimate the number of disease onsets from sentinel reporting data.

### param_estim_by_lm(xdata,ydata)
- Parameters estimation with a linear regression model. 
The function accepts time series data as the explanatory variable $xdata$ in the first argument and time series data as the objective variable $ydata$ in the second argument. It performs a linear regression analysis of the input values converted to ordinary logarithms. The returned value of the function is a list of some statistics in the regression analysis -- $a$: the intercept, $b$: the slope, the $t$-score $tval$, the number of data $num.d$, the mean of $xdata$ $xmean$, the sum of the squares of the difference between each $xdata$ and the $xmean$ $sxx$, the unbiased variance of $ydata$ $uv$.

### epi.prediction_by_lm(xdata,pa,pb,pr)
- Prediction of the number of disease onsets using the estimated parameters with a linear regression model. This function applies a regression equation used the estimated parameters $a$ and $b$ to the input time series data by the function “param_estim_by_lm()” and returns the log predicted values. The residual standard deviation $r$ is used to evaluate the confidence and the prediction intervals of the predicted values.

### vi.param_estim_by_em(sewagedata,onsetdata)
- Parameter estimation by Shedding profile model. The function accepts as its first argument time-series data on virus concentration and as its second argument time-series data on the number of disease onsets. It estimates the parameters by optimization fitting (the maximum likelihood prediction) based on the shedding profile model.
 The results are returned as the parameters of the viral shedding curve for a person with the onset, $\nu$:the maximum rates of virus shedding (GC/L), $\omega$:the disease onset rate (1/day), $\gamma$:the recovery rate (1/day) and $r$:residual standard deviation of the number of disease onsets.

### vi.prediction_by_em_rev(sewagedata,v,omega,gamma,pr)
- Prediction of the number of disease onsets using estimated parameters by Shedding profile model.
 This function applies the estimated parameters $\nu$, $\gamma$ and $\omega$ to the input time series data of virus concentration and returns a log prediction based on the Shedding curve model. The fifth argument, the residual standard deviation $r$, is used to calculate confidence and prediction intervals for the predicted values. 

### iv.param_estim_by_em(onsetdata,sewagedata)
- Parameter estimation by Shedding profile model.

### iv.prediction_by_em(onsetdata,v,omega,gamma,pr)
- Prediction of the viral consentration using estimated parameters by Shedding profile model.

## Requirement 
"COVIVIS-Functions.R" requires the following libraries to run:
- dplyr 
- magrittr
- nleqslv
<!--
- ggplot2
- DT
- nleqslv
-->

## Note

## Author
Akiko Ohtsuki <br>
ohtsuki.t2.work@gmail.com

## License
The source code is licensed MIT License, see LICENSE for details.<br>
https://opensource.org/licenses/mit-license.php
