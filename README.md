# COVIVIS functions on R

[![GPL-3.0](https://custom-icon-badges.herokuapp.com/badge/license-GPL%203.0-8BB80A.svg?logo=law&logoColor=white)]()
[![R](https://custom-icon-badges.herokuapp.com/badge/R-198CE7.svg?logo=R&logoColor=white)]()
[![Static Badge](https://img.shields.io/badge/COVIVIS-v1.0.0-CCCCCC?link=https%3A%2F%2Fcovivis.soken.ac.jp%2F)](https://covivis.soken.ac.jp/)


## Overview
[COVIVIS](https://covivis.soken.ac.jp/) is an online tool that predicts the number of infected people from time-series data of the virus concentrations in waste water. It was developed and provided by Prof. Akira Sasaki at RCIES, SOKENDAI. This R script provides functions for estimating and predicting the number of disease cases, and COVIVIS uses some of the functions. The simulation results may not perfectly match the output from COVIVIS due to the use of random numbers and the specification of internal functions, but they are qualitatively identical. In addition, the script does not include the codes for data formating　or processing.

## Functions
### find_me_from_sd(m,sd)
- Calculation shape and scale parameters of Weibull distribution from the mean and the standard deviation of data.

### generate_onset(reporteddata,xmu,xsd)
- Estimate the number of disease onset from the number of reported cases with Monte Carlo method. The distribution of time delay from disease onset to reporting follows a Weibull distribution with mean $m$ and standard deviation $SD$.

### param_estim_by_lm(xdata,ydata)
- Parameters estimation with a linear regression model. 
The function accepts time series data as the explanatory variable in the first argument and time series data as the objective variable in the second argument. It performs a linear regression analysis of the input values converted to ordinary logarithms and returns the regression coefficients, the intercept $a$ and the slope $b$, and its residual standard deviation $r$ between the explanatory and the objective values.

### epi.prediction_by_lm(xdata,pa,pb,pr)
- Prediction using the estimated parameters with a linear regression model. This function applies a regression equation used the estimated parameters $a$ and $b$ to the input time series data by the function “param_estim_by_lm()” and returns the log predicted values. The residual standard deviation $r$ is used to evaluate the confidence and the prediction intervals of the predicted values.

### vi.param_estim_by_em_rev(sewagedata,onsetdata)
- Parameter estimation by Shedding profile model. The function accepts as its first argument time-series data on virus concentration and as its second argument time-series data on the number of disease onset. It estimates the parameters by optimization fitting (the maximum likelihood prediction) based on the shedding profile model.
 The results are returned as the parameters of the viral shedding curve for a person with the onset, $\nu$:the maximum rates of virus shedding (GC/L), $\omega$:the disease onset rate (1/day), $\gamma$:the recovery rate (1/day) and $r$:residual standard deviation of the number of disease onset.

### vi.prediction_by_em_rev(sewagedata,v,omega,gamma,pr)
- Prediction of the number of disease onset using estimated parameters by Shedding profile model.
 This function applies the estimated parameters $\nu$, $\gamma$ and $\omega$ to the input time series data of virus concentration and returns a log prediction based on the Shedding curve model. The fifth argument, the residual standard deviation $r$, is used to calculate confidence and prediction intervals for the predicted values. 

### iv.param_estim_by_em(onsetdata,sewagedata)
### iv.prediction_by_em(onsetdata,v,omega,gamma,pr)

## Requirement 
"COVIVIS-Functions.R" and "COVIVIS-Simulations.Rmd" require the following libraries to run:
- dplyr 
- magrittr
- nleqslv
- ggplot2
- DT
- nleqslv

## Note

## Author
Akiko Ohtsuki <br>
ohtsuki.t2.work@gmail.com

## License
The source code is licensed GPL v3 License, see LICENSE.<br>
COVIVIS is released under the GPL v3 License, see LICENSE.<br>
https://www.gnu.org/licenses/gpl-3.0
