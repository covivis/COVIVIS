# COVIVIS functions on R

[![GPL-3.0](https://custom-icon-badges.herokuapp.com/badge/license-GPL%203.0-8BB80A.svg?logo=law&logoColor=white)]()
[![R](https://custom-icon-badges.herokuapp.com/badge/R-198CE7.svg?logo=R&logoColor=white)]()
[![Static Badge](https://img.shields.io/badge/COVIVIS-v1.0-CCCCCC?link=https%3A%2F%2Fcovivis.soken.ac.jp%2F)](https://covivis.soken.ac.jp/)


## Overview
[COVIVIS](https://covivis.soken.ac.jp/) is an online tool that predicts the number of infectious from time-series data of the virus concentrations in waste water. It is developed and provided by Pf. Akira Sasaki at RCIES, SOKENDAI. This script provides functions for estimating and predicting the number of disease cases, and it is integrated into COVIVIS. The simulation results may not perfectly match the output from COVIVIS due to the use of random numbers and updates in parameters, but they are qualitatively identical. In addition, codes for data format are not included.

## Functions
### find_me_from_sd(m,sd)
 Calculation shape and scale parameters of Weibull distribution from the mean and the standard deviation of data.

### generate_onset(reporteddata,xmu,xsd)
 Estimate the number of disease onset from the number of reported cases with Monte Carlo method. The distribution of time delay from disease onset to reporting follows a Weibull distribution with mean $m$ and standard deviation $SD$.

### param_estim_by_lm(xdata,ydata)
 Parameters estimation with a linear regression model. 
The function accepts time series data as the explanatory variable in the first argument and time series data as the objective variable in the second argument. It performs a linear regression analysis of the input values converted to ordinary logarithms and returns the regression coefficients, the intercept $a$ and the slope $b$, and its residual standard deviation $r$ between the explanatory and the objective values.
### epi.prediction_by_lm(xdata,pa,pb,pr)
### vi.param_estim_by_em_rev(sewagedata,onsetdata)
### vi.prediction_by_em_rev(sewagedata,v,omega,gamma,pr)
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
Research Center for Integrative Evolutionary Sciense, SOKENDAI<br>
ohtsuki_akiko@soken.ac.jp

## License
The source code is licensed GPL v3 License, see LICENSE.<br>
COVIVIS is released under the GPL v3 License, see LICENSE.<br>
https://www.gnu.org/licenses/gpl-3.0
