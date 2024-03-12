# COVIVIS functions on R

[![GPL-3.0](https://custom-icon-badges.herokuapp.com/badge/license-GPL%203.0-8BB80A.svg?logo=law&logoColor=white)]()

## Overview
COVIVIS is an online tool that predicts the number of infectious from time-series data of the virus concentrations in waste water. It is developed and provided by Pf. Akira Sasaki at RCIES, SOKENDAI. This script provides functions for estimating and predicting the number of disease cases, and it is integrated into COVIVIS. The simulation results may not perfectly match the output from COVIVIS due to the use of random numbers and updates in parameters, but they are qualitatively identical. In addition, codes for data format are not included.

## Functions
### find_me_from_sd(m,sd)
### generate_onset(reporteddata,xmu,xsd)
### param_estim_by_lm(xdata,ydata)
### epi.prediction_by_lm(xdata,pa,pb,pr)
### vi.param_estim_by_em_rev(sewagedata,onsetdata)
### vi.prediction_by_em_rev(sewagedata,v,omega,gamma,pr)
### iv.param_estim_by_em(onsetdata,sewagedata)
### iv.prediction_by_em(onsetdata,v,omega,gamma,pr)

## Requirement 
"COVIVISFunction.ver.1.0.R" and "COVIVISSimulations---.Rmd" require the following libraries to run:
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
This R source codes are released under the GPL v3 License, see LICENSE.<br>
https://www.gnu.org/licenses/gpl-3.0
