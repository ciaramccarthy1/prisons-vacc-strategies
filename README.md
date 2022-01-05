# prisons-vacc-strategies

Analysis of the impact of seven different vaccination scenarios on cases, QALY loss and deaths from COVID-19 in a category B local prison in England and Wales. 

This involved adapting CovidM, a transmission-dynamic model for SARS-CoV-2 transmission  developed by Davies et al., which has previously been used to evaluate the impact of vaccination in the UK at the population level. 

### Quick start guide
1. Install all required packages using `R/packages`. The version of CovidM used and instructions for installation can be found in the `covidm_for_fitting` folder. 

2. Set file paths within set-up script
 cm_path defines where covidm is stored e.g.
 ```{r eval=FALSE}
 cm_path = "C:/Users/CiaraMcCarthy/covidm-twodose/covidm_for_fitting/"
 ```
 pris_path defines where the contents of the prisons-july2021 Github directory are stored e.g.
 ```{r eval=FALSE}
 pris_path <- "~/prisons-tidy/Prisons-july2021"
 ```
3. Run `scripts/set-up-upgrading.R`. This loads all dependencies, including covidm. It sources the following other scripts in the repository:
      * `R/functions_may21.R` - required functions
      * `R/sensitivity.R` - assigns values for vaccine- and prison-related parameters; scales susceptibility to achieve desired R0.

4. `scripts/final-plots.R` includes code for all plots in the order that they appear in the manuscript. The script sources the following other scripts in the repository:
      * `scripts/vacc-scenarios.R` - runs model for each of the seven vaccination scenarios
      * `scripts/psa_define.R` - defines values for all parameters varied in probabilistic sensitivity analysis, based on values generated from Latin hypercube sampling
      * `scripts/vacc_short.R` - shortened version of `scripts/vacc-scenarios.R` script
      * `scripts/sensitivity.R` - as above

The sources for the QALY values used in `data/qalycalc-prisoners.xlsx` and `data/qalycalc-staff.xlsx` are described in McCarthy et al. [will include link to manuscript here]

### Citation
McCarthy et al., Impact of COVID-19 vaccination in prisons in England and Wales: a metapopulation model [will include link to manuscript here]

##### Citing CovidM:
Davies NG, Kucharski AJ, Eggo RM, Gimma A, Edmunds WJ, Jombart T, et al. Effects of non-pharmaceutical interventions on COVID-19 cases, deaths, and demand for hospital services in the UK: a modelling study. Lancet Public Health 2020. 

Sandmann FG, Davies NG, Vassall A, Edmunds WJ, Jit M, Sun FY, et al. The potential health and economic value of SARS-CoV-2 vaccination alongside physical distancing in the UK: a transmission model-based future scenario analysis and economic evaluation. Lancet Infect Dis. 2021.
