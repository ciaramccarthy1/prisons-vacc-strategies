# Value of COVID-19 vaccine - prison modelling adaptation
# original code of "covidm" by ND, code adapted 31.03.2020

# clear global environment
# rm(list = ls())

## BASIC SET UP - DOESN'T CHANGE ##
library(ggplot2)
library(tidyverse)
library(data.table)
library(Rcpp)
library(rlang)
library(stringr)
library(cowplot)
library(magrittr)
library(RcppGSL)
library(qs)
library(viridis)
library(colorblindr)
library(readxl)
library(epiR)
library(crayon)
library(lubridate)
library(HDInterval)
library(ggthemes)
library(ggpubr)
library(gridExtra)
library(tictoc)
library(lhs)
library(dplyr)
library(Rcpp)

# load epi model (covidm)
cm_path = "C:/Users/CiaraMcCarthy/Documents/prisons-vacc-strategies/covidm_for_fitting/"
cm_force_rebuild = T;
cm_build_verbose = T;
cm_version = 2;
source(paste0(cm_path, "/R/covidm.R"))

pris_path <- "~/Documents/prisons-vacc-strategies"
save_path <- "C:/Users/CiaraMcCarthy/Documents/prisons-vacc-strategies/outputs/"

setwd("~/Documents/prisons-vacc-strategies")

## FUNCTIONS ##
# Contains lockdown function, labelling for plots, beta and gamma moment functions for PSA

source(paste0(pris_path, "/R/functions_may21.R"))

## QALY values 
qalycalc.prisoners <- read_xlsx(paste0(pris_path,"/data/qalycalc-prisoners.xlsx"))
aefi.minor <- (c(rep(0,3), rep(0.87, 8), rep(0.75, 3), rep(0.63, 2)))*(1/365.25)*0.5 # frequency*qaly loss
aefi.fatal <- c(rep(0,3), rep(1/50000, 5), rep(1/100000, 8))*as.vector(qalycalc.prisoners$qaly.value[qalycalc.prisoners$compartment=="death_o"]) # frequency*qaly.loss
qalycalc.prisoners <- qalycalc.prisoners %>% select(-alpha, -beta)
qalycalc.prisoners <- rbind(qalycalc.prisoners, data.frame(
                    group=unique(qalycalc.prisoners$group),
                    compartment="aefi.minor",
                    qaly.value=aefi.minor))
qalycalc.prisoners <- rbind(qalycalc.prisoners, data.frame(
  group=unique(qalycalc.prisoners$group),
  compartment="aefi.fatal",
  qaly.value=aefi.fatal))

qalycalc.staff <- read_xlsx(paste0(pris_path,"/data/qalycalc-staff.xlsx"))
aefi.minor <- (c(rep(0,3), rep(0.87, 8), rep(0.75, 3), rep(0.63, 2)))*(1/365.25)*0.5 # frequency*qaly loss
aefi.fatal <- c(rep(0,3), rep(1/50000, 5), rep(1/100000, 8))*as.vector(qalycalc.staff$qaly.value[qalycalc.staff$compartment=="death_o"]) # frequency*qaly.loss
qalycalc.staff <- qalycalc.staff %>% select(-alpha, -beta)
qalycalc.staff <- rbind(qalycalc.staff, data.frame(
  group=unique(qalycalc.staff$group),
  compartment="aefi.minor",
  qaly.value=aefi.minor))
qalycalc.staff <- rbind(qalycalc.staff, data.frame(
  group=unique(qalycalc.staff$group),
  compartment="aefi.fatal",
  qaly.value=aefi.fatal))

qalycalc <- qalycalc.prisoners %>% mutate(population="(C) prisoners")
qalycalc <- rbind(qalycalc, qalycalc.staff %>% mutate(population="(A) non-prisoner-facing staff"))
qalycalc <- rbind(qalycalc, qalycalc.staff %>% mutate(population="(B) prisoner-facing staff"))

# By scenario


# loading lockdown function
# fun_lockdown(lockdown_ini  = 5)

# Copied the data folder from covid-uk:
cm_matrices     <- readRDS(paste0(cm_path, "/data/all_matrices.rds"));
cm_populations  <- readRDS(paste0(cm_path, "/data/wpp2019_pop2020.rds"));
cm_structure_UK <- readRDS(paste0(cm_path, "/data/structure_UK.rds"));

nr_pops = 3

params = cm_parameters_SEI3R(rep(cm_uk_locations("UK", 0), nr_pops),
                                 date_start = "2020-02-01",
                                 deterministic = TRUE)
timehorizon <- 1
delay <- 365.25
params$date1 <- as.Date(365.25+round(timehorizon*365.25), origin=params$date0)
params$time1 <- as.Date(365.25+round(timehorizon*365.25), origin=params$date0)

## Mixing within populations:
for(i in 1:3){params$pop[[i]]$contact <- c(1,1,1,1)}
for(i in 1:nr_pops) {
params$pop[[i]]$matrices$work[4:16,4:16] <- 1
params$pop[[i]]$matrices$work[1:3,1:3] <- 0
params$pop[[i]]$matrices$school[1:16, 1:16] <- 0
params$pop[[i]]$matrices$other[1:16, 1:16] <- 0
params$pop[[i]]$matrices$home[1:16,1:16] <- 0
}



## Mixing between populations:
params$travel[1,] <- c(0.8, 0.2, 0)
params$travel[2,] <- c(0.2, 0.4, 0.4)
params$travel[3,] <- c(0, 0.4, 0.6)

# relabel populations
params$pop[[1]]$name = "(A) non-prisoner-facing staff"
params$pop[[2]]$name = "(B) prisoner-facing staff"
params$pop[[3]]$name = "(C) prisoners"

# Staff number and age structure:

Total_staff <- sum(10347, 12279, 11502, 14195, 4115)
age_staff <- c(0,0,0,
               10347/Total_staff/3, 10347/Total_staff/3, 10347/Total_staff/3,
               12279/Total_staff/2, 12279/Total_staff/2,
               11502/Total_staff/2, 11502/Total_staff/2,
               14195/Total_staff/2, 14195/Total_staff/2,
               4115/Total_staff/4, 4115/Total_staff/4, 4115/Total_staff/4, 4115/Total_staff/4)

params$pop[[1]]$size <- round(70 * age_staff) # local NPF
params$pop[[2]]$size <- round(315 * age_staff) # local PF

# Prisoner number and age structure:

prisoner_pop <- 824 # local prison
age_prisoners <- c(0,0,0,
                   53/824,# 15-20
                   102/824, # 21-24
                   154/824, # 25-29
                   274/2/824, # 30-39 divided by 2
                   274/2/824, # 30-39 divided by 2
                   147/2/824, # 40-49 divided by 2
                   147/2/824, # 40-49 divided by 2
                   65/2/824, # 50-59 divided by 2
                   65/2/824, # 50-59 divided by 2
                   19/2/824, # 60-69 divided by 2
                   19/2/824, # 60-69 divided by 2
                   10/2/824, # 70+ divided by 2
                   10/2/824) # 70+ divided by 2

# assume split up by average age distribution
params$pop[[3]]$size <- round(prisoner_pop * age_prisoners)

##### DEATHS #####

# Health burden processes
probs = fread(
  "Age,Prop_symptomatic,IFR,Prop_inf_hosp,Prop_inf_critical,Prop_critical_fatal,Prop_noncritical_fatal,Prop_symp_hospitalised,Prop_hospitalised_critical
10,0.66,8.59E-05,0.002361009,6.44E-05,0.32,0,0,0.17
20,0.66,0.000122561,0.003370421,9.19E-05,0.32,9.47E-04,0.007615301,0.17
30,0.66,0.000382331,0.010514103,0.000286748,0.32,0.001005803,0.008086654,0.17
40,0.66,0.000851765,0.023423527,0.000638823,0.32,0.001231579,0.009901895,0.17
50,0.66,0.001489873,0.0394717,0.001117404,0.32,0.002305449,0.018535807,0.17
60,0.66,0.006933589,0.098113786,0.005200192,0.32,0.006754596,0.054306954,0.17
70,0.66,0.022120421,0.224965092,0.016590316,0.32,0.018720727,0.150514645,0.17
80,0.66,0.059223786,0.362002579,0.04441784,0.32,0.041408882,0.332927412,0.17
100,0.66,0.087585558,0.437927788,0.065689168,0.32,0.076818182,0.617618182,0.17")

# Source for updated Prop_critical_fatal and Prop_hospitalised_critical: https://www.bmj.com/content/369/bmj.m1985/

P.death_nonicu <- c(0.00003,
                    0.00001, 0.00001, 
                    0.00003, 0.00006, 
                    0.00013, 0.00024, 0.00040, 0.00075, 
                    0.00121, 0.00207, 0.00323, 0.00456, 
                    0.01075, 0.01674, 
                    0.03203, 
                    0.08292)
# Ensemble estimates, Table S3, p. 22/23 https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2918-0/MediaObjects/41586_2020_2918_MOESM1_ESM.pdf

# limit to 16 age groups, and replace the last group with slightly higher IFRs informed by REACT3, Table 3
P.death_nonicu <- P.death_nonicu[1:16]
P.death_nonicu[16] <- 0.116

# update prop symptomatic with ratio admisssion/cases over deaths/cases in UK on 15.07.2020
P.symp_hospitalised <- c(0.010, 0.002, 0.002, 0.003, 
                         0.005, 0.011, 0.015, 0.020, 
                         0.027, 0.044, 0.062, 0.073, 
                         0.080, 0.084, 0.109, 0.454)



reformat = function(P)
{
  # 70-74,3388.488  75-79,2442.147  80-84,1736.567  85-89,1077.555  90-94,490.577  95-99,130.083  100+,15.834
  x = c(P[1:7], weighted.mean(c(P[8], P[9]), c(3388.488 + 2442.147, 1736.567 + 1077.555 + 490.577 + 130.083 + 15.834)));
  return (rep(x, each = 2))
}

P.icu_symp     = P.symp_hospitalised * reformat(probs[, Prop_hospitalised_critical]);
P.nonicu_symp  = P.symp_hospitalised * reformat(probs[, (1 - Prop_hospitalised_critical)]);
P.death_icu    = reformat(probs[, Prop_critical_fatal]);
P.death_nonicu = reformat(probs[, Prop_noncritical_fatal]);
hfr = probs[, Prop_noncritical_fatal / Prop_symp_hospitalised]

burden_processes = list(
  list(source = "Ip", type = "multinomial", names = c("to_icu", "to_nonicu", "null"), report = c("", "", ""),
       prob = matrix(c(P.icu_symp, P.nonicu_symp, 1 - P.icu_symp - P.nonicu_symp), nrow = 3, ncol = 16, byrow = T),
       delays = matrix(c(cm_delay_gamma(7, 7, 60, 0.25)$p, cm_delay_gamma(7, 7, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 3, byrow = T)),
  
  list(source = "to_icu", type = "multinomial", names = "icu", report = "p",
       prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
       delays = matrix(cm_delay_gamma(10, 10, 60, 0.25)$p, nrow = 1, byrow = T)),
  
  list(source = "to_nonicu", type = "multinomial", names = "nonicu", report = "p",
       prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
       delays = matrix(cm_delay_gamma(8, 8, 60, 0.25)$p, nrow = 1, byrow = T)),
  
  list(source = "E", type = "multinomial", names = c("death", "null"), report = c("o", ""),
       prob = matrix(c(P.death_nonicu, 1 - P.death_nonicu), nrow = 2, ncol = 16, byrow = T),
       delays = matrix(c(cm_delay_gamma(22, 22, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),
  
  list(source = "V", type="multinomial", names="onedose", report="i",
       prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
       delays = matrix(cm_delay_gamma(8, 8, 60, 0.25)$p, nrow = 1, byrow = T)),
  
  list(source = "V2", type="multinomial", names="twodose", report="i",
       prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
       delays = matrix(cm_delay_gamma(8, 8, 60, 0.25)$p, nrow = 1, byrow = T))
   )

params$processes = burden_processes

## SET UP FOR PSA ## 

# Efficacy values are from Pouwels et al. 

# Efficacy against infection - first dose #
se.inf1 <- (0.63-0.50)/3.92
mean.inf1 <- 0.57 
infeff1.parms <- beta_mom(mean.inf1, se.inf1)

# Efficacy against disease - first dose # (efficacy against self-reported symptoms)
se.dis1 <- ((0.64-0.50)/(1-0.50))/3.92 # would be - 0 because efficacy against disease given infection can't go below zero? (I guess technically maybe it could - would mean that once infected, vaccinated ppl. were more likely to experience symptoms than unvaccinated people)
mean.dis1 <- (0.58-0.57)/(1-0.57)
diseff1.parms <- beta_mom(mean.dis1, se.dis1)

# Efficacy against infection - second dose #
se.inf2 <- (0.83-0.77)/3.92
mean.inf2 <- 0.80
infeff2.parms <- beta_mom(mean.inf2, se.inf2)

# Efficacy against disease - second dose #
se.dis2 <- ((0.86-0.77)/(1-0.77) - 0 )/3.92
mean.dis2 <- (0.84-0.8)/(1-0.8)
diseff2.parms <- beta_mom(mean.dis2, se.dis2)

# Duration of immunity - natural and vaccine-induced #
# Log-normal for immunity
a <- 365.25*(45/52)
b <- 365.25*5
m <- 365.25*3
mean.imm <- 365.25*3
nom <- a - 2*m + b
imm.parms.mean <- (a + 2*m + b)/4
imm.parms.sd <- sqrt((1/12)*(((a - 2*m + b)^2)/4 + (b-a)^2))
imm.parms <- lognorm_mom(imm.parms.mean, imm.parms.sd)


# Staff turnover rate #
mean.staff <- 0.084/365.25 # annual rate 
se.staff <- sqrt((mean.staff*(1-mean.staff))/52757)
staff.parms <- beta_mom(mean.staff, se.staff)

# Prisoner turnover rate #
mean.pris <- 0.006 # daily rate
se.pris <- 0.0017
pris.parms <- beta_mom(mean.pris, se.pris)

# QALY loss #
# per case
sym.parms <- list(alpha=1.349, beta=167.3)
# per non-icu
nonicu.parms <- beta_mom(0.018, 0.0018)
# per icu
icu.parms <- beta_mom(0.154, 0.0304)

# LFD sensitivity
mean.LFD.sens <- 0.8
se.LFD.sens <- 0.4/6
lfd.parms <- beta_mom(mean.LFD.sens, se.LFD.sens)

# LFD uptake
mean.LFD.uptake <- 0.508
se.LFD.uptake <- (0.92-0.28)/6
lfd.uptake.parms <- beta_mom(mean.LFD.uptake, se.LFD.uptake)

# Vaccine uptake
mean.vac.uptake <- 0.675
vac.parms.sd <- sqrt((1/12)*(((0.36 - 2*0.675 + 0.87)^2)/4 + (0.87-0.36)^2))
vac.parms <- beta_mom(mean.vac.uptake, vac.parms.sd)


# R0

a.R0 <- 3.2
b.R0 <- 8
mean.R0 <- 5.08
R0.parms.sd <- sqrt((1/12)*(((a.R0 - 2*mean.R0 + b.R0)^2)/4 + (b.R0-a.R0)^2))
R0.parms <- lognorm_mom(mean.R0, R0.parms.sd)

# Percentage reduction in R0 resulting from cohorting + shielding.

## SET AT BASE CASE VALUES ##
eff_inf1 <- mean.inf1
eff_inf2 <- mean.inf2
eff_dis1 <- mean.dis1
eff_dis2 <- mean.dis2
target_R0 <- 5
n_vacc_daily <- 20 # increase to 40 for prison size > 1000
imm_nat <- 1/mean.imm
imm_vac <- 1/mean.imm
staff_to <- mean.staff
b_rate <- mean.pris

## Set up for dose calculations
pop1 <- sum(params$pop[[1]]$size) * (1+(staff_to*(365.25*timehorizon)))
pop2 <- sum(params$pop[[2]]$size) * (1+(staff_to*(365.25*timehorizon)))
pop3 <- sum(params$pop[[3]]$size) * (1+(b_rate*(365.25*timehorizon)))
older.vector <- c(sum(params$pop[[1]]$size[11:16])*(1+(staff_to*(365.25*timehorizon))), sum(params$pop[[2]]$size[11:16])*(1+(staff_to*(365.25*timehorizon))), sum(params$pop[[3]]$size[11:16])*(1+(b_rate*(365.25*timehorizon))))
older <- sum(older.vector)

## Introduction via NPF staff:

## COMMUNITY TRANSMISSION ##

# Susceptibility from Davies et al. Nat Med.
u <-  c(0, 0, 0, 0.3815349, 0.7859512,
                               0.7859512, 0.8585759, 0.8585759, 0.7981468, 0.7981468,
                               0.8166960, 0.8166960, 0.8784811, 0.8784811, 0.7383189, 0.7383189)

# Number of contacts from CoMix

# June-July 2020
period2 <- c(0, # 0-4
             0, # 5-9
             0, # 10-14
             3.469002154, # 15-19
             3.469002154, # 20-24
             3.469002154, # 25-29
             4.744094276, # 30-34
             4.744094276, # 35-39
             4.160068547, # 40-44
             4.160068547, # 45-49
             3.405925893, # 50-54
             3.405925893, # 55-59
             2.815015262, # 60-64
             2.815015262, # 65-69
             2.468282967, # 70-74
             2.468282967) # 75+

# Dec 2020-Jan 2021
period7 <- c(0, # 0-4
             0, # 5-9
             0, # 10-14
             3.161434379, # 15-19
             3.161434379, # 20-24
             3.161434379, # 25-29
             3.683240043, # 30-34
             3.683240043, # 35-39
             4.160068547, # 40-44
             4.160068547, # 45-49
             3.1583651, # 50-54
             3.1583651, # 55-59
             2.8713022, # 60-64
             2.8713022, # 65-69
             2.188806371, # 70-74
             2.188806371) # 75+

# Rest of the parameters #
source(paste0(pris_path, "/scripts/sensitivity.R"))


# Vaccination before outbreak #
immune1 <- as.Date((params$time0)+28, origin=params$date0)
immune2 <- as.Date((params$time0+(7*12)+14), origin=params$date0)

