## Vaccine parameters

## VACCINE PARAMETERS ##

for(i in 1:nr_pops){params$pop[[i]]$ei_v <- rep(mean.inf1,16) # vaccine efficacy against infection, first dose
params$pop[[i]]$ei_v2 <- rep(mean.inf2,16) # vaccine efficacy against infection, second dose
params$pop[[i]]$ed_vi <- rep(mean.dis1,16) # vaccine efficacy against disease given infection, first dose (i.e. no additional protection)
params$pop[[i]]$ed_vi2 <- rep(mean.dis2,16) # vaccine efficacY against disease given infection, second dose (0.78-0.67)/(1-0.67)
# params$pop[[i]]$v2 = rep(0,16) # nobody going straight from S to V2
}


for(i in 1:nr_pops){
params$pop[[i]]$wn <- rep(imm_nat,16)
params$pop[[i]]$wv <- rep(imm_vac,16)
params$pop[[i]]$wv2 <- rep(imm_vac,16)} # But that means the same rate of waning in people with one or two doses which doesn't really make sense

## Turnover

params$pop[[1]]$B <- staff_to * age_staff
params$pop[[2]]$B <- staff_to * age_staff
params$pop[[1]]$D <- c(rep(0,3), rep(staff_to,13))
params$pop[[2]]$D <- c(rep(0,3), rep(staff_to, 13))

params$pop[[3]]$B <- b_rate * age_prisoners # daily
params$pop[[3]]$D <- c(rep(0,3), rep(b_rate, 13))


Births <- sum(params$pop[[1]]$B*sum(params$pop[[1]]$size)) # age-stratified birth rates*(population size)
Deaths <- sum(params$pop[[1]]$size * params$pop[[1]]$D) # age-stratified population*(constant death rate)

## Scaling u based on desired R0
for(i in 1:nr_pops){
  # age-dependent ratio of infection:cases (based on Davies et al, Nature paper)
  params$pop[[i]]$y <- c(0.2904047, 0.2904047, 0.2070468, 0.2070468, 0.2676134,
                         0.2676134, 0.3284704, 0.3284704, 0.3979398, 0.3979398,
                         0.4863355, 0.4863355, 0.6306967, 0.6306967, 0.6906705, 0.6906705)
  
  # susceptibility (based on Davies et al, Nature paper)
   params$pop[[i]]$u <- c(0, 0, 0, 0.3815349, 0.7859512,
                         0.7859512, 0.8585759, 0.8585759, 0.7981468, 0.7981468,
                         0.8166960, 0.8166960, 0.8784811, 0.8784811, 0.7383189, 0.7383189)
  
  
  # scale u (susceptibility) to achieve desired R0
  current_R0 = cm_calc_R0(params, i); # calculate R0 in population i of params
  params$pop[[i]]$u = params$pop[[i]]$u * target_R0 / current_R0
  } 

## Community prevalence and contacts ##

uptake <- mean.vac.uptake

# Proportion of cases not detected through regular testing
LFD.uptake <- mean.LFD.uptake # LFD uptake 21-28 feb (COVID-19 transmission in prison settings report)
LFD.sens <- mean.LFD.sens

undetect <- 1 - (LFD.uptake*LFD.sens)

# Prevalence * number of contacts * susceptibility * proportion undetected 
prev <- (0.0003 + 0.0206)/2
prev <- 0.0047
contacts <- (period2+period7)/2

lambda <- prev*contacts*u*undetect

# lambda <- incidence*undetect

# 0.003 and 0.0206 = high and low rolling 14-day average proportion of PCR tests coming back positive -> if a prevalence, does it matter if daily or fortnightly?

import.total.NPF <- sum(lambda*params$pop[[1]]$size)
import.total.PF <- sum(lambda*params$pop[[2]]$size)

seed_freq.NPF <- 1/import.total.NPF
seed_freq.PF <- 1/import.total.PF

# Imported infections

params$pop[[1]]$seed_times <- round(c(seq(from=365.25, to=365.25+round(365.25*timehorizon), by=seed_freq.NPF)), 0)
params$pop[[2]]$seed_times <- round(c(seq(from=365.25, to=365.25+round(365.25*timehorizon), by=seed_freq.PF)), 0)
params$pop[[3]]$seed_times <- 365.25
for(i in nr_pops){params$pop[[i]]$dist_seed_ages <- c(rep(0,3), rep(1,13))}



#### VACCINATION RATE - NEW PRISONERS/STAFF
new1 <- n_vacc_daily/(sum(params$pop[[1]]$size)*staff_to*uptake)
new2 <- n_vacc_daily/(sum(params$pop[[2]]$size)*staff_to*uptake)
new3 <- n_vacc_daily/(sum(params$pop[[3]]$size)*b_rate*uptake)

prop.pop1 <- sum(params$pop[[1]]$size)/sum(params$pop[[1]]$size+params$pop[[2]]$size+params$pop[[3]]$size)
prop.pop2 <- sum(params$pop[[2]]$size)/sum(params$pop[[1]]$size+params$pop[[2]]$size+params$pop[[3]]$size)
prop.pop3 <- sum(params$pop[[3]]$size)/sum(params$pop[[1]]$size+params$pop[[2]]$size+params$pop[[3]]$size)

vacc_vals7.pop1 <- c((n_vacc_daily-5)*prop.pop1*params$pop[[1]]$size/sum(params$pop[[1]]$size))
vacc_vals7.pop2 <- c((n_vacc_daily-5)*prop.pop2*params$pop[[2]]$size/sum(params$pop[[2]]$size))
vacc_vals7.pop3 <- c((n_vacc_daily-5)*prop.pop3*params$pop[[3]]$size/sum(params$pop[[3]]$size)) + 5*params$pop[[3]]$size/sum(params$pop[[3]]$size)

new1.7 <- sum(vacc_vals7.pop1)/(sum(params$pop[[1]]$size)*staff_to*uptake)
new2.7 <- sum(vacc_vals7.pop2)/(sum(params$pop[[2]]$size)*staff_to*uptake)
new3.7 <- sum(vacc_vals7.pop3)/(sum(params$pop[[3]]$size)*b_rate*uptake)





