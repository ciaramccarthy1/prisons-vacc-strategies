#################################################
####### FIGURE 2 - new variant scenario #########
#################################################

# Incidence of new clinical cases over time (one-year time horizon)
rm(list = ls())

source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")
source('~/Documents/prisons-vacc-strategies/R/psa-newvar.R')

# # Change time horizon so that I have two-year burn-in period
# 
# timehorizon <- 3
# params$date1 <- as.Date(round(timehorizon*365.25), origin=params$date0)
# params$time1 <- as.Date(round(timehorizon*365.25), origin=params$date0)
# 
# # Update vaccine introduction timing - after one year 
# 
# immune1 <- as.Date((params$time0)+365.25+28, origin=params$date0)
# immune2 <- as.Date((params$time0+(7*12)+365.25+14), origin=params$date0)

## Run analysis ##
n_psa <- 10
n_parms <- 16
seeds <- c(123, 456, 789, 101, 999)
tic()
save_path <- "C:/Users/CiaraMcCarthy/Documents/prisons-vacc-strategies/outputs/"
for(z in seeds){
  psa <- createpsa(z, n_psa)
  #for(i in 1:nr_pops){
  #  params$pop[[i]]$imm0 <- c(rep(0.5,16))
  #}
  psa_temp <- psa_newvar(n_psa, psa) # modified function that allows for reinfection and earlier seed times
  write.csv(psa_temp[[1]], paste0(save_path,"om-case-total", z, ".csv"))
  write.csv(psa_temp[[2]], paste0(save_path,"om-case-time", z, ".csv"))
  write.csv(psa_temp[[3]], paste0(save_path,"om-qaly-total", z, ".csv"))
  write.csv(psa_temp[[4]], paste0(save_path,"om-death-total", z, ".csv"))
  write.csv(psa_temp[[5]], paste0(save_path,"om-cases-prcc", z, ".csv"))
  write.csv(psa_temp[[6]], paste0(save_path,"om-qalys-prcc", z, ".csv"))
  write.csv(psa_temp[[7]], paste0(save_path,"om-dose-total", z, ".csv"))
}
toc()

cases_om <- data.frame()
for(i in 1:5){
  case.time <- read.csv(paste0(save_path,"om-case-time", seeds[i], ".csv")) %>% select(scenario_run, t, run, total)
  case.time <- case.time %>% mutate(run=run + (i-1)*n_psa)
  cases_om <- rbind(cases_om, case.time)
}

cases_om <- cases_om %>% group_by(scenario_run, t) %>% summarise(median.total=median(total), lower=quantile(total, prob=0.025), upper=quantile(total, prob=0.975))
cases_om <- labelforplots(cases_om)
cases_om <- cases_om %>% mutate(t=(t-365)) %>% filter(t>=0)

fig_om <- ggplot(cases_om, aes(x = t, y = median.total, color=scenario)) +
  geom_line() + labs(x="Time in days", y="Incidence of new clinical cases") +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=scenario), alpha=0.3, color=NA) +
  facet_wrap(~scenario, ncol=2) +
  theme_bw() + theme(text=element_text(size=8)) +
  labs(color="Scenario", fill="Scenario") + scale_fill_viridis(discrete=TRUE) + scale_colour_viridis(discrete=TRUE)

ggsave(paste0(save_path,"fig_newvar.png"), fig_om, width=170, height=150, units="mm")

######## ATTEMPT TWO #########
rm(list = ls())

source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")

n_psa <- 10
n_parms <- 16
seeds <- c(123, 456, 789, 101, 999)
tic()
save_path <- "C:/Users/CiaraMcCarthy/Documents/prisons-vacc-strategies/outputs/"
for(z in seeds){
  psa <- createpsa(z, n_psa)
  for(i in 1:nr_pops){
    params$pop[[i]]$imm0 <- c(rep(0.5,16))
  }
  psa_temp <- psafunc(n_psa, psa) # modified function that allows for reinfection and earlier seed times
  write.csv(psa_temp[[1]], paste0(save_path,"om-case-total", z, ".csv"))
  write.csv(psa_temp[[2]], paste0(save_path,"om-case-time", z, ".csv"))
  write.csv(psa_temp[[3]], paste0(save_path,"om-qaly-total", z, ".csv"))
  write.csv(psa_temp[[4]], paste0(save_path,"om-death-total", z, ".csv"))
  write.csv(psa_temp[[5]], paste0(save_path,"om-cases-prcc", z, ".csv"))
  write.csv(psa_temp[[6]], paste0(save_path,"om-qalys-prcc", z, ".csv"))
  write.csv(psa_temp[[7]], paste0(save_path,"om-dose-total", z, ".csv"))
}
toc()

cases_om <- data.frame()
for(i in 1:5){
  case.time <- read.csv(paste0(save_path,"om-case-time", seeds[i], ".csv")) %>% select(scenario_run, t, run, total)
  case.time <- case.time %>% mutate(run=run + (i-1)*n_psa)
  cases_om <- rbind(cases_om, case.time)
}

cases_om <- cases_om %>% group_by(scenario_run, t) %>% summarise(median.total=median(total), lower=quantile(total, prob=0.025), upper=quantile(total, prob=0.975))
cases_om <- labelforplots(cases_om)
cases_om <- cases_om %>% mutate(t=(t-365)) %>% filter(t>=0)

fig_om <- ggplot(cases_om, aes(x = t, y = median.total, color=scenario)) +
  geom_line() + labs(x="Time in days", y="Incidence of new clinical cases") +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=scenario), alpha=0.3, color=NA) +
  facet_wrap(~scenario, ncol=2) +
  theme_bw() + theme(text=element_text(size=8)) +
  labs(color="Scenario", fill="Scenario") + scale_fill_viridis(discrete=TRUE) + scale_colour_viridis(discrete=TRUE)

ggsave(paste0(save_path,"fig_newvar2.png"), fig_om, width=170, height=150, units="mm")


