#################################
###### SUPPLEMENTARY FIG 3 ######
#################################

# Incidence of clinical cases over time (five-year time horizon)

rm(list = ls())

n_psa <- 10
n_parms <- 16
seeds <- c(123, 456, 789, 101, 999)
source("~/Documents/prisons-vacc-strategies/scripts/psa5yr.R")

for(z in seeds){
  source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")
  timehorizon <- 5
  params$time1 <- "2026-02-01"
  psa <- createpsa(z, n_psa)
  psa_temp <- psafive(n_psa, psa)
  write.csv(psa_temp[[2]], paste0(save_path,"case-5yr", z, ".csv"))
  rm(list=setdiff(ls(), c("n_psa", "n_parms", "seeds", "psa", "psafive")))
}

source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")

cases_fig10 <- data.frame()
for(i in 1:5){
  case.5yr <- read.csv(paste0(save_path,"case-5yr", seeds[i], ".csv")) %>% select(scenario_run, t, run, total)
  case.5yr<- case.5yr %>% mutate(run=run + (i-1)*n_psa)
  cases_fig10 <- rbind(cases_fig10, case.5yr)
}

# cases_fig10 <- results_psa %>% group_by(scenario_run, t, run) %>% summarise(total=sum(total))
cases_fig10 <- cases_fig10 %>% group_by(scenario_run, t) %>% summarise(mean.total=median(total), lower=quantile(total, prob=0.025), upper=quantile(total, prob=0.975))
cases_fig10 <- labelforplots(cases_fig10)
cases_fig10 <- cases_fig10 %>% mutate(t=(t-delay)) %>% filter(t>=0)
cases_fig10 <- cases_fig10 %>% mutate(t=t/365.25)

fig10 <- ggplot(cases_fig10, aes(x = t, y = mean.total, color=scenario)) +
  geom_line() + labs(x="Time in years", y="Incidence of new clinical cases") +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=scenario), alpha=0.3, color=NA) +
  facet_wrap(~scenario, ncol=2) +
  theme_bw() + theme(text=element_text(size=8)) +
  labs(color="Scenario", fill="Scenario") + scale_fill_viridis(discrete=TRUE) + scale_colour_viridis(discrete=TRUE)

ggsave(paste0(save_path,"supfig3.png"), fig10, width=170, height=150, units="mm")

# Back to base case
target_R0 <- 5
source(paste0(pris_path, "/scripts/sensitivity.R"))
params$pop[[1]]$seed_times <- round(c(seq(from=365.25, to=365.25*(timehorizon+1), by=seed_freq.NPF)), 0)
params$pop[[2]]$seed_times <- round(c(seq(from=365.25, to=365.25*(timehorizon+1), by=seed_freq.PF)), 0)
params$pop[[3]]$seed_times <- 365.25
immune1 <- ((params$time0)+28)
immune2 <- ((params$time0+(7*12)+14))
time0 <- 0
