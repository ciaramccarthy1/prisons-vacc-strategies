##################################
###### SUPPLEMENTARY FIG 1 #######
##################################

## Varying timing - with PSA

rm(list = ls())

## NEED TO SET:
n_psa <- 100
n_parms <- 16
seeds <- c(123, 456, 789, 101, 999)

for(z in seeds){
  source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")
  source("~/Documents/prisons-vacc-strategies/R/psa-simul.R")
  params$pop[[1]]$seed_times <- round(c(seq(from=365.25, to=365.25*(timehorizon+1), by=seed_freq.NPF)), 0)
  params$pop[[2]]$seed_times <- round(c(seq(from=365.25, to=365.25*(timehorizon+1), by=seed_freq.PF)), 0)
  params$pop[[3]]$seed_times <- 365.25
  immune1 <- (params$time0+28+365.25)
  immune2 <- (params$time0+365.25+(7*12)+14)
  print("check")
  psa <- createpsa(z, n_psa)
  psa_temp <- psa.simul(n_psa, psa)
  write.csv(psa_temp, paste0(save_path,"case-simul", z, ".csv"))
  rm(list=setdiff(ls(), c("seeds2","n_psa", "n_parms", "seeds")))
}

source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")

cases_simul <- data.frame()
for(i in 1:5){
  simul_temp <- read.csv(paste0(save_path,"case-simul", seeds[i], ".csv")) %>% select(scenario_run, t, run, total)
  simul_temp <- simul_temp %>% mutate(run=run + (i-1)*n_psa)
  cases_simul <- rbind(cases_simul, simul_temp)
}

cases_simul <- cases_simul %>% group_by(scenario_run, t) %>% summarise(mean.total=median(total), lower=quantile(total, prob=0.025), upper=quantile(total, prob=0.975))
cases_simul <- labelforplots(cases_simul)
cases_simul <- cases_simul %>% mutate(t=(t-365)) %>% filter(t>=0)

fig9 <- ggplot(cases_simul, aes(x = t, y = mean.total, color=scenario)) +
  geom_line() + labs(x="Time in days", y="Incidence of new clinical cases") +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=scenario), alpha=0.3, color=NA) +
  facet_wrap(~scenario, ncol=2) +
  theme_bw() + theme(text=element_text(size=8)) +
  labs(color="Scenario", fill="Scenario") + scale_fill_viridis(discrete=TRUE) + scale_colour_viridis(discrete=TRUE)

ggsave(paste0(save_path,"supfig1.png"), fig9, width=170, height=150, units="mm")

# Back to base case:
immune1 <- as.Date((params$time0)+28, origin=params$date0)
immune2 <- as.Date((params$time0+(7*12)+14), origin=params$date0)
