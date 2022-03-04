######################################
#### Extra plot - prior immunity #####
######################################


rm(list = ls())

## NEED TO SET:
n_psa <- 10
n_parms <- 16
seeds <- c(123, 456, 789, 101, 999)


for(z in seeds){
  source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")
  params$pop[[1]]$seed_times <- round(c(seq(from=delay, to=365.25*timehorizon+delay, by=seed_freq.NPF)), 0)
  params$pop[[2]]$seed_times <- round(c(seq(from=delay, to=365.25*timehorizon+delay, by=seed_freq.PF)), 0)
  params$pop[[3]]$seed_times <- delay
  immune1 <- (params$time0+28)
  immune2 <- (params$time0+(7*12)+14)
  print("check")
  psa <- createpsa(z, n_psa)
  psa_imm <- data.frame()
  for(s in c(0.1,0.3,0.5)){
    for(i in 1:nr_pops){
      params$pop[[i]]$imm0 <- rep(s,16)
    }
    psa_temp <- psafunc(n_psa, psa)
    psa_temp <- psa_temp[[2]] %>% mutate(imm=s)
    psa_imm <- rbind(psa_temp, psa_imm)
  }
  write.csv(psa_imm, paste0(save_path,"case-imm", z, ".csv"))
  rm(list=setdiff(ls(), c("seeds2","n_psa", "n_parms", "seeds")))
}

source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")

cases_simul <- data.frame()
for(i in 1:5){
  simul_temp <- read.csv(paste0(save_path,"case-imm", seeds[i], ".csv")) %>% select(scenario_run, t, run, total, imm)
  simul_temp <- simul_temp %>% mutate(run=run + (i-1)*n_psa)
  cases_simul <- rbind(cases_simul, simul_temp)
}

cases_simul <- cases_simul %>% group_by(scenario_run, t, imm) %>% summarise(mean.total=median(total), lower=quantile(total, prob=0.025), upper=quantile(total, prob=0.975))
cases_simul <- labelforplots(cases_simul)
cases_simul <- cases_simul %>% mutate(t=(t-delay)) %>% filter(t>=0)
cases_simul$imm <- as.factor(cases_simul$imm)

fig9 <- ggplot(cases_simul, aes(x = t, y = mean.total, color=scenario, alpha=imm)) +
  geom_line() + labs(x="Time in days", y="Incidence of new clinical cases") +
  # geom_ribbon(aes(ymin=lower, ymax=upper, fill=scenario), alpha=0.3, color=NA) +
  facet_wrap(~scenario, ncol=2) +
  scale_alpha_discrete(range=c(0.4,1)) +
  theme_bw() + theme(text=element_text(size=8)) +
  labs(color="Scenario", fill="Scenario", alpha="Proportion with prior immunity") + scale_fill_viridis(discrete=TRUE) + scale_colour_viridis(discrete=TRUE)

ggsave(paste0(save_path,"supfig-imm.png"), fig9, width=170, height=150, units="mm")

# Back to base case:
immune1 <- as.Date((params$time0)+28, origin=params$date0)
immune2 <- as.Date((params$time0+(7*12)+14), origin=params$date0)