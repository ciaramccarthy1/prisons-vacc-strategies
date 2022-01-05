##########################
####### FIGURE 2 #########
##########################

# Incidence of new clinical cases over time (one-year time horizon)
rm(list = ls())

## NEED TO SET:
n_psa <- 100
n_parms <- 16
seeds <- c(123, 456, 789, 101, 999)
source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")

cases_fig5 <- data.frame()
for(i in 1:5){
  case.time <- read.csv(paste0(save_path,"case-time", seeds[i], ".csv")) %>% select(scenario_run, t, run, total)
  case.time <- case.time %>% mutate(run=run + (i-1)*n_psa)
  cases_fig5 <- rbind(cases_fig5, case.time)
}

cases_fig5 <- cases_fig5 %>% group_by(scenario_run, t) %>% summarise(median.total=median(total), lower=quantile(total, prob=0.025), upper=quantile(total, prob=0.975))
cases_fig5 <- labelforplots(cases_fig5)
cases_fig5 <- cases_fig5 %>% mutate(t=(t-365)) %>% filter(t>=0)

fig5 <- ggplot(cases_fig5, aes(x = t, y = median.total, color=scenario)) +
  geom_line() + labs(x="Time in days", y="Incidence of new clinical cases") +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=scenario), alpha=0.3, color=NA) +
  facet_wrap(~scenario, ncol=2) +
  theme_bw() + theme(text=element_text(size=8)) +
  labs(color="Scenario", fill="Scenario") + scale_fill_viridis(discrete=TRUE) + scale_colour_viridis(discrete=TRUE)

ggsave(paste0(save_path,"fig2.png"), fig5, width=170, height=150, units="mm")

