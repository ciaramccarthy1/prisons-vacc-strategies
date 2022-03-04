#####################
##### FIGURE 7 ######
#####################

# Incidence of new clinical cases, with and without cohorting

rm(list = ls())

## NEED TO SET:
n_psa <- 10
n_parms <- 16
seeds <- c(123, 456, 789, 101, 999)

## Reduction in contacts due to shielding and cohorting ##
for(z in seeds){
  source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")
  psa.npi <- createpsa(z, n_psa)
  psa.npi$target_R0 <- psa.npi$target_R0*0.6
  psa_temp <- psafunc(n_psa, psa.npi)
  write.csv(psa_temp[[2]], paste0(save_path,"case-npi", z, ".csv"))
  rm(list=setdiff(ls(), c("seeds2","n_psa", "n_parms", "seeds")))
}

source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")

case_npis <- data.frame()
for(i in 1:5){
  case.npi_temp <- read.csv(paste0(save_path,"case-npi", seeds[i], ".csv")) %>% select(scenario_run, t, run, total)
  case.npi_temp <- case.npi_temp %>% mutate(run=run + (i-1)*n_psa)
  case_npis <- rbind(case_npis, case.npi_temp)
}

case_npis <- case_npis %>% group_by(scenario_run, t) %>% summarise(median.total=median(total), lower=quantile(total, prob=0.025), upper=quantile(total, prob=0.975))
case_npis <- labelforplots(case_npis)
case_npis <- case_npis %>% mutate(t=(t-delay)) %>% filter(t>=0)

### Also need output from w/o NPIs:
cases_fig5 <- data.frame()
for(i in 1:5){
  case.time <- read.csv(paste0(save_path,"case-time", seeds[i], ".csv")) %>% select(scenario_run, t, run, total)
  case.time <- case.time %>% mutate(run=run + (i-1)*n_psa)
  cases_fig5 <- rbind(cases_fig5, case.time)
}

cases_fig5 <- cases_fig5 %>% group_by(scenario_run, t) %>% summarise(median.total=median(total), lower=quantile(total, prob=0.025), upper=quantile(total, prob=0.975))
cases_fig5 <- labelforplots(cases_fig5)
cases_fig5 <- cases_fig5 %>% mutate(t=(t-delay)) %>% filter(t>=0)

## Putting them together:

case_npis <- case_npis %>% mutate(NPIs="With shielding and cohorting")
cases_fig5 <- cases_fig5 %>% mutate(NPIs="Without shielding and cohorting")

case_npis <- rbind(case_npis, cases_fig5)

fig.npis <- ggplot() +
  geom_line(case_npis, mapping=aes(x = t, y = median.total, color=NPIs)) +
  geom_ribbon(case_npis, mapping=aes(x=t, ymin=lower, ymax=upper, fill=NPIs, alpha=NPIs), color=NA) +
  facet_wrap(~scenario, ncol=2) + scale_color_manual(values =c("#0072B2","#E69F00")) + scale_fill_manual(values=c("#0072B2","#E69F00")) +
  scale_alpha_discrete(range=c(0.4,0.2)) + 
  theme_bw() + theme(text=element_text(size=8)) +
  labs(x="Time in days", y="Incidence of new clinical cases")
 

ggsave(paste0(save_path,"fig7.png"), fig.npis, width=170, height=150, units="mm")

