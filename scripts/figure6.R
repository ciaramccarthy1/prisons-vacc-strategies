#######################
###### FIGURE 6 #######
#######################

rm(list = ls())
# Deterministic sensitivity analysis

source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")

##### Prisoner turnover rate ######
prisoner_to_vector <- c(0, mean.pris, mean.pris*2)
names(prisoner_to_vector) <- c("lower", "base", "upper")
results_prisoner <- data.frame()
for(j in prisoner_to_vector){
  b_rate <- j
  source(paste0(pris_path, "/scripts/sensitivity.R"))
  source(paste0(pris_path, "/scripts/vacc-scenarios.R"))
  results_temp <- results_df %>% group_by(run, population, scenario_run) %>% summarize(total=sum(total)) %>% labelforplots()
  results_temp <- results_temp %>% filter(scenario == "(1) no vaccination") %>% group_by(scenario_run, scenario, scenario_nr, population) %>%
    summarize(total_base=total) %>% ungroup %>% select(population, total_base) %>%
    full_join(results_temp, by="population")
  results_temp <- results_temp %>% group_by(scenario_run, scenario) %>% summarise(total_base=sum(total_base), total=sum(total))
  results_temp <- results_temp %>% mutate(value = round((1-total/total_base)*100, 1), case.avert=total_base - total)
  results_prisoner <- rbind(results_prisoner, results_temp %>% mutate(prisoner.turnover=b_rate, value=names(prisoner_to_vector)[prisoner_to_vector==j]))}
results_prisoner <- results_prisoner %>% mutate(daily.pris=prisoner.turnover*sum(params$pop[[3]]$size))
prisoner_plot <- ggplot(results_prisoner, aes(x=daily.pris, y=case.avert, color=scenario)) + geom_line() + geom_point() + labs(x="New residents per day", y="Cases averted over one year", color="Scenario") + theme_bw() + scale_color_viridis(discrete=TRUE) +theme(text=element_text(size=8))

# Back to base value
b_rate <- mean.pris
params$pop[[3]]$B <- b_rate * age_prisoners # daily
params$pop[[3]]$D <- c(rep(0,3), rep(b_rate, 13))

source(paste0(pris_path, "/scripts/sensitivity.R"))

##### Duration of vaccine immunity #####

imm_vac <- log((mean.inf2-mean.imm.vac)/mean.inf2)/-140
low_imm_vac <- log((mean.inf2-0.34)/mean.inf2)/-140
high_imm_vac <- log((mean.inf2-0.08)/mean.inf2)/-140

imm_vector <- c(low_imm_vac, imm_vac, high_imm_vac)
names(imm_vector) <- c("lower", "base", "upper")
results_imm_vac <- data.frame()
for(j in imm_vector){
  imm_value <- j
  params$pop[[1]]$wv <- rep(imm_value,16)
  params$pop[[2]]$wv <- rep(imm_value,16)
  params$pop[[3]]$wv <- rep(imm_value,16)
  params$pop[[1]]$wv2 <- rep(imm_value,16)
  params$pop[[2]]$wv2 <- rep(imm_value,16)
  params$pop[[3]]$wv2 <- rep(imm_value,16)
  source(paste0(pris_path, "/scripts/vacc-scenarios.R"))
  results_temp <- results_df %>% group_by(run, population, scenario_run) %>% summarize(total=sum(total)) %>% labelforplots()
  results_temp <- results_temp %>% filter(scenario == "(1) no vaccination") %>% group_by(scenario_run, scenario, scenario_nr, population) %>%
    summarize(total_base=total) %>% ungroup %>% select(population, total_base) %>%
    full_join(results_temp, by="population")
  results_temp <- results_temp %>% group_by(scenario_run, scenario) %>% summarise(total_base=sum(total_base), total=sum(total))
  results_temp <- results_temp %>% mutate(case.avert=total_base - total)
  results_imm_vac <- rbind(results_imm_vac, results_temp %>% mutate(imm.vac=j, value=names(imm_vector)[imm_vector==j]))}
results_imm_vac <- results_imm_vac %>% mutate(loss=-(mean.inf2*exp(-imm.vac*140) - mean.inf2))
vac_imm_plot <- ggplot(results_imm_vac, aes(x=loss*100, y=case.avert, color=scenario)) + geom_line() + geom_point() + labs(x="Loss of vaccine immunity over 5 months (%)", y="Cases averted over one year") + theme_bw() + scale_color_viridis(discrete=TRUE)  + theme(text=element_text(size=8))

# Back to base value 
for(i in 1:nr_pops){
  params$pop[[i]]$wv <- rep(imm_vac, 16)
  params$pop[[i]]$wv2 <- rep(imm_vac, 16)}

source(paste0(pris_path, "/scripts/sensitivity.R"))

## DURATION OF NATURAL IMMUNITY ##
results_imm_nat <- data.frame()

imm_nat <- log(1-mean.imm.nat)/-365.25
low_imm_nat <- log(1-0.8)/-365.25
high_imm_nat <- log(1)/-365.25

nat_vector <- c(low_imm_nat, imm_nat, high_imm_nat)

for(j in nat_vector){
  params$pop[[1]]$wn <- rep(j, 16)
  params$pop[[2]]$wn <- rep(j, 16)
  params$pop[[3]]$wn <- rep(j, 16)
  source(paste0(pris_path, "/scripts/vacc-scenarios.R"))
  results_temp <- results_df %>% group_by(run, population, scenario_run) %>% summarize(total=sum(total)) %>% labelforplots()
  results_temp <- results_temp %>% filter(scenario == "(1) no vaccination") %>% group_by(scenario_run, scenario, scenario_nr, population) %>%
    summarize(total_base=total) %>% ungroup %>% select(population, total_base) %>%
    full_join(results_temp, by="population")
  results_temp <- results_temp %>% group_by(scenario_run, scenario) %>% summarise(total_base=sum(total_base), total=sum(total))
  results_temp <- results_temp %>% mutate(case.avert=total_base - total)
  results_imm_nat <- rbind(results_imm_nat, results_temp %>% mutate(imm.nat=j, value=names(nat_vector)[nat_vector==j]))}

results_imm_nat <- results_imm_nat %>% mutate(loss=-(exp(-imm.nat*365.25) - 1))
nat_imm_plot <- ggplot(results_imm_nat, aes(x=loss*100, y=case.avert, color=scenario)) + geom_line() + geom_point() + labs(x="Loss of natural immunity over one year (%)", y="Cases averted over one year") + theme_bw() + scale_color_viridis(discrete=TRUE)  + theme(text=element_text(size=8))

# Back to base value 
for(i in 1:nr_pops){
  params$pop[[i]]$wn <- rep(imm_nat,16)}

source(paste0(pris_path, "/scripts/sensitivity.R"))

### Vaccine efficacy against disease - second dose ###

eff_vector <- c(0.5, mean.dis2, 0.95)
names(eff_vector) <- c("lower", "base", "upper")
results_eff <- data.frame()
for(j in eff_vector){
  params$pop[[1]]$ed_vi2 <- rep(VEdis(j, mean.inf2),16) 
  params$pop[[2]]$ed_vi2 <- rep(VEdis(j, mean.inf2),16)  
  params$pop[[3]]$ed_vi2 <- rep(VEdis(j, mean.inf2),16) 
  source(paste0(pris_path, "/scripts/vacc-scenarios.R"))
  results_temp <- results_df %>% group_by(run, population, scenario_run) %>% summarize(total=sum(total)) %>% labelforplots()
  results_temp <- results_temp %>% filter(scenario == "(1) no vaccination") %>% group_by(scenario_run, scenario, scenario_nr, population) %>%
    summarize(total_base=total) %>% ungroup %>% select(population, total_base) %>%
    full_join(results_temp, by="population")
  results_temp <- results_temp %>% group_by(scenario_run, scenario) %>% summarise(total_base=sum(total_base), total=sum(total))
  results_temp <- results_temp %>% mutate(case.avert=total_base - total)
  results_eff <- rbind(results_eff, results_temp %>% mutate(eff=j, value=names(eff_vector)[eff_vector==j]))}
eff_plot <- ggplot(results_eff, aes(x=eff, y=case.avert, color=scenario)) + geom_line() + geom_point() + labs(x="VE against disease (2nd dose)", y="Cases averted over one year") + theme_bw() + scale_color_viridis(discrete=TRUE) + theme(text=element_text(size=8))

# Back to base value 
for(i in 1:nr_pops){
  params$pop[[i]]$ed_vi2 <- rep(VEdis(mean.dis2, mean.inf2),16) 
}

source(paste0(pris_path, "/scripts/sensitivity.R"))


sens <- plot_grid(prisoner_plot + theme(legend.position = "none"),
                  eff_plot + theme(legend.position = "none"),
                  vac_imm_plot + theme(legend.position = "none"), 
                  nat_imm_plot + theme(legend.position = "none"),
                  ncol=2, align = 'vh')
legend_sens <- get_legend(prisoner_plot + 
                            theme(legend.position = "right", text=element_text(size=8)))

fig.sens <- plot_grid(sens, legend_sens, ncol = 2, rel_widths = c(1, 0.4))

ggsave(paste0(save_path,"fig6.png"), fig.sens, width=170, height=110, units="mm")
print(fig.sens)
