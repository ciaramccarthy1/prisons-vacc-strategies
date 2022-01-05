##########################
####### FIGURE 4 #########
##########################

rm(list = ls())
# Increasing IFR in prisoners by shifting age-specific IFRs by 10 years

## NEED TO SET:

n_psa <- 100
n_parms <- 16
seeds <- c(123, 456, 789, 101, 999)

source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")

# Prep burden_processes.old
P.icu_symp     = c(P.icu_symp[3:16],rep(P.icu_symp[16],2));
P.nonicu_symp  = c(P.nonicu_symp[3:16],rep(P.nonicu_symp[16],2));
P.death_icu    = c(P.death_icu[3:16], rep(P.death_icu[16],2));
P.death_nonicu = c(P.death_nonicu[3:16], rep(P.death_nonicu[16],2));
hfr = c(hfr[3:16],rep(hfr[16],2))

burden_processes.old = list(
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

for(z in seeds){
  source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")
  psa <- createpsa(z, n_psa)
  psa_temp <- psa.age(n_psa, psa)
  write.csv(psa_temp[[1]], paste0(save_path,"case-staff", z, ".csv"))
  write.csv(psa_temp[[2]], paste0(save_path,"qaly-staff", z, ".csv"))
  write.csv(psa_temp[[3]], paste0(save_path,"death-staff", z, ".csv"))
  write.csv(psa_temp[[7]], paste0(save_path,"dose-staff", z, ".csv"))
  params$processes = burden_processes.old
  psa_temp <- psa.age(n_psa, psa)
  write.csv(psa_temp[[4]], paste0(save_path,"case-prisoners", z, ".csv"))
  write.csv(psa_temp[[5]], paste0(save_path,"qaly-prisoners", z, ".csv"))
  write.csv(psa_temp[[6]], paste0(save_path,"death-prisoners", z, ".csv"))
  rm(list=setdiff(ls(), c("n_psa", "n_parms", "seeds", "psa", "burden_processes.old")))
}

source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")
n_psa <- 100
n_parms <- 16
seeds <- c(123, 456, 789, 101, 999)

results_vacc <- data.frame()
for(z in 1:5){
  dose.total <- read.csv(paste0(save_path,"dose-staff", seeds[z], ".csv")) %>% select(scenario_run, run, vacc.count)
  dose.total <- dose.total %>% mutate(run=run + (z-1)*n_psa)
  results_vacc <- rbind(dose.total, results_vacc)
}

results_vacc <- labelforplots(results_vacc)

case.staff <- data.frame()
for(z in 1:5){
  case.staff_temp <- read.csv(paste0(save_path,"case-staff", seeds[z], ".csv")) %>% select(scenario_run, run, total, population)
  case.staff_temp <- case.staff_temp %>% mutate(run=run + (z-1)*n_psa)
  case.staff <- rbind(case.staff_temp, case.staff)
}

case.pris <- data.frame()
for(z in 1:5){
  case.pris_temp <- read.csv(paste0(save_path,"case-prisoners", seeds[z], ".csv")) %>% select(scenario_run, run, total, population)
  case.pris_temp <- case.pris_temp %>% mutate(run=run + (z-1)*n_psa)
  case.pris <- rbind(case.pris_temp, case.pris)
}

case.avert <- rbind(case.staff, case.pris) %>% group_by(run, scenario_run) %>% summarise(total=sum(total))

case.avert <- case.avert %>% dplyr::filter(scenario_run=="run1") %>% group_by(run) %>% summarise(total_base=total) %>% ungroup %>% dplyr::select(run, total_base) %>%
  full_join(case.avert, by=c("run"))
case.avert <- case.avert %>% mutate(case.avert=total_base-total)
# case.avert <- case.avert %>% group_by(scenario_run) %>% summarise(mean.case.avert=mean(case.avert), lower.case=quantile(case.avert, prob=0.05), upper.case=quantile(case.avert, prob=0.95))
case.avert <- labelforplots(case.avert)

dose.calc <- left_join(results_vacc, case.avert, by=c("scenario_nr", "scenario", "scenario_run", "run"))

# Deaths averted #

#death.avert <- deaths_psa %>% filter(t>365.25) %>% group_by(scenario_run, run) %>% summarise(total=sum(total))

death.avert <- data.frame()
for(i in 1:5){
  death.staff <- read.csv(paste0(save_path,"death-staff", seeds[i], ".csv")) %>% select(scenario_run, run, total, population)
  death.staff <- death.staff %>% mutate(run=run + (i-1)*n_psa)
  death.pris <- read.csv(paste0(save_path,"death-prisoners", seeds[i], ".csv")) %>% select(scenario_run, run, total, population)
  death.pris <- death.pris %>% mutate(run=run + (i-1)*n_psa)
  death.total <- rbind(death.staff, death.pris)
  death.avert <- rbind(death.total, death.avert)
}
death.avert <- death.avert %>% group_by(run, scenario_run) %>% summarise(total=sum(total))

death.avert <- death.avert %>% dplyr::filter(scenario_run=="run1") %>% group_by(run) %>% summarise(total_base=total) %>% ungroup %>% dplyr::select(run, total_base) %>%
  full_join(death.avert, by=c("run"))
death.avert <- death.avert %>% mutate(death.avert=total_base-total)
# death.avert <- death.avert %>% group_by(scenario_run) %>% summarise(mean.death.avert=mean(death.avert), lower.death=quantile(death.avert, prob=0.05), upper.death=quantile(death.avert, prob=0.95))
death.avert <- labelforplots(death.avert)

dose.calc <- left_join(dose.calc, death.avert, by=c("scenario_nr", "scenario", "scenario_run", "run"))

# QALYs averted #

# qaly.avert <- qaly_psa %>% group_by(scenario_run, run) %>% summarise(total=sum(qaly.loss))
qaly.avert <- data.frame()
for(i in 1:5){
  qaly.staff <- read.csv(paste0(save_path,"qaly-staff", seeds[i], ".csv")) %>% select(scenario_run, run, total, population)
  qaly.staff <- qaly.staff %>% mutate(run=run + (i-1)*n_psa)
  qaly.pris <- read.csv(paste0(save_path,"qaly-prisoners", seeds[i], ".csv")) %>% select(scenario_run, run, total, population)
  qaly.pris <- qaly.pris %>% mutate(run=run + (i-1)*n_psa)
  qaly.total <- rbind(qaly.staff, qaly.pris)
  qaly.avert <- rbind(qaly.total, qaly.avert)
}

qaly.avert <- qaly.avert %>% group_by(run, scenario_run) %>% summarise(total=sum(total))

qaly.avert <- qaly.avert %>% dplyr::filter(scenario_run=="run1") %>% group_by(run) %>% summarise(total_base=total) %>% ungroup %>% dplyr::select(run, total_base) %>%
  full_join(qaly.avert, by=c("run"))
qaly.avert <- qaly.avert %>% mutate(qaly.avert=total_base-total)
# qaly.avert <- qaly.avert %>% group_by(scenario_run) %>% summarise(mean.qaly.avert=mean(qaly.avert), lower.qaly=quantile(qaly.avert, prob=0.05), upper.qaly=quantile(qaly.avert, prob=0.95))
qaly.avert <- labelforplots(qaly.avert)
dose.calc <- left_join(dose.calc, qaly.avert, by=c("scenario_nr", "scenario", "scenario_run", "run"))

dose.calc <- dose.calc %>% mutate(vacc.per.caseavert= ifelse(scenario_nr=="(1)", 0, vacc.count/case.avert),
                                  vacc.per.deathavert= ifelse(scenario_nr=="(1)", 0, vacc.count/death.avert),
                                  vacc.per.qalyavert= ifelse(scenario_nr=="(1)", 0, vacc.count/qaly.avert))

dose.test <- dose.calc %>% filter(scenario_nr=="(3)"| scenario_nr=="(4)")

dose.calc <- dose.calc %>% group_by(scenario_run) %>% summarise(median.case=median(vacc.per.caseavert),
                                                                median.qaly=median(vacc.per.qalyavert),
                                                                median.death=median(vacc.per.deathavert),
                                                                lower.case=quantile(vacc.per.caseavert, prob=0.025),
                                                                upper.case=quantile(vacc.per.caseavert, prob=0.975),
                                                                lower.qaly=quantile(vacc.per.qalyavert, prob=0.025),
                                                                upper.qaly=quantile(vacc.per.qalyavert, prob=0.975),
                                                                lower.death=quantile(vacc.per.deathavert, prob=0.025),
                                                                upper.death=quantile(vacc.per.deathavert, prob=0.975))

dose.calc <- labelforplots(dose.calc)

case <- ggplot(dose.calc, aes(x=scenario_nr, y=median.case, fill=scenario)) + geom_col(alpha = 0.7) + 
  xlab("Vaccination scenario") + ylab("Vaccinations per case averted") + scale_fill_viridis(discrete=TRUE, name="Scenario") + theme_bw() + coord_cartesian(ylim = c(0, 20)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + theme(text=element_text(size=8.5)) + geom_errorbar(mapping=aes(ymin=lower.case, ymax=upper.case), width=0.2)

death <- ggplot(dose.calc, aes(x=scenario_nr, y=median.death, fill=scenario)) + geom_col(alpha = 0.7) + 
  xlab("Vaccination scenario") + ylab("Vaccinations per death averted") + scale_fill_viridis(discrete=TRUE, name="Scenario") + theme_bw() + coord_cartesian(ylim = c(0, 2000)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + theme(text=element_text(size=8.5)) + geom_errorbar(mapping=aes(ymin=lower.death, ymax=upper.death), width=0.2)

qaly <- ggplot(dose.calc, aes(x=scenario_nr, y=median.qaly, fill=scenario)) + geom_col(alpha = 0.7) + 
  xlab("Vaccination scenario") + ylab("Vaccinations per QALY loss averted") + 
  scale_fill_viridis(discrete=TRUE, name="Scenario") + theme_bw() + coord_cartesian(ylim = c(0, 250)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + theme(text=element_text(size=8.5)) +
  geom_errorbar(mapping=aes(ymin=lower.qaly, ymax=upper.qaly), width=0.2)

dose.per.outcome <- plot_grid(case + theme(legend.position = "none"),
                              qaly + theme(legend.position = "none"),
                              death + theme(legend.position = "none"), ncol=1, labels=c("D", "E", "F"), vjust=1, label_size=10)

### DETERMINISTIC RESULTS ##
source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")
# Take previous results for staff:
# Deterministic:
source(paste0(pris_path, "/scripts/vacc-scenarios.R"))
results_df_staff <- results_df %>% filter(population=="(A) non-prisoner-facing staff" | population=="(B) prisoner-facing staff")
results_deaths_staff <- results_deaths %>% filter(population=="(A) non-prisoner-facing staff" | population=="(B) prisoner-facing staff")
results_qalys_staff <- results_qalys %>% filter(population=="(A) non-prisoner-facing staff" | population=="(B) prisoner-facing staff")

# But new results for prisoners:

params$processes <- burden_processes.old

source(paste0(pris_path, "/scripts/vacc-scenarios.R"))
results_df_pris <- results_df %>% filter(population=="(C) prisoners")
results_deaths_pris <- results_deaths %>% filter(population=="(C) prisoners")
results_qalys_pris <- results_qalys %>%  filter(population=="(C) prisoners")

results_df <- rbind(results_df_staff, results_df_pris)
results_deaths <- rbind(results_deaths_staff, results_deaths_pris)
results_qalys <- rbind(results_qalys_staff, results_qalys_pris)

## CASES over 1 year ##

cases_fig7 <- results_df %>% group_by(population, scenario_run) %>% summarize(total=sum(total))
cases_fig7 <- labelforplots(cases_fig7)

cases_fig7 = cases_fig7 %>% dplyr::filter(scenario == "(1) no vaccination") %>% group_by(scenario_run, population, scenario, scenario_nr) %>%
  summarize(total_base=total) %>% ungroup %>% dplyr::select(population, total_base) %>%
  full_join(cases_fig7, by=c("population")) %>%
  mutate(value = round((1-total/total_base)*100, 1))

cases_fig7 <- cases_fig7 %>% mutate(case.avert=total_base - total)
cases_stack <- cases_fig7 %>% group_by(scenario) %>% mutate(mean.bar=sum(total), total.value=round((100-sum(total)/sum(total_base)*100),1), overall.base=sum(total_base)) %>% ungroup()

case.stack.plot <- ggplot(cases_stack, aes(x = scenario, y = total, fill=scenario, alpha=population)) +
  geom_bar(stat = "identity", show.legend = T) +
  scale_alpha_discrete(range=c(0.4,1)) +
  coord_cartesian(ylim = c(0, 860)) +
  scale_x_discrete(labels=unique(cases_stack$scenario_nr)) +
  labs(x="Vaccination scenario", y="Clinical cases over one year", legend.title="Scenario") +
  theme_bw() +
  scale_color_OkabeIto() +  
  theme(text=element_text(size=8.5)) +
  scale_fill_viridis(discrete=TRUE) +
  geom_text(aes(x = scenario, y = mean.bar,
                label = ifelse(total.value==0,
                               paste0(format(total.value, nsmall = 1),"%"),
                               paste0("-",format(total.value, nsmall = 1),"%")) ),
            angle  = 0, vjust  = -0.25, hjust  = 0.5, size   = 2, lineheight = 0.8, colour = "black") +
  geom_hline(aes(yintercept=overall.base), linetype="dashed", alpha=0.5, colour = "red") +
  guides(
    alpha = guide_legend(
      title = "Population",
      override.aes = aes(label = "")))
# + geom_errorbar(aes(ymin=mean.bar-error.bar, ymax = mean.bar+error.bar), width=0.2)

## DEATHS ##

deaths_fig7 <- results_deaths %>% group_by(run, population, scenario_run) %>% summarize(total=sum(total))
deaths_fig7 <- labelforplots(deaths_fig7)
deaths_fig7 <- deaths_fig7 %>% dplyr::filter(scenario == "(1) no vaccination") %>% group_by(scenario_run, population, scenario, scenario_nr) %>%
  summarize(total_base=total) %>% ungroup %>% dplyr::select(population, total_base) %>%
  full_join(deaths_fig7, by=c("population")) %>%
  mutate(value = round((1-total/total_base)*100, 1))

deaths_stack <- deaths_fig7 %>%
  group_by(scenario) %>% 
  mutate(mean.bar=sum(total), total.value=round((100-sum(total)/sum(total_base)*100),1), overall.base=sum(total_base)) %>% ungroup()

deaths.stack.plot <- ggplot(deaths_stack, aes(x = scenario, y = total, fill = scenario, alpha=population)) +
  geom_bar(stat = "identity", show.legend = T) +
  scale_alpha_discrete(range=c(0.4,1)) +
  coord_cartesian(ylim = c(0, 15)) +
  scale_x_discrete(labels=unique(deaths_stack$scenario_nr)) +
  labs(x="Vaccination scenario", y="Deaths over one year", legend.title="Scenario") +
  theme_bw() +
  theme(text=element_text(size=8.5)) +
  scale_fill_viridis(discrete=TRUE) + 
  geom_text(aes(x = scenario, y = mean.bar,
                label = ifelse(total.value==0,
                               paste0(format(total.value, nsmall = 1),"%"),
                               paste0("-",format(total.value, nsmall = 1),"%")) ),
            angle  = 0, vjust  = -0.25, hjust  = 0.5, size   = 2, lineheight = 0.8, colour = "black") +
  geom_hline(aes(yintercept=overall.base), linetype="dashed", alpha=0.5, colour = "red")+
  guides(
    alpha = guide_legend(
      title = "Population",
      override.aes = aes(label = "")))

## QALYS ##

qalys_fig7 <- labelforplots(results_qalys)

qalys_fig7 <- qalys_fig7 %>% ungroup() %>%
  dplyr::filter(scenario == "(1) no vaccination") %>%
  dplyr::rename(total_base = qaly.loss) %>%
  dplyr::select(-scenario_run, -scenario, -scenario_nr) %>%
  full_join(qalys_fig7 %>% ungroup() %>%
              dplyr::select(-scenario_run),
            by=c("population")) %>%
  mutate(value = round((1-qaly.loss/total_base)*100, 1))

qaly_stack <- qalys_fig7 %>% dplyr::filter(population!="(A-C) all prisoners and staff")
qaly_stack <- qaly_stack %>% group_by(scenario) %>% mutate(mean.bar=sum(qaly.loss), total.value=round((100-sum(qaly.loss)/sum(total_base)*100),1), overall.base=sum(total_base)) %>% ungroup()

qaly.stack.plot <- ggplot(qaly_stack, aes(x = scenario, y = qaly.loss, fill = scenario, alpha=population)) +
  geom_bar(stat = "identity", show.legend = T) +
  scale_alpha_discrete(range=c(0.4,1)) +
  coord_cartesian(ylim = c(0, 200)) +
  scale_x_discrete(labels=unique(qaly_stack$scenario_nr)) +
  labs(x="Vaccination scenario", y="QALYs lost over one year", legend.title="Scenario") +
  theme_bw() +
  theme(text=element_text(size=8.5)) +
  scale_fill_viridis(discrete = TRUE) + 
  geom_text(aes(x = scenario, y = mean.bar,
                label = ifelse(total.value==0,
                               paste0(format(total.value, nsmall = 1),"%"),
                               paste0("-",format(total.value, nsmall = 1),"%")) ),
            angle  = 0, vjust  = -0.25, hjust  = 0.5, size   = 2, lineheight = 0.8, colour = "black") +
  geom_hline(aes(yintercept=overall.base), linetype="dashed", alpha=0.5, colour = "red")+
  guides(
    alpha = guide_legend(
      title = "Population",
      override.aes = aes(label = "")))

# Stacked bar chart #
basecase <- plot_grid(case.stack.plot + theme(legend.position = "none"),
                      qaly.stack.plot + theme(legend.position = "none"),
                      deaths.stack.plot + theme(legend.position = "none"),
                      ncol=1,
                      align = 'vh',
                      labels = c("A", "B", "C"), label_size=10)
legend_7 <- get_legend(case.stack.plot + 
                         theme(legend.position = "right"))

fig7 <- plot_grid(basecase, legend_7, ncol = 2, rel_widths = c(1, 0.25))

## Doses per case, qaly lost and death averted, introduction pre-outbreak. Time horizon = 1 year ##


fig7a <- plot_grid(basecase, dose.per.outcome, legend_7, ncol = 3, rel_widths = c(1, 1, 0.75))
ggsave(paste0(save_path,"fig4.png"), fig7a, width=170, height=190, units="mm")

# Back to base case:
params$processes = burden_processes

