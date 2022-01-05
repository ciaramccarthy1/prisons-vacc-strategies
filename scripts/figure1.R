#########################
####### FIGURE 1 ########
#########################

# Cases, QALYs and deaths per vaccination strategy + vaccination courses per case, QALY and death averted

rm(list = ls())

## NEED TO SET:
seeds <- c(123, 456, 789, 101, 999)
n_psa <- 100
n_parms <- 16

####################
#### BASECASE ######
####################


# Barchart with cases, QALYs lost and deaths, introducing vaccination prior to outbreak. Time horizon = 1 year #
source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")
source(paste0(pris_path, "/scripts/vacc-scenarios.R"))

## CASES over 1 year ##

cases_fig1 <- results_df %>% group_by(population, scenario_run) %>% summarize(total=sum(total))
cases_fig1 <- labelforplots(cases_fig1)

cases_fig1 = cases_fig1 %>% dplyr::filter(scenario == "(1) no vaccination") %>% group_by(scenario_run, population, scenario, scenario_nr) %>%
  summarize(total_base=total) %>% ungroup %>% dplyr::select(population, total_base) %>%
  full_join(cases_fig1, by=c("population")) %>%
  mutate(value = round((1-total/total_base)*100, 1))

cases_fig1 <- cases_fig1 %>% mutate(case.avert=total_base - total)
cases_stack <- cases_fig1 %>% group_by(scenario) %>% mutate(mean.bar=sum(total), total.value=round((100-sum(total)/sum(total_base)*100),1), overall.base=sum(total_base)) %>% ungroup()

case.stack.plot <- ggplot(cases_stack, aes(x = scenario, y = total, fill=scenario, alpha=population), show_guide=F) +
  geom_bar(stat = "identity", show.legend = T) +
  scale_alpha_discrete(range=c(0.4,1)) +
  scale_x_discrete(labels=unique(cases_stack$scenario_nr)) +
  coord_cartesian(ylim = c(0, 860)) +
  labs(x="Vaccination scenario", y="Clinical cases over one year", fill="Scenario", alpha="Population") +
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  guides(
    alpha = guide_legend(
      title = "Population",
      override.aes = aes(label = "")))

## DEATHS ##

deaths_fig1 <- results_deaths %>% group_by(run, population, scenario_run) %>% summarize(total=sum(total))
deaths_fig1 <- labelforplots(deaths_fig1)
deaths_fig1 <- deaths_fig1 %>% dplyr::filter(scenario == "(1) no vaccination") %>% group_by(scenario_run, population, scenario, scenario_nr) %>%
  summarize(total_base=total) %>% ungroup %>% dplyr::select(population, total_base) %>%
  full_join(deaths_fig1, by=c("population")) %>%
  mutate(value = round((1-total/total_base)*100, 1)) 

deaths_stack <- deaths_fig1 %>%
  group_by(scenario) %>% 
  mutate(mean.bar=sum(total), total.value=round((100-sum(total)/sum(total_base)*100),1), overall.base=sum(total_base)) %>% ungroup()

deaths.stack.plot <- ggplot(deaths_stack, aes(x = scenario, y = total, fill = scenario, alpha=population)) +
  geom_bar(stat = "identity", show.legend = T) +
  scale_alpha_discrete(range=c(0.4,1)) +
  scale_x_discrete(labels=unique(deaths_stack$scenario_nr)) +
  coord_cartesian(ylim = c(0, 15)) +
  labs(x="Vaccination scenario", y="Deaths over one year", legend.title="Scenario") +
  theme_bw() +
  theme(text=element_text(size=8.5)) +
  scale_fill_viridis(discrete=TRUE) + 
  geom_text(aes(x = scenario, y = mean.bar,
                label = ifelse(total.value==0,
                               paste0(format(total.value, nsmall = 1),"%"),
                               paste0("-",format(total.value, nsmall = 1),"%")) ),
            angle  = 0, vjust  = -0.25, hjust  = 0.5, size   = 2, lineheight = 0.8, colour = "black") +
  geom_hline(aes(yintercept=overall.base), linetype="dashed", alpha=0.5, colour = "red") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  guides(
    alpha = guide_legend(
      title = "Population",
      override.aes = aes(label = "")))

## QALYS ##

qalys_fig1 <- labelforplots(results_qalys)

qalys_fig1 <- qalys_fig1 %>% ungroup() %>%
  dplyr::filter(scenario == "(1) no vaccination") %>%
  dplyr::rename(total_base = qaly.loss) %>%
  dplyr::select(-scenario_run, -scenario, -scenario_nr) %>%
  full_join(qalys_fig1 %>% ungroup() %>%
              dplyr::select(-scenario_run),
            by=c("population")) %>%
  mutate(value = round((1-qaly.loss/total_base)*100, 1))

qaly_stack <- qalys_fig1 %>% dplyr::filter(population!="(A-C) all prisoners and staff")
qaly_stack <- qaly_stack %>% group_by(scenario) %>% mutate(mean.bar=sum(qaly.loss), total.value=round((100-sum(qaly.loss)/sum(total_base)*100),1), overall.base=sum(total_base)) %>% ungroup()

qaly.stack.plot <- ggplot(qaly_stack, aes(x = scenario, y = qaly.loss, fill = scenario, alpha=population)) +
  geom_bar(stat = "identity", show.legend = T) +
  scale_alpha_discrete(range=c(0.4,1)) +
  scale_x_discrete(labels=unique(qaly_stack$scenario_nr)) +
  labs(x="Vaccination scenario", y="QALYs lost over one year", legend.title="Scenario") +
  coord_cartesian(ylim = c(0, 200)) +
  theme_bw() +
  theme(text=element_text(size=8.5)) +
  scale_fill_viridis(discrete = TRUE) + 
  geom_text(aes(x = scenario, y = mean.bar,
                label = ifelse(total.value==0,
                               paste0(format(total.value, nsmall = 1),"%"),
                               paste0("-",format(total.value, nsmall = 1),"%")) ),
            angle  = 0, vjust  = -0.25, hjust  = 0.5, size = 2, lineheight = 0.8, colour = "black") +
  geom_hline(aes(yintercept=overall.base), linetype="dashed", alpha=0.5, colour = "red") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
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
                      labels=c("A", "B", "C"),
                      vjust=1, label_size=10)
legend_1 <- get_legend(case.stack.plot + 
                         theme(legend.position = "right", legend.box.margin = margin(0, 0, 0, 6))
)

fig1 <- plot_grid(basecase, legend_1, ncol = 2, rel_widths = c(1, 0.25))

###################################
##### DOSES WITH ERROR BARS #######
###################################

# Vaccine doses #

results_vacc <- data.frame()
for(z in 1:5){
  dose.total <- read.csv(paste0(save_path,"dose-total", seeds[z], ".csv"))
  dose.total <- dose.total %>% mutate(run=run + (z-1)*n_psa)
  results_vacc <- rbind(dose.total, results_vacc)
}

results_vacc <- labelforplots(results_vacc)

# baseline_vacc <- data.frame(
#  scenario_run=c("run1", "run2", "run3", "run4", "run5", "run6", "run7"),
#  dose=c(0, sum(params$pop[[1]]$size), 
#         sum(params$pop[[2]]$size),
#         sum(params$pop[[1]]$size)+sum(params$pop[[2]]$size),
#         sum(params$pop[[3]]$size),
#        sum(params$pop[[1]]$size[11:16])+sum(params$pop[[2]]$size[11:16])+sum(params$pop[[3]]$size[11:16]),
#         sum(params$pop[[1]]$size)+sum(params$pop[[2]]$size)+sum(params$pop[[3]]$size)
#  )
#)

# results_vacc <- results_vacc %>% group_by(scenario_run) %>% summarise(mean.vacc.count=mean(vacc.count), lower=quantile(vacc.count, prob=0.05), upper=quantile(vacc.count, prob=0.95)) 

# Cases averted #

# case.avert <- results_psa %>% filter(t>365.25) %>% group_by(scenario_run, run) %>% summarise(total=sum(total))
case.avert <- data.frame()
for(z in 1:5){
  case.total <- read.csv(paste0(save_path,"case-total", seeds[z], ".csv")) %>% select(scenario_run, run, total)
  case.total <- case.total %>% mutate(run=run + (z-1)*n_psa)
  case.avert <- rbind(case.total, case.avert)
}

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
  death.total <- read.csv(paste0(save_path,"death-total", seeds[i], ".csv")) %>% select(scenario_run, run, total)
  death.total <- death.total %>% mutate(run=run + (i-1)*n_psa)
  death.avert <- rbind(death.total, death.avert)
}

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
  qaly.total <- read.csv(paste0(save_path,"qaly-total", seeds[i], ".csv")) %>% select(scenario_run, run, total)
  qaly.total <- qaly.total %>% mutate(run=run + (i-1)*n_psa)
  qaly.avert <- rbind(qaly.total, qaly.avert)
}

qaly.avert <- qaly.avert %>% dplyr::filter(scenario_run=="run1") %>% group_by(run) %>% summarise(total_base=total) %>% ungroup %>% dplyr::select(run, total_base) %>%
  full_join(qaly.avert, by=c("run"))
qaly.avert <- qaly.avert %>% mutate(qaly.avert=total_base-total)
# qaly.avert <- qaly.avert %>% group_by(scenario_run) %>% summarise(mean.qaly.avert=mean(qaly.avert), lower.qaly=quantile(qaly.avert, prob=0.05), upper.qaly=quantile(qaly.avert, prob=0.95))
qaly.avert <- labelforplots(qaly.avert)
dose.calc <- left_join(dose.calc, qaly.avert, by=c("scenario_nr", "scenario", "scenario_run", "run"))


dose.calc <- dose.calc %>% mutate(vacc.per.caseavert= ifelse(scenario_nr=="(1)", 0, vacc.count/case.avert),
                                  vacc.per.deathavert= ifelse(scenario_nr=="(1)", 0, vacc.count/death.avert),
                                  vacc.per.qalyavert= ifelse(scenario_nr=="(1)", 0, vacc.count/qaly.avert))

dose.calc <- dose.calc %>% group_by(scenario_run) %>% summarise(mean.case=median(vacc.per.caseavert),
                                                                mean.qaly=median(vacc.per.qalyavert),
                                                                mean.death=median(vacc.per.deathavert),
                                                                lower.case=quantile(vacc.per.caseavert, prob=0.025),
                                                                upper.case=quantile(vacc.per.caseavert, prob=0.975),
                                                                lower.qaly=quantile(vacc.per.qalyavert, prob=0.025),
                                                                upper.qaly=quantile(vacc.per.qalyavert, prob=0.975),
                                                                lower.death=quantile(vacc.per.deathavert, prob=0.025),
                                                                upper.death=quantile(vacc.per.deathavert, prob=0.975))

dose.calc <- labelforplots(dose.calc)

case <- ggplot(dose.calc, aes(x=scenario_nr, y=mean.case, fill=scenario)) + geom_col(alpha = 0.7) + 
  xlab("Vaccination scenario") + ylab("Vaccinations per case averted") + scale_fill_viridis(discrete=TRUE, name="Scenario") + theme_bw() + coord_cartesian(ylim = c(0, 20)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + theme(text=element_text(size=8.5)) + geom_errorbar(mapping=aes(ymin=lower.case, ymax=upper.case), width=0.2)

death <- ggplot(dose.calc, aes(x=scenario_nr, y=mean.death, fill=scenario)) + geom_col(alpha = 0.7) + 
  xlab("Vaccination scenario") + ylab("Vaccinations per death averted") + scale_fill_viridis(discrete=TRUE, name="Scenario") + theme_bw() + coord_cartesian(ylim = c(0, 2000)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + theme(text=element_text(size=8.5)) + geom_errorbar(mapping=aes(ymin=lower.death, ymax=upper.death), width=0.2)

qaly <- ggplot(dose.calc, aes(x=scenario_nr, y=mean.qaly, fill=scenario)) + geom_col(alpha = 0.7) + 
  xlab("Vaccination scenario") + ylab("Vaccinations per QALY loss averted") + 
  scale_fill_viridis(discrete=TRUE, name="Scenario") + theme_bw() + coord_cartesian(ylim = c(0, 250)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + theme(text=element_text(size=8.5)) +
  geom_errorbar(mapping=aes(ymin=lower.qaly, ymax=upper.qaly), width=0.2)

dose.per.outcome <- plot_grid(case + theme(legend.position = "none"),
                              qaly + theme(legend.position = "none"),
                              death + theme(legend.position = "none"), ncol=1, labels=c("D", "E", "F"), vjust=1, label_size=10)

fig1a <- plot_grid(basecase, dose.per.outcome, legend_1, ncol = 3, rel_widths = c(1, 1, 0.75))
ggsave(paste0(save_path,"fig1.png"), fig1a, width=170, height=190, units="mm")

