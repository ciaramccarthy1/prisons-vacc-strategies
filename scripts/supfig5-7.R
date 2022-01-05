#####################################
##### SUPPLEMENTARY FIGURES 5-7 #####
#####################################

rm(list = ls())
# Barchart with cases, QALYs lost and deaths, introducing vaccination prior to outbreak. Time horizon = 1 year #
source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")
source(paste0(pris_path, "/scripts/vacc-scenarios.R"))

cases_facet <- results_df %>% group_by(population, scenario_run) %>% summarise(total=sum(total))
cases_facet <- rbind(cases_facet, cases_facet %>% group_by(scenario_run) %>% summarise(total=sum(total)) %>% mutate(population="(A-C) all prisoners and staff"))
cases_facet <- labelforplots(cases_facet)

cases_facet <- cases_facet %>% dplyr::filter(scenario == "(1) no vaccination") %>% group_by(scenario_run, population, scenario, scenario_nr) %>%
  summarize(total_base=total) %>% ungroup %>% dplyr::select(population, total_base) %>%
  full_join(cases_facet, by=c("population")) %>%
  mutate(value = round((1-total/total_base)*100, 1))

cases_facet <- cases_facet %>% mutate(case.avert=total_base - total)

cases_facet$population <- factor(cases_facet$population, levels=c(c("(A) non-prisoner-facing staff", "(B) prisoner-facing staff", "(C) prisoners", "(A-C) all prisoners and staff")))

# Plot by population
case_bar <- ggplot(cases_facet, aes(x = scenario, y = total, fill = scenario)) +
  geom_bar(stat = "identity", alpha = 0.7, show.legend=T) +
  facet_wrap(~population, ncol=2, scales = "free") +
  scale_x_discrete(labels=unique(cases_facet$scenario_nr)) +
  labs(x="Vaccination scenario", y="Clinical cases over one year", fill="Scenario") +
  theme_bw()+ theme(text = element_text(size = 8.5)) +
  scale_fill_viridis(discrete=T) + 
  geom_text(aes(x = scenario, y = total,
                label = ifelse(value==0,
                               paste0(format(value, nsmall = 1),"%"),
                               ifelse(value<0, paste0(format(value*-1, nsmall = 1),"%"),
                                      paste0("-",format(value, nsmall = 1),"%")))  ),
            angle  = 0, vjust  = -0.25, hjust  = 0.5, size   = 2, lineheight = 0.8, colour = "black") +  geom_hline(aes(yintercept=total_base), linetype="dashed", alpha=0.5, colour = "red") 

## DEATHS ##

# 

deaths_facet <- results_deaths %>% group_by(run, population, scenario_run) %>% summarize(total=sum(total))
deaths_facet <- rbind(deaths_facet, deaths_facet %>% group_by(scenario_run) %>% summarise(total=sum(total)) %>% mutate(population="(A-C) all prisoners and staff"))
deaths_facet <- labelforplots(deaths_facet)
deaths_facet$population <- factor(deaths_facet$population, levels=c(c("(A) non-prisoner-facing staff", "(B) prisoner-facing staff", "(C) prisoners", "(A-C) all prisoners and staff")))
deaths_facet <- deaths_facet %>% dplyr::filter(scenario == "(1) no vaccination") %>% group_by(scenario_run, population, scenario, scenario_nr) %>%
  summarize(total_base=total) %>% ungroup %>% dplyr::select(population, total_base) %>%
  full_join(deaths_facet, by=c("population")) %>%
  mutate(value = round((1-total/total_base)*100, 1)) 

death_bar <- ggplot(deaths_facet, aes(x = scenario, y = total, fill = scenario)) +
  geom_bar(stat = "identity", alpha = 0.7, show.legend=T) +
  facet_wrap(~population, ncol=2, scales = "free") +
  scale_x_discrete(labels=unique(deaths_facet$scenario_nr)) +
  labs(x="Vaccination scenario", y="Deaths over one year", fill="Scenario") +
  theme_bw()+ theme(text = element_text(size = 8.5)) +
  scale_fill_viridis(discrete=T) + 
  geom_text(aes(x = scenario, y = total,
                label = ifelse(value==0,
                               paste0(format(value, nsmall = 1),"%"),
                               ifelse(value<0, paste0(format(value*-1, nsmall = 1),"%"),
                                      paste0("-",format(value, nsmall = 1),"%")))  ),
            angle  = 0, vjust  = -0.25, hjust  = 0.5, size   = 2, lineheight = 0.8, colour = "black") +
  geom_hline(aes(yintercept=total_base), linetype="dashed", alpha=0.5, colour = "red")

# QALYS ##

qalys_facet <- labelforplots(results_qalys)
qalys_facet$population <- factor(qalys_facet$population, levels=c(c("(A) non-prisoner-facing staff", "(B) prisoner-facing staff", "(C) prisoners", "(A-C) all prisoners and staff")))
 
qalys_facet <- qalys_facet %>% ungroup() %>%
  dplyr::filter(scenario == "(1) no vaccination") %>%
  dplyr::rename(total_base = qaly.loss) %>%
  dplyr::select(-scenario_run, -scenario, -scenario_nr) %>%
  full_join(qalys_facet %>% ungroup() %>%
              dplyr::select(-scenario_run),
            by=c("population")) %>%
  mutate(value = round((1-qaly.loss/total_base)*100, 1))

qaly_bar <- ggplot(qalys_facet, aes(x = scenario, y = qaly.loss, fill = scenario)) +
  geom_bar(stat = "identity", alpha = 0.7, show.legend=T) +
  facet_wrap(~population, ncol=2, scales = "free") +
  scale_x_discrete(labels=unique(qalys_facet$scenario_nr)) +
  labs(x="Vaccination scenario", y="QALYs lost over one year", fill="Scenario") +
  scale_fill_viridis(discrete=T) + 
  theme_bw()+ theme(text = element_text(size = 8.5)) +
  geom_text(aes(x = scenario, y = qaly.loss,
                label = ifelse(value==0,
                               paste0(format(value, nsmall = 1),"%"),
                               ifelse(value<0, paste0(format(value*-1, nsmall = 1),"%"),
                                      paste0("-",format(value, nsmall = 1),"%")))  ),
            angle  = 0, vjust  = -0.25, hjust  = 0.5, size   = 2, lineheight = 0.8, colour = "black") +
  geom_hline(aes(yintercept=total_base), linetype="dashed", alpha=0.5, colour = "red")

ggsave(paste0(save_path,"supfig5.png"), case_bar, width=170, height=125, units="mm")
ggsave(paste0(save_path,"supfig6.png"), qaly_bar, width=170, height=125, units="mm")
ggsave(paste0(save_path,"supfig7.png"), death_bar, width=170, height=125, units="mm")

