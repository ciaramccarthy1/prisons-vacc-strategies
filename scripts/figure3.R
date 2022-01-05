######################
###### FIGURE 3 ######
######################
rm(list = ls())


# Barchart with cases, QALYs lost and deaths, introducing vaccination prior to outbreak. Time horizon = 1 year, Varying contact patterns. #
source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")

##############################################################################
###### Scenario 1: Homogeneous staff; no contact with prisoners ##############
##############################################################################

## Mixing between populations:
params$travel[1,] <- c(0.5, 0.5, 0)
params$travel[2,] <- c(0.5, 0.5, 0)
params$travel[3,] <- c(0, 0, 1)

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
            angle  = 0, vjust  = -0.25, hjust  = 0.5, size   = 2,lineheight = 0.8, colour = "black") +
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
            angle  = 0, vjust  = -0.25, hjust  = 0.5, size = 2, lineheight = 0.8, colour = "black") +
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
scen1 <- plot_grid(case.stack.plot + theme(legend.position = "none"),
                      qaly.stack.plot + theme(legend.position = "none"),
                      deaths.stack.plot + theme(legend.position = "none"),
                      ncol=1,
                      align = 'vh')
legend_1 <- get_legend(case.stack.plot + 
                         theme(legend.position = "right", legend.box.margin = margin(0, 0, 0, 6))
)

fig1 <- plot_grid(scen1, legend_1, ncol = 2, rel_widths = c(1, 0.25))

##################################################################
#### Scenario 2: NPF none; homogenous PF and prisoners ############
###################################################################

## Mixing between populations:
params$travel[1,] <- c(1, 0, 0)
params$travel[2,] <- c(0, 0.5, 0.5)
params$travel[3,] <- c(0, 0.5, 0.5)

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
            angle  = 0, vjust  = -0.25, hjust  = 0.5, size = 2, lineheight = 0.8, colour = "black") +
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
            angle  = 0, vjust  = -0.25, hjust  = 0.5, size = 2, lineheight = 0.8, colour = "black") +
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
scen2 <- plot_grid(case.stack.plot + theme(legend.position = "none"),
                      qaly.stack.plot + theme(legend.position = "none"),
                      deaths.stack.plot + theme(legend.position = "none"),
                      ncol=1,
                      align = 'vh')

legend_7 <- get_legend(case.stack.plot + 
                         theme(legend.box.margin = margin(0, 0, 0, 6)))

fig.mix <- plot_grid(scen1, scen2, legend_1, ncol=3, rel_widths = c(1,1,0.75), labels=c("1", "2"), vjust=1, label_size=10)
ggsave(paste0(save_path,"fig3.png"), fig.mix, width=170, height=190, units="mm")
