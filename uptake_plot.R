rm(list = ls())

## NEED TO SET:
seeds <- c(123, 456, 789, 101, 999)
n_psa <- 10
n_parms <- 16

####################
#### BASECASE ######
####################

# Barchart with cases, QALYs lost and deaths, introducing vaccination prior to outbreak. Time horizon = 1 year #
source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")
params$pop[[3]]$seed_times <- 99999

cases_uptake <- data.frame()
uptake.staff <- 0.3
for(j in c(0.3, 0.6, 0.9)){
  uptake.pris <- j
  source(paste0(pris_path, "/scripts/vacc-scenarios-uptake.R"))
  results_df <- results_df %>% group_by(population, scenario_run) %>% summarize(total=sum(total)) %>% 
    mutate(uptake.staff="Staff uptake = 0.3", uptake.pris=paste0("Prisoner uptake = ", j))
  cases_uptake <- cases_uptake %>% rbind(results_df)
}
  
uptake.staff <- 0.6  
for(j in c(0.3, 0.6, 0.9)){
  uptake.pris <- j
  source(paste0(pris_path, "/scripts/vacc-scenarios-uptake.R"))
  results_df <- results_df %>% group_by(population, scenario_run) %>% summarize(total=sum(total)) %>% 
    mutate(uptake.staff="Staff uptake = 0.6", uptake.pris=paste0("Prisoner uptake = ", j))
  cases_uptake <- cases_uptake %>% rbind(results_df)
}

uptake.staff <- 0.9
for(j in c(0.3, 0.6, 0.9)){
  uptake.pris <- j
  source(paste0(pris_path, "/scripts/vacc-scenarios-uptake.R"))
  results_df <- results_df %>% group_by(population, scenario_run) %>% summarize(total=sum(total)) %>% 
    mutate(uptake.staff="Staff uptake = 0.9", uptake.pris=paste0("Prisoner uptake = ", j))
  cases_uptake <- cases_uptake %>% rbind(results_df)
}

cases_uptake <- labelforplots(cases_uptake)

cases_fig <- cases_uptake %>% dplyr::filter(scenario == "(1) no vaccination") %>% group_by(scenario_run, population, scenario, scenario_nr, uptake.staff, uptake.pris) %>%
    summarize(total_base=total) %>% ungroup %>% dplyr::select(population, total_base, uptake.staff, uptake.pris) %>%
    full_join(cases_uptake, by=c("population", "uptake.staff", "uptake.pris")) %>%
    mutate(value = round((1-total/total_base)*100, 1))
  
  cases_fig <- cases_fig %>% mutate(case.avert=total_base - total)
  cases_fig <- popnchange(cases_fig)
  
  
  cases_stack <- cases_fig %>% group_by(scenario, uptake.staff, uptake.pris) %>% mutate(mean.bar=sum(total), total.value=round((100-sum(total)/sum(total_base)*100),1), overall.base=sum(total_base)) %>% ungroup()
  
  case.stack.plot <- ggplot(cases_stack, aes(x = scenario, y = total, fill=scenario, alpha=population), show_guide=F) +
    geom_bar(stat = "identity", show.legend = T) +
    scale_alpha_discrete(range=c(0.4,1)) +
    scale_x_discrete(labels=unique(cases_stack$scenario_nr)) +
    coord_cartesian(ylim = c(0, 860)) +
    labs(x="Vaccination scenario", y="Clinical cases over one year", fill="Scenario", alpha="Population") +
    theme_bw() +
    scale_color_OkabeIto() +  
    theme(text=element_text(size=10), legend.text=element_text(size=12)) +
    scale_fill_viridis(discrete=TRUE) +
    geom_text(aes(x = scenario, y = mean.bar,
                  label = ifelse(total.value==0,
                                 paste0(format(total.value, nsmall = 1),"%"),
                                 paste0("-",format(total.value, nsmall = 1),"%")) ),
              angle  = 0, vjust  = -0.25, hjust  = 0.5, size   = 2.5, lineheight = 0.8, colour = "black") +
    geom_hline(aes(yintercept=overall.base), linetype="dashed", alpha=0.5, colour = "red") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    guides(
      alpha = guide_legend(
        title = "Population",
        override.aes = aes(label = ""))) +
    facet_grid(uptake.staff~uptake.pris)
  
  ggsave(paste0(save_path,"fig-uptake.png"), case.stack.plot, width=260, height=210, units="mm")
  
  

