
if(exists("run")){rm(run)}
n <- 1 # number of runs = 1 (currently deterministic)

new1 <- params$pop[[1]]$size*staff_to*uptake.staff
new2 <- params$pop[[2]]$size*staff_to*uptake.staff
new3 <- params$pop[[3]]$size*b_rate*uptake.pris

dose_gap <- 7*12

lambda <- prev*contacts*u*undetect

# 0.003 and 0.0206 = high and low rolling 14-day average proportion of PCR tests coming back positive -> if a prevalence, does it matter if daily or fortnightly?

import.total.NPF <- sum(lambda*params$pop[[1]]$size)
import.total.PF <- sum(lambda*params$pop[[2]]$size)

seed_freq.NPF <- 1/import.total.NPF
seed_freq.PF <- 1/import.total.PF

# Imported infections

params$pop[[1]]$seed_times <- round(c(seq(from=delay, to=delay+round(365.25+timehorizon), by=seed_freq.NPF)), 0)
params$pop[[2]]$seed_times <- round(c(seq(from=delay, to=delay+round(365.25*timehorizon), by=seed_freq.PF)), 0)
params$pop[[3]]$seed_times <- delay
for(i in nr_pops){params$pop[[i]]$dist_seed_ages <- c(rep(0,3), rep(1,13))}

# scenarios
  # (1) no vaccination
  for(i in 1:nr_pops){
    params$pop[[i]]$v = rep(0, 16)
    params$pop[[i]]$v12 = rep(0,16) # (re)set to no vaccines being administered
  }
  params$schedule <- list()         # no scheduled changes to parameters
  run1 <- cm_simulate(params, 1)
  results_run1 = run1$dynamics[compartment %in% c("death_o", "cases", "icu_p", "nonicu_p", "onedose_i", "twodose_i"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run1")
  
  plot_check1 <- run1$dynamics[compartment %in% c("S", "R", "V", "V2"),] %>% group_by(run, t, compartment) %>% summarise(daily_sum = sum(value)) %>%
    rbind(run1$dynamics[compartment %in% c("Is", "Ip", "Ia"),] %>% group_by(run,t) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="I")) %>%
    rbind(run1$dynamics[compartment %in% c("E", "Ea"),] %>% group_by(run,t) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="E")) %>% ungroup()
  plot_check1$compartment <- factor(plot_check1$compartment, levels=c("S", "E", "I", "R", "V", "V2"))
  
  total_pop1 <- run1$dynamics[compartment %in% c("S", "E", "R", "V", "V2", "Ip", "Is", "Ia", "Ea"),] %>% group_by(run, t) %>% summarise(total_popn = sum(value)) %>% ungroup() %>% ggplot(aes(x=t,y=total_popn)) + geom_line()
  all_compartments1 <- run1$dynamics[compartment %in% c("S", "E", "R", "V", "V2", "Ip", "Is", "Ia", "Ea"),] %>% group_by(run, t, compartment) %>% summarise(daily_sum = sum(value)) %>% ungroup() %>% ggplot(aes(x=t,y=daily_sum,colour=compartment,linetype=compartment, group=compartment)) + geom_line(show.legend = T)
  rm(run1)
  
  # (2) Just non-prisoner facing (non-operational) staff
  
  ## Number of doses needed for this strategy:
  scen2 <- sum(params$pop[[1]]$size)*uptake.staff
  ## Number of days needed to administer this number of doses
  done2 <- scen2/n_vacc_daily
  ## Accounting for new people arriving whilst vaccination campaign is happening:
  while((sum(params$pop[[1]]$size)*staff_to*done2 + sum(params$pop[[1]]$size))*uptake.staff > n_vacc_daily*done2){
    done2 <- done2+0.01
  }
  if(done2>12*7){cat(red("WARNING: Administration of first and second dose overlap - need to account for in vaccination rates"))}
  
  vacc_vals2 <- c(n_vacc_daily * params$pop[[1]]$size/sum(params$pop[[1]]$size))
  
  for(i in 1:nr_pops){
    params$pop[[i]]$v = rep(0, 16)
    params$pop[[i]]$v12 = rep(0,16) # (re)set to no vaccines being administered
  }
  params$schedule <- list() 
  params$schedule = list(
    list(
      parameter = "v",
      pops = 0,
      mode = "assign",
      values = list(vacc_vals2, new1),
      times = c(immune1,immune1+done2)),
    list(
      parameter = "v12",
      pops = 0,
      mode = "assign",
      values = list(vacc_vals2, new1),
      times = c(immune2, immune2+done2))
    )
  # Both first and second doses only administered for a short period required to vaccinate up to assumed uptake.
  
  run2 = cm_simulate(params, 1)
 #  results_run2 = run2$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run2")
  results_run2 = run2$dynamics[compartment %in% c("death_o", "cases", "icu_p", "nonicu_p", "onedose_i", "twodose_i"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run2")
  
  plot_check2 <- run2$dynamics[compartment %in% c("S", "R", "V", "V2"),] %>% group_by(run, t, compartment) %>% summarise(daily_sum = sum(value)) %>%
    rbind(run2$dynamics[compartment %in% c("Is", "Ip", "Ia"),] %>% group_by(run,t) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="I")) %>%
    rbind(run2$dynamics[compartment %in% c("E", "Ea"),] %>% group_by(run,t) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="E")) %>% ungroup()
  plot_check2$compartment <- factor(plot_check2$compartment, levels=c("S", "E", "I", "R", "V", "V2"))

  rm(run2)
  # (3) Just prisoner facing staff
  
  ## Number of doses needed for this strategy:
  ## Number of doses needed for this strategy:
  scen3 <- sum(params$pop[[2]]$size)*uptake.staff
  ## Number of days needed to administer this number of doses
  done3 <- scen3/n_vacc_daily
  ## Accounting for new people arriving whilst vaccination campaign is happening:
  while((sum(params$pop[[2]]$size)*staff_to*done3 + sum(params$pop[[2]]$size))*uptake.staff > n_vacc_daily*done3){
    done3 <- done3+0.01
  }
  if(done3>12*7){cat(red("Administration of first and second dose overlap - need to account for in vaccination rates"))}
  
  
  vacc_vals3 <-   c(n_vacc_daily * params$pop[[2]]$size/sum(params$pop[[2]]$size))
  # vacc_vals3 <- list
  for(i in 1:nr_pops){
  params$pop[[i]]$v <- rep(0,16)
  params$pop[[i]]$v12 <- rep(0,16)
  }
  params$schedule <- list() 
  params$schedule = list(
    list(
      parameter = "v",
      pops = 1,
      mode = "assign",
      values = list(vacc_vals3, new2),
      times = c(immune1, immune1+done3)),
    list(
      parameter = "v12",
      pops = 1,
      mode = "assign",
      values = list(vacc_vals3, new2),
      times = c(immune2, immune2+done3)))
  
  run3 = cm_simulate(params, 1)
 #  results_run3 = run3$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run3")
  results_run3 = run3$dynamics[compartment %in% c("death_o", "cases", "icu_p", "nonicu_p", "onedose_i", "twodose_i"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run3")
  
  plot_check3 <- run3$dynamics[compartment %in% c("S", "R", "V", "V2"),] %>% group_by(run, t, compartment) %>% summarise(daily_sum = sum(value)) %>%
    rbind(run3$dynamics[compartment %in% c("Is", "Ip", "Ia"),] %>% group_by(run,t) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="I")) %>%
    rbind(run3$dynamics[compartment %in% c("E", "Ea"),] %>% group_by(run,t) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="E")) %>% ungroup()
  plot_check3$compartment <- factor(plot_check3$compartment, levels=c("S", "E", "I", "R", "V", "V2"))

  total_pop3 <- run3$dynamics[compartment %in% c("S", "E", "R", "V", "V2", "Ip", "Is", "Ia", "Ea"),] %>% group_by(run, t) %>% summarise(total_popn = sum(value)) %>% ungroup() %>% ggplot(aes(x=t,y=total_popn)) + geom_line()
  all_compartments3 <- run3$dynamics[compartment %in% c("S", "E", "R", "V", "V2", "Ip", "Is", "Ia", "Ea"),] %>% group_by(run, t, compartment) %>% summarise(daily_sum = sum(value)) %>% ungroup() %>% ggplot(aes(x=t,y=daily_sum,colour=compartment,linetype=compartment, group=compartment)) + geom_line(show.legend = T)

  
  # Plot by population
  bypopn_3 <- run3$dynamics[compartment %in% c("S", "R", "V", "V2"),] %>% group_by(run, t, compartment, population) %>% summarise(daily_sum = sum(value)) %>%
    rbind(run3$dynamics[compartment %in% c("Is", "Ip", "Ia"),] %>% group_by(run,t, population) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="I")) %>%
    rbind(run3$dynamics[compartment %in% c("E", "Ea"),] %>% group_by(run,t, population) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="E")) %>% ungroup()
  bypopnplot_3 <- bypopn_3 %>% filter(compartment %in% c("S", "V", "V2")) %>% ggplot(aes(x=t,y=daily_sum,colour=compartment,linetype=population)) + geom_line(show.legend = T, lwd=0.5, alpha=0.7) + theme_bw()
  
  rm(run3)

  # (4) All prison staff
  
  ## Number of doses needed for this strategy:
  scen4 <- sum(params$pop[[1]]$size+params$pop[[2]]$size)*uptake.staff
  ## Number of days needed to administer this number of doses
  done4 <- scen4/n_vacc_daily
  ## Accounting for the new people arriving during vaccination campaign:
  while((sum(params$pop[[1]]$size)*staff_to*done4 + sum(params$pop[[2]]$size)*staff_to*done4 + sum(params$pop[[1]]$size + params$pop[[2]]$size))*uptake.staff > n_vacc_daily*done4){
    done4 <- done4+0.01
  }
  
  if(done4>12*7){cat(red("WARNING: Administration of first and second dose overlap - need to account for in vaccination rates"))}
  prop.pop1 <- sum(params$pop[[1]]$size)/sum(params$pop[[1]]$size+params$pop[[2]]$size)
  prop.pop2 <- sum(params$pop[[2]]$size)/sum(params$pop[[1]]$size+params$pop[[2]]$size)
  
  vacc_vals4.pop1 <- c(n_vacc_daily*prop.pop1*params$pop[[1]]$size/sum(params$pop[[1]]$size))
  vacc_vals4.pop2 <- c(n_vacc_daily*prop.pop2*params$pop[[2]]$size/sum(params$pop[[2]]$size))
  
  # vacc_vals4 <- list(c(rep(0,3), rep(1, 13)))
  params$schedule <- list() 
  for(i in 1:nr_pops){
    params$pop[[i]]$v = rep(0, 16)
    params$pop[[i]]$v12 = rep(0,16) # (re)set to no vaccines being administered
  }
  params$schedule = list(
    list(
      parameter = "v",
      pops = c(0), # populations are from 0 to N-1
      mode = "assign",
      values = list(vacc_vals4.pop1, new1),
      times = c(immune1, immune1+done4)),
    list(
      parameter = "v12",
      pops = c(0),
      mode = "assign",
      values = list(vacc_vals4.pop1, new1),
      times = c(immune2, immune2+done4)),
    list(
      parameter = "v",
      pops = c(1), # populations are from 0 to N-1
      mode = "assign",
      values = list(vacc_vals4.pop2, new2),
      times = c(immune1, immune2+done4)),
    list(
      parameter = "v12",
      pops = c(1),
      mode = "assign",
      values = list(vacc_vals4.pop2, new2),
      times = c(immune2, immune2+done4)))

  
  run4 = cm_simulate(params, 1)
  # results_run4 = run4$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run4")
  results_run4 = run4$dynamics[compartment %in% c("death_o", "cases", "icu_p", "nonicu_p", "onedose_i", "twodose_i"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run4")
  
  plot_check4 <- run4$dynamics[compartment %in% c("S", "R", "V", "V2"),] %>% group_by(run, t, compartment) %>% summarise(daily_sum = sum(value)) %>%
    rbind(run4$dynamics[compartment %in% c("Is", "Ip", "Ia"),] %>% group_by(run,t) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="I")) %>%
    rbind(run4$dynamics[compartment %in% c("E", "Ea"),] %>% group_by(run,t) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="E")) %>% ungroup()
  plot_check4$compartment <- factor(plot_check4$compartment, levels=c("S", "E", "I", "R", "V", "V2"))
  
  # Plot by population
  bypopn_4 <- run4$dynamics[compartment %in% c("S", "R", "V", "V2"),] %>% group_by(run, t, compartment, population) %>% summarise(daily_sum = sum(value)) %>%
    rbind(run4$dynamics[compartment %in% c("Is", "Ip", "Ia"),] %>% group_by(run,t, population) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="I")) %>%
    rbind(run4$dynamics[compartment %in% c("E", "Ea"),] %>% group_by(run,t, population) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="E")) %>% ungroup()
  bypopnplot_4 <- bypopn_4 %>% filter(compartment %in% c("S", "V", "V2")) %>% ggplot(aes(x=t,y=daily_sum,colour=compartment,linetype=population)) + geom_line(show.legend = T, lwd=0.5, alpha=0.7) + theme_bw()
  rm(run4)

  gc()
  
  # (5) prison population
  
  # (5) prison population
  
  scen5 <- sum(params$pop[[3]]$size)*uptake.pris
  ## Number of days needed to administer this number of doses
  done5 <- scen5/n_vacc_daily
  
  while((sum(params$pop[[3]]$size)*b_rate*done5 + sum(params$pop[[3]]$size))*uptake.pris > n_vacc_daily*done5){
    done5 <- done5+0.01
  }
  if(done5>12*7){cat(red("WARNING: Administration of first and second dose overlap - need to account for in vaccination rates"))}
  
  
  
  ## Number of doses needed for this strategy:
  params$schedule <- list() 
  for(i in 1:nr_pops){
    params$pop[[i]]$v = rep(0, 16)
    params$pop[[i]]$v12 = rep(0,16) # (re)set to no vaccines being administered
  }
  vacc_vals5 <- c(n_vacc_daily * params$pop[[3]]$size/sum(params$pop[[3]]$size))
  params$schedule = list(
    list(
      parameter = "v",
      pops = 2,
      mode = "assign",
      values = list(vacc_vals5, new3),
      times = c(immune1, immune1+done5)),
    list(
      parameter = "v12",
      pops = 2,
      mode = "assign",
      values = list(vacc_vals5, new3),
      times = c(immune2, immune2+done5)))
  
  run5 = cm_simulate(params, 1)
  # results_run5 = run5$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run5")
  results_run5 = run5$dynamics[compartment %in% c("death_o", "cases", "icu_p", "nonicu_p", "onedose_i", "twodose_i"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run5")
  
  plot_check5 <- run5$dynamics[compartment %in% c("S", "R", "V", "V2"),] %>% group_by(run, t, compartment) %>% summarise(daily_sum = sum(value)) %>%
    rbind(run5$dynamics[compartment %in% c("Is", "Ip", "Ia"),] %>% group_by(run,t) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="I")) %>%
    rbind(run5$dynamics[compartment %in% c("E", "Ea"),] %>% group_by(run,t) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="E")) %>% ungroup()
  plot_check5$compartment <- factor(plot_check5$compartment, levels=c("S", "E", "I", "R", "V", "V2"))
  
  gc()
  
  #6 Just vulnerable (over 50) population
  total.popn <-  sum(params$pop[[1]]$size[11:16]+params$pop[[2]]$size[11:16]+params$pop[[3]]$size[11:16])
  ## Number of doses needed for this strategy:
  scen6 <- sum(params$pop[[1]]$size[11:16]+params$pop[[2]]$size[11:16])*uptake.staff + sum(params$pop[[3]]$size[11:16])*uptake.pris
  ## Number of days needed to administer this number of doses
  done6 <- scen6/n_vacc_daily
  ## Accounting for new arrivals over period of vaccine programme
  while(((sum(params$pop[[1]]$size[11:16])*staff_to*done6 + sum(params$pop[[2]]$size[11:16])*staff_to*done6)*uptake.staff + sum(params$pop[[3]]$size[11:16])*b_rate*done6*uptake.pris + scen6) > n_vacc_daily*done6){
    done6 <- done6+0.01
  }
  

  if(done6>12*7){cat(red("WARNING: Administration of first and second dose overlap - need to account for in vaccination rates"))}
  prop.pop1 <- sum(params$pop[[1]]$size[11:16])/sum(params$pop[[1]]$size[11:16]+params$pop[[2]]$size[11:16]+params$pop[[3]]$size[11:16])
  prop.pop2 <- sum(params$pop[[2]]$size[11:16])/sum(params$pop[[1]]$size[11:16]+params$pop[[2]]$size[11:16]+params$pop[[3]]$size[11:16])
  prop.pop3 <- sum(params$pop[[3]]$size[11:16])/sum(params$pop[[1]]$size[11:16]+params$pop[[2]]$size[11:16]+params$pop[[3]]$size[11:16])
  
  
  vacc_vals6.pop1 <-c(rep(0,10), n_vacc_daily*prop.pop1*params$pop[[1]]$size[11:16]/sum(params$pop[[1]]$size[11:16]))
  vacc_vals6.pop2 <-c(rep(0,10), n_vacc_daily*prop.pop2*params$pop[[2]]$size[11:16]/sum(params$pop[[2]]$size[11:16]))
  vacc_vals6.pop3 <-c(rep(0,10), n_vacc_daily*prop.pop3*params$pop[[3]]$size[11:16]/sum(params$pop[[3]]$size[11:16]))

  params$schedule <- list() 
  for(i in 1:nr_pops){
    params$pop[[i]]$v = rep(0, 16)
    params$pop[[i]]$v12 = rep(0,16) # (re)set to no vaccines being administered
  }
  params$schedule = list(
    list(
      parameter = "v",
      pops = c(0),
      mode = "assign",
      values = list(vacc_vals6.pop1, c(rep(0,10), new1[11:16])),
      times = c(immune1, immune1+done6)),
    list(
      parameter = "v12",
      pops = c(0),
      mode = "assign",
      values = list(vacc_vals6.pop1, c(rep(0,10), new1[11:16])),
      times = c(immune2, immune2+done6)),
    list(
      parameter = "v",
      pops = c(1),
      mode = "assign",
      values = list(vacc_vals6.pop2, c(rep(0,10), new2[11:16])),
      times = c(immune1, immune1+done6)),
    list(
      parameter = "v12",
      pops = c(1),
      mode = "assign",
      values = list(vacc_vals6.pop2, c(rep(0,10), new2[11:16])),
      times = c(immune2, immune2+done6)),
    list(
      parameter = "v",
      pops = c(2),
      mode = "assign",
      values = list(vacc_vals6.pop3, c(rep(0,10), new3[11:16])),
      times = c(immune1, immune1+done6)),
    list(
      parameter = "v12",
      pops = c(2),
      mode = "assign",
      values = list(vacc_vals6.pop3, c(rep(0,10), new3[11:16])),
      times = c(immune2, immune2+done6)))

  run6 = cm_simulate(params, 1)
  # results_run6 = run6$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run6")
  results_run6 = run6$dynamics[compartment %in% c("death_o", "cases", "icu_p", "nonicu_p", "onedose_i", "twodose_i"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run6")

  plot_check6 <- run6$dynamics[compartment %in% c("S", "R", "V", "V2"),] %>% group_by(run, t, compartment) %>% summarise(daily_sum = sum(value)) %>%
    rbind(run6$dynamics[compartment %in% c("Is", "Ip", "Ia"),] %>% group_by(run,t) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="I")) %>%
    rbind(run6$dynamics[compartment %in% c("E", "Ea"),] %>% group_by(run,t) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="E")) %>% ungroup()
  plot_check6$compartment <- factor(plot_check6$compartment, levels=c("S", "E", "I", "R", "V", "V2"))

#  rm(run6)

  gc()
  
  #7 All
  
  ## Number of doses needed for this strategy:
  total.popn <- sum(params$pop[[1]]$size+params$pop[[2]]$size+params$pop[[3]]$size)
  scen7 <- sum(params$pop[[1]]$size+params$pop[[2]]$size)*uptake.staff + sum(params$pop[[3]]$size)*uptake.pris
  done7 <- scen7/n_vacc_daily
  while(((sum(params$pop[[1]]$size)*staff_to*done7 + sum(params$pop[[2]]$size)*staff_to*done7)*uptake.staff + sum(params$pop[[3]]$size)*b_rate*done7*uptake.pris + scen7) > n_vacc_daily*done7){
    done7 <- done7+0.01
  }


  # vacc_vals7 <- list(c(rep(0,3), rep(1, 13)))
  params$schedule <- list() 
  for(i in 1:nr_pops){
    params$pop[[i]]$v = rep(0, 16)
    params$pop[[i]]$v12 = rep(0,16) # (re)set to no vaccines being administered
  }
  params$schedule = list(
    list(
      parameter = "v",
      pops = c(0),
      mode = "assign",
      values = list(vacc_vals7.pop1, new1),
      times = c(immune1, immune1+done7)),
    list(
      parameter = "v12",
      pops = c(0),
      mode = "assign",
      values = list(vacc_vals7.pop1),
      times = c(immune2)),
    list(
        parameter = "v",
        pops = c(1),
        mode = "assign",
        values = list(vacc_vals7.pop2, new2),
        times = c(immune1, immune1+done7)),
      list(
        parameter = "v12",
        pops = c(1),
        mode = "assign",
        values = list(vacc_vals7.pop2),
        times = c(immune2)),
      list(
        parameter = "v",
        pops = c(2),
        mode = "assign",
        values = list(vacc_vals7.pop3, new3),
        times = c(immune1, immune1+done7)),
      list(
        parameter = "v12",
        pops = c(2),
        mode = "assign",
        values = list(vacc_vals7.pop3),
        times = c(immune2)))
  
  gc()
  run7 = cm_simulate(params, 1)
  # results_run7 = run7$dynamics[c(compartment == "cases_i"), .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run7")
  results_run7 = run7$dynamics[compartment %in% c("death_o", "cases", "icu_p", "nonicu_p", "onedose_i", "twodose_i"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run7")
  
  plot_check7 <- run7$dynamics[compartment %in% c("S", "R", "V", "V2"),] %>% group_by(run, t, compartment) %>% summarise(daily_sum = sum(value)) %>%
    rbind(run7$dynamics[compartment %in% c("Is", "Ip", "Ia"),] %>% group_by(run,t) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="I")) %>%
    rbind(run7$dynamics[compartment %in% c("E", "Ea"),] %>% group_by(run,t) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="E")) %>% ungroup()
  plot_check7$compartment <- factor(plot_check7$compartment, levels=c("S", "E", "I", "R", "V", "V2"))
  
  # Plot by population
  bypopn_7 <- run7$dynamics[compartment %in% c("S", "R", "V", "V2"),] %>% group_by(run, t, compartment, population) %>% summarise(daily_sum = sum(value)) %>%
    rbind(run7$dynamics[compartment %in% c("Is", "Ip", "Ia"),] %>% group_by(run,t, population) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="I")) %>%
    rbind(run7$dynamics[compartment %in% c("E", "Ea"),] %>% group_by(run,t, population) %>% summarise(daily_sum=sum(value)) %>% mutate(compartment="E")) %>% ungroup()
  bypopnplot_7 <- bypopn_7 %>% filter(compartment %in% c("S", "V", "V2")) %>% ggplot(aes(x=t,y=daily_sum,colour=compartment,linetype=population)) + geom_line(show.legend = T, lwd=0.5, alpha=0.7) + theme_bw()
  
  rm(run7)
  gc()
  
  # get all scenarios
  
  results_df = do.call("rbind", list(results_run1, results_run2, results_run3, results_run4, results_run5, results_run6, results_run7)) %>% dplyr::filter(compartment=="cases")
  
  # same for qalys 
  
  results_qalys = do.call("rbind", list(results_run1, results_run2, results_run3, results_run4, results_run5, results_run6, results_run7)) %>%
    group_by(scenario_run, run, group, population, compartment) %>% summarise(total=sum(total))
  
  results_qalys <- results_qalys %>%
    group_by(scenario_run, population, group, compartment) %>% 
    summarise(total=sum(total)) %>% ungroup()
  
  results_dose = do.call("rbind", list(results_run1, results_run2, results_run3, results_run4, results_run5, results_run6, results_run7)) %>% dplyr::filter(compartment %in% c("onedose_i", "twodose_i"))
  totaldose <- results_dose %>% group_by(scenario_run, group, population) %>% summarise(total=sum(total))
  
  results_qalys <- rbind(results_qalys, totaldose %>% mutate(compartment="aefi.minor"))
  results_qalys <- rbind(results_qalys, totaldose %>% mutate(compartment="aefi.fatal"))
  
  results_qalys <- left_join(results_qalys, qalycalc, by=c("compartment", "group", "population")) %>% filter(compartment!=c("onedose_i", "twodose_i"))
  results_qalys <- results_qalys %>% mutate(qaly.loss=total*qaly.value) %>% 
    group_by(scenario_run, population) %>% 
    summarise(qaly.loss=sum(qaly.loss)) %>% 
    ungroup()
  

  results_qalys <- rbind(results_qalys, results_qalys %>% group_by(scenario_run) %>% summarise(qaly.loss=sum(qaly.loss), population="(A-C) all prisoners and staff"))
  
  ## same for deaths
  
  results_deaths = do.call("rbind", list(results_run1, results_run2, results_run3, results_run4, results_run5, results_run6, results_run7)) %>% dplyr::filter(compartment=="death_o")
  
  
  ## Doses
  
  results_dose = do.call("rbind", list(results_run1, results_run2, results_run3, results_run4, results_run5, results_run6, results_run7)) %>% dplyr::filter(compartment %in% c("onedose_i", "twodose_i"))
  totaldose <- results_dose %>% group_by(scenario_run, group) %>% summarise(total=sum(total))
  
  rm(results_run1, results_run2, results_run3, results_run4, results_run5, results_run6, results_run7)
  