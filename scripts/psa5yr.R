###############################################
##### Running scenarios using psa values ######
###############################################

psafive <- function(n_psa, psa){
  n_psa <- n_psa
  total_psa <- data.frame()
  time_psa <- data.frame()
  psa <- psa
  
  tic()
  params_list <- list()
  for(k in 1:n_psa){
    psa.values <- psa[k,]
    for(i in 1:nr_pops){
      params$pop[[i]]$wn <- rep(psa.values$imm_nat, 16)
      params$pop[[i]]$wv <- rep(psa.values$imm_vac, 16)
      params$pop[[i]]$wv2 <- rep(psa.values$imm_vac, 16)
      params$pop[[i]]$ei_v <- rep(psa.values$eff_inf1,16) # vaccine efficacy against infection, first dose
      params$pop[[i]]$ei_v2 <- rep(psa.values$eff_inf2,16) # vaccine efficacy against infection, second dose
      params$pop[[i]]$ed_vi <- rep(psa.values$eff_dis1,16) # vaccine efficacy against disease given infection, first dose (i.e. no additional protection)
      params$pop[[i]]$ed_vi2 <- rep(psa.values$eff_dis2,16) # vaccine efficacY against disease given infection, second dose (0.78-0.67)/(1-0.67)
    }
    
    params$pop[[1]]$B <- psa.values$staff_to * age_staff
    params$pop[[2]]$B <- psa.values$staff_to * age_staff
    params$pop[[1]]$D <- rep(psa.values$staff_to,16)
    params$pop[[2]]$D <- rep(psa.values$staff_to,16)
    
    params$pop[[3]]$B <- psa.values$b_rate * age_prisoners # daily
    params$pop[[3]]$D <- rep(psa.values$b_rate, 16)
    
    
    ## QALYs ##
    
    qalycalc$qaly.value[qalycalc$compartment=="cases"] <- psa.values$qaly.sym
    qalycalc$qaly.value[qalycalc$compartment=="icu_p"] <- psa.values$qaly.icu
    qalycalc$qaly.value[qalycalc$compartment=="nonicu_p"] <- psa.values$qaly.nonicu
    
    
    ## R0 - scaling ##
    
    for(i in 1:nr_pops){
      # age-dependent ratio of infection:cases (based on Davies et al, Nature paper)
      params$pop[[i]]$y <- c(0.2904047, 0.2904047, 0.2070468, 0.2070468, 0.2676134,
                             0.2676134, 0.3284704, 0.3284704, 0.3979398, 0.3979398,
                             0.4863355, 0.4863355, 0.6306967, 0.6306967, 0.6906705, 0.6906705)
      
      # susceptibility (based on Davies et al, Nature paper)
      params$pop[[i]]$u <- c(0, 0, 0, 0.3815349, 0.7859512,
                             0.7859512, 0.8585759, 0.8585759, 0.7981468, 0.7981468,
                             0.8166960, 0.8166960, 0.8784811, 0.8784811, 0.7383189, 0.7383189)
      
      
      # scale u (susceptibility) to achieve desired R0
      current_R0 = cm_calc_R0(params, i); # calculate R0 in population i of params
      params$pop[[i]]$u = params$pop[[i]]$u * psa.values$target_R0 / current_R0
    } 
    
    ## Community prevalence and contacts ##
    
    # Prevalence * number of contacts * susceptibility * proportion undetected 
    
    psa.lambda <- psa.values$prev*contacts*u*(1-psa.values$LFD.uptake)*(1-psa.values$LFD.sens)*(1-(psa.values$vac.uptake*psa.values$eff_inf2))
    
    # 0.003 and 0.0206 = high and low rolling 14-day average proportion of PCR tests coming back positive -> if a prevalence, does it matter if daily or fortnightly?
    
    psa.seed_freq.NPF <- 1/sum(psa.lambda*params$pop[[1]]$size)
    psa.seed_freq.PF <- 1/sum(psa.lambda*params$pop[[2]]$size)
    
    # Imported infections
    
    params$pop[[1]]$seed_times <- round(c(seq(from=365.25, to=365.25+round(365.25*timehorizon), by=psa.seed_freq.NPF)), 0)
    params$pop[[2]]$seed_times <- round(c(seq(from=365.25, to=365.25+round(365.25*timehorizon), by=psa.seed_freq.PF)), 0)
    params$pop[[3]]$seed_times <- 365.25
    for(i in nr_pops){params$pop[[i]]$dist_seed_ages <- c(rep(0,3), rep(1,13))}
    
    uptake <- psa.values$vac.uptake
    
    #### VACCINATION RATE - NEW PRISONERS/STAFF
    new1 <- n_vacc_daily/(sum(params$pop[[1]]$size)*psa.values$staff_to*uptake)
    new2 <- n_vacc_daily/(sum(params$pop[[2]]$size)*psa.values$staff_to*uptake)
    new3 <- n_vacc_daily/(sum(params$pop[[3]]$size)*psa.values$b_rate*uptake)
    
    prop.pop1 <- sum(params$pop[[1]]$size)/sum(params$pop[[1]]$size+params$pop[[2]]$size+params$pop[[3]]$size)
    prop.pop2 <- sum(params$pop[[2]]$size)/sum(params$pop[[1]]$size+params$pop[[2]]$size+params$pop[[3]]$size)
    prop.pop3 <- sum(params$pop[[3]]$size)/sum(params$pop[[1]]$size+params$pop[[2]]$size+params$pop[[3]]$size)
    
    vacc_vals7.pop1 <- c((n_vacc_daily-5)*prop.pop1*params$pop[[1]]$size/sum(params$pop[[1]]$size))
    vacc_vals7.pop2 <- c((n_vacc_daily-5)*prop.pop2*params$pop[[2]]$size/sum(params$pop[[2]]$size))
    vacc_vals7.pop3 <- c((n_vacc_daily-5)*prop.pop3*params$pop[[3]]$size/sum(params$pop[[3]]$size)) + 5*params$pop[[3]]$size/sum(params$pop[[3]]$size)
    
    new1.7 <- sum(vacc_vals7.pop1)/(sum(params$pop[[1]]$size)*psa.values$staff_to*uptake)
    new2.7 <- sum(vacc_vals7.pop2)/(sum(params$pop[[2]]$size)*psa.values$staff_to*uptake)
    new3.7 <- sum(vacc_vals7.pop3)/(sum(params$pop[[3]]$size)*psa.values$b_rate*uptake)
    # Vacc short #
    
    if(exists("run")){rm(run)}
    n <- 1 # number of runs = 1 (currently deterministic)
    
    dose_gap <- 7*12
    
    # scenarios
    # (1) no vaccination
    for(i in 1:nr_pops){
      params$pop[[i]]$v = rep(0, 16)
      params$pop[[i]]$v12 = rep(0,16) # (re)set to no vaccines being administered
    }
    params$schedule <- list()         # no scheduled changes to parameters
    run1 <- cm_simulate(params, 1)
    results_run1 = run1$dynamics[compartment %in% c("cases"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run1")
    rm(run1)
    
    # (2) Just non-prisoner facing (non-operational) staff
    
    ## Number of doses needed for this strategy:
    scen2 <- sum(params$pop[[1]]$size)*uptake
    ## Number of days needed to administer this number of doses
    done2 <- scen2/n_vacc_daily
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
        values = list(vacc_vals2, vacc_vals2/new1),
        times = c(immune1,immune1+done2)),
      list(
        parameter = "v12",
        pops = 0,
        mode = "assign",
        values = list(vacc_vals2),
        times = c(immune2))
    )
    # Both first and second doses only administered for a short period required to vaccinate up to assumed uptake.
    
    run2 = cm_simulate(params, 1)
    #  results_run2 = run2$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run2")
    results_run2 = run2$dynamics[compartment %in% c("cases"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run2")
    rm(run2)
    # (3) Just prisoner facing staff
    
    ## Number of doses needed for this strategy:
    scen3 <- sum(params$pop[[2]]$size)*uptake
    ## Number of days needed to administer this number of doses
    done3 <- scen3/n_vacc_daily
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
        values = list(vacc_vals3, vacc_vals3/new2),
        times = c(immune1, immune1+done3)),
      list(
        parameter = "v12",
        pops = 1,
        mode = "assign",
        values = list(vacc_vals3),
        times = c(immune2)))
    
    run3 = cm_simulate(params, 1)
    #  results_run3 = run3$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run3")
    results_run3 = run3$dynamics[compartment %in% c("cases"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run3")
    rm(run3)
    gc()
    
    # (4) All prison staff
    
    ## Number of doses needed for this strategy:
    scen4 <- sum(params$pop[[1]]$size+params$pop[[2]]$size)*uptake
    ## Number of days needed to administer this number of doses
    done4 <- scen4/n_vacc_daily
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
        values = list(vacc_vals4.pop1, vacc_vals4.pop1/new1),
        times = c(immune1, immune1+done4)),
      list(
        parameter = "v12",
        pops = c(0),
        mode = "assign",
        values = list(vacc_vals4.pop1),
        times = c(immune2)),
      list(
        parameter = "v",
        pops = c(1), # populations are from 0 to N-1
        mode = "assign",
        values = list(vacc_vals4.pop2, vacc_vals4.pop2/new2),
        times = c(immune1, immune1+done4)),
      list(
        parameter = "v12",
        pops = c(1),
        mode = "assign",
        values = list(vacc_vals4.pop2),
        times = c(immune2)))
    
    
    run4 = cm_simulate(params, 1)
    # results_run4 = run4$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run4")
    results_run4 = run4$dynamics[compartment %in% c("cases"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run4")
    rm(run4)
    gc()
    
    # (5) prison population
    
    scen5 <- sum(params$pop[[3]]$size)*uptake
    ## Number of days needed to administer this number of doses
    done5 <- scen5/n_vacc_daily
    
    while((sum(params$pop[[3]]$size)*b_rate*done5 + sum(params$pop[[3]]$size))*0.9 > 20*done5){
      done5 <- done5+1
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
        values = list(vacc_vals5, vacc_vals5/new3),
        times = c(immune1, immune1+done5)),
      list(
        parameter = "v12",
        pops = 2,
        mode = "assign",
        values = list(vacc_vals5),
        times = c(immune2)))
    
    run5 = cm_simulate(params, 1)
    # results_run5 = run5$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run5")
    results_run5 = run5$dynamics[compartment %in% c("cases"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run5")
    gc()
    rm(run5)
    
    #6 Just vulnerable (over 50) population
    
    total.popn <-  sum(params$pop[[1]]$size[11:16]+params$pop[[2]]$size[11:16]+params$pop[[3]]$size[11:16])
    ## Number of doses needed for this strategy:
    scen6 <-total.popn*uptake
    ## Number of days needed to administer this number of doses
    done6 <- scen6/n_vacc_daily
    ## Accounting for new arrivals over period of vaccine programme
    while((sum(params$pop[[1]]$size[11:16])*staff_to*done6 + sum(params$pop[[2]]$size[11:16])*staff_to*done6 + sum(params$pop[[3]]$size[11:16])*b_rate*done6 + total.popn)*0.9 > 20*done6){
      done6 <- done6+1
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
        values = list(vacc_vals6.pop1, vacc_vals6.pop1/new1),
        times = c(immune1, immune1+done6)),
      list(
        parameter = "v12",
        pops = c(0),
        mode = "assign",
        values = list(vacc_vals6.pop1),
        times = c(immune2)),
      list(
        parameter = "v",
        pops = c(1),
        mode = "assign",
        values = list(vacc_vals6.pop2, vacc_vals6.pop2/new2),
        times = c(immune1, immune1+done6)),
      list(
        parameter = "v12",
        pops = c(1),
        mode = "assign",
        values = list(vacc_vals6.pop2),
        times = c(immune2)),
      list(
        parameter = "v",
        pops = c(2),
        mode = "assign",
        values = list(vacc_vals6.pop3, vacc_vals6.pop3/new3),
        times = c(immune1, immune1+done6)),
      list(
        parameter = "v12",
        pops = c(2),
        mode = "assign",
        values = list(vacc_vals6.pop3),
        times = c(immune2)))
    
    run6 = cm_simulate(params, 1)
    # results_run6 = run6$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run6")
    results_run6 = run6$dynamics[compartment %in% c("cases"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run6")
    rm(run6)
    gc()
    
    #7 All
    
    ## Number of doses needed for this strategy:
    total.popn <- sum(params$pop[[1]]$size+params$pop[[2]]$size+params$pop[[3]]$size)
    scen7 <- total.popn*uptake
    done7 <- scen7/n_vacc_daily
    while((prisoner_pop*b_rate*done7 + 70*staff_to*done7 + 315*staff_to*done7 + total.popn)*0.9 > 20*done7){
      done7 <- done7+1
    }
    
    ## Number of days needed to administer this number of doses
    
    if(done7>12*7){cat(red("WARNING: Administration of first and second dose overlap - need to account for in vaccination rates"))}
    
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
        values = list(vacc_vals7.pop1, vacc_vals7.pop1/new1.7),
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
        values = list(vacc_vals7.pop2, vacc_vals7.pop2/new2.7),
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
        values = list(vacc_vals7.pop3, vacc_vals7.pop3/new3.7),
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
    results_run7 = run7$dynamics[compartment %in% c("cases"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run7")
    rm(run7)
    gc()
    
    # get all scenarios
    
    results_df = do.call("rbind", list(results_run1, results_run2, results_run3, results_run4, results_run5, results_run6, results_run7)) %>% dplyr::filter(compartment=="cases")

    
    ### END OF VACC SHORT ###
    
    results_df <- results_df %>% dplyr::filter(compartment=="cases") %>% 
      mutate(run=k)
    case.total_psa <- results_df %>% group_by(scenario_run, run) %>% summarise(total=sum(total))
    case.time_psa <- results_df %>% group_by(scenario_run, run, t) %>% summarise(total=sum(total))
    total_psa <- rbind(total_psa, case.total_psa)
    time_psa <- rbind(time_psa, case.time_psa)
    rm(results_df)
    rm(case.total_psa)
    rm(case.time_psa)
  }
  toc()

  return(list(total_psa, time_psa))}
