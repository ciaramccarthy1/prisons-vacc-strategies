## Modified psa function that includes seed times over three years and changes reinfection parameter at t= 2 years

psa_newvar <- function(n_psa, psa){
  results_psa <- data.frame()
  qaly_psa <- data.frame()
  deaths_psa <- data.frame()
  dose_psa <- data.frame()
  
  tic()
  params_list <- list()
  for(k in 1:n_psa){
    psa.values <- psa[k,]
    
  #  # Waning vaccine immunity 
  #  psa.values$imm_vac <- log((psa.values$eff_inf2-psa.values$vac_decay)/psa.values$eff_inf2)/-140
  #  # Waning natural immunity
  #  psa.values$imm_nat <- log(1-psa.values$nat_decay)/-365.25
  #  
  #  # Efficacy against disease
  #  psa.values$eff_dis1 <- VEdis(psa.values$eff_tot1, psa.values$eff_inf1)
  #  psa.values$eff_dis2 <- VEdis(psa.values$eff_tot2, psa.values$eff_inf2)
    
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
    
    ## Changing waning of natural immunity to duration of e.g. 6 months:
    for(i in 1:nr_pops){
    params$pop[[i]]$wn <- rep((1/(365.25/2)), 16)}
    
    
    ## QALYs ##
    
    qalycalc$qaly.value[qalycalc$compartment=="cases"] <- psa.values$qaly.sym
    qalycalc$qaly.value[qalycalc$compartment=="to_icu_i"] <- psa.values$qaly.icu
    qalycalc$qaly.value[qalycalc$compartment=="to_nonicu_i"] <- psa.values$qaly.nonicu
    
    
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
    
    # CHANGED - imported infections
    
    params$pop[[1]]$seed_times <- round(c(seq(from=365.25, to=365.25+round(365.25*timehorizon), by=psa.seed_freq.NPF)), 0)
    params$pop[[2]]$seed_times <- round(c(seq(from=365.25, to=365.25+round(365.25*timehorizon), by=psa.seed_freq.PF)), 0)
    params$pop[[3]]$seed_times <- 365.25
    for(i in nr_pops){params$pop[[i]]$dist_seed_ages <- c(rep(0,3), rep(1,13))}
    
    uptake <- psa.values$vac.uptake
    
    #### VACCINATION RATE - NEW PRISONERS/STAFF
    new1 <- params$pop[[1]]$size*psa.values$staff_to*uptake
    new2 <- params$pop[[2]]$size*psa.values$staff_to*uptake
    new3 <- params$pop[[3]]$size*psa.values$b_rate*uptake
    
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
    params$schedule <- list(
     list(
       parameter = "pi_r",
       pops = c(0,1,2),
       mode = "assign",
       values = list(c(rep(0.8,16))),
       times = c(365.25)
     )
    )
    run1 <- cm_simulate(params, 1)
    results_run1 = run1$dynamics[compartment %in% c("death_o", "cases", "to_icu_i", "to_nonicu_i", "onedose_i", "twodose_i"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run1")
    rm(run1)
    
    # (2) Just non-prisoner facing (non-operational) staff
    
    ## Number of doses needed for this strategy:
    scen2 <- sum(params$pop[[1]]$size)*uptake
    ## Number of days needed to administer this number of doses
    done2 <- scen2/n_vacc_daily
    ## Accounting for new people arriving whilst vaccination campaign is happening:
    while((sum(params$pop[[1]]$size)*staff_to*done2 + sum(params$pop[[1]]$size))*uptake > n_vacc_daily*done2){
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
        values = list(vacc_vals2),
        times = c(immune2)),
      list(
        parameter = "pi_r",
        pops = c(0,1,2),
        mode = "assign",
        values = list(c(rep(0.8,16))),
        times = c(365.25)
      )
    )
    # Both first and second doses only administered for a short period required to vaccinate up to assumed uptake.
    
    run2 = cm_simulate(params, 1)
    #  results_run2 = run2$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run2")
    results_run2 = run2$dynamics[compartment %in% c("death_o", "cases", "to_icu_i", "to_nonicu_i", "onedose_i", "twodose_i"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run2")
    rm(run2)
    # (3) Just prisoner facing staff
    
    ## Number of doses needed for this strategy:
    scen3 <- sum(params$pop[[2]]$size)*uptake
    ## Number of days needed to administer this number of doses
    done3 <- scen3/n_vacc_daily
    ## Accounting for new people arriving whilst vaccination campaign is happening:
    while((sum(params$pop[[2]]$size)*staff_to*done3 + sum(params$pop[[2]]$size))*uptake > n_vacc_daily*done3){
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
        values = list(vacc_vals3),
        times = c(immune2)),
     list(
       parameter = "pi_r",
       pops = c(0,1,2),
       mode = "assign",
       values = list(c(rep(0.8,16))),
       times = c(365.25)))
    
    run3 = cm_simulate(params, 1)
    #  results_run3 = run3$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run3")
    results_run3 = run3$dynamics[compartment %in% c("death_o", "cases", "to_icu_i", "to_nonicu_i", "onedose_i", "twodose_i"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run3")
    rm(run3)
    gc()
    
    # (4) All prison staff
    
    ## Number of doses needed for this strategy:
    scen4 <- sum(params$pop[[1]]$size+params$pop[[2]]$size)*uptake
    ## Number of days needed to administer this number of doses
    done4 <- scen4/n_vacc_daily
    ## Accounting for the new people arriving during vaccination campaign:
    while((sum(params$pop[[1]]$size)*staff_to*done4 + sum(params$pop[[2]]$size)*staff_to*done4 + sum(params$pop[[1]]$size + params$pop[[2]]$size))*uptake > n_vacc_daily*done4){
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
        values = list(vacc_vals4.pop1),
        times = c(immune2)),
      list(
        parameter = "v",
        pops = c(1), # populations are from 0 to N-1
        mode = "assign",
        values = list(vacc_vals4.pop2, new2),
        times = c(immune1, immune1+done4)),
      list(
        parameter = "v12",
        pops = c(1),
        mode = "assign",
        values = list(vacc_vals4.pop2),
        times = c(immune2)),
      list(
        parameter = "pi_r",
        pops = c(0,1,2),
        mode = "assign",
        values = list(c(rep(0.8,16))),
        times = c(365.25))
      )
    
    
    run4 = cm_simulate(params, 1)
    # results_run4 = run4$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run4")
    results_run4 = run4$dynamics[compartment %in% c("death_o", "cases", "to_icu_i", "to_nonicu_i", "onedose_i", "twodose_i"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run4")
    rm(run4)
    gc()
    
    # (5) prison population
    
    scen5 <- sum(params$pop[[3]]$size)*uptake
    ## Number of days needed to administer this number of doses
    done5 <- scen5/n_vacc_daily
    
    while((sum(params$pop[[3]]$size)*b_rate*done5 + sum(params$pop[[3]]$size))*uptake > n_vacc_daily*done5){
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
        values = list(vacc_vals5),
        times = c(immune2)),
      list(
        parameter = "pi_r",
        pops = c(0,1,2),
        mode = "assign",
        values = list(c(rep(0.8,16))),
        times = c(365.25)))
    
    run5 = cm_simulate(params, 1)
    # results_run5 = run5$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run5")
    results_run5 = run5$dynamics[compartment %in% c("death_o", "cases", "to_icu_i", "to_nonicu_i", "onedose_i", "twodose_i"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run5")
    gc()
    rm(run5)
    
    #6 Just vulnerable (over 50) population
    
    total.popn <-  sum(params$pop[[1]]$size[11:16]+params$pop[[2]]$size[11:16]+params$pop[[3]]$size[11:16])
    ## Number of doses needed for this strategy:
    scen6 <-total.popn*uptake
    ## Number of days needed to administer this number of doses
    done6 <- scen6/n_vacc_daily
    ## Accounting for new arrivals over period of vaccine programme
    while((sum(params$pop[[1]]$size[11:16])*staff_to*done6 + sum(params$pop[[2]]$size[11:16])*staff_to*done6 + sum(params$pop[[3]]$size[11:16])*b_rate*done6 + total.popn)*uptake > n_vacc_daily*done6){
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
        values = list(vacc_vals6.pop1, new1),
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
        values = list(vacc_vals6.pop2, new2),
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
        values = list(vacc_vals6.pop3, new3),
        times = c(immune1, immune1+done6)),
      list(
        parameter = "v12",
        pops = c(2),
        mode = "assign",
        values = list(vacc_vals6.pop3),
        times = c(immune2)),
      list(
        parameter = "pi_r",
        pops = c(0,1,2),
        mode = "assign",
        values = list(c(rep(0.8,16))),
        times = c(365.25)))
    
    run6 = cm_simulate(params, 1)
    # results_run6 = run6$dynamics[compartment == "cases_i", .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run6")
    results_run6 = run6$dynamics[compartment %in% c("death_o", "cases", "to_icu_i", "to_nonicu_i", "onedose_i", "twodose_i"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run6")
    rm(run6)
    gc()
    
    #7 All
    
    ## Number of doses needed for this strategy:
    total.popn <- sum(params$pop[[1]]$size+params$pop[[2]]$size+params$pop[[3]]$size)
    scen7 <- total.popn*uptake
    done7 <- scen7/n_vacc_daily
    while((prisoner_pop*b_rate*done7 + sum(params$pop[[1]]$size)*staff_to*done7 + sum(params$pop[[2]]$size)*staff_to*done7 + total.popn)*uptake > n_vacc_daily*done7){
      done7 <- done7+0.01
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
        times = c(immune2)),
      list(
        parameter = "pi_r",
        pops = c(0,1,2),
        mode = "assign",
        values = list(c(rep(0.8,16))),
        times = c(365.25)))
    
    gc()
    run7 = cm_simulate(params, 1)
    # results_run7 = run7$dynamics[c(compartment == "cases_i"), .(total = sum(value)), by = .(run, population, t)] %>% mutate(scenario_run="run7")
    results_run7 = run7$dynamics[compartment %in% c("death_o", "cases", "to_icu_i", "to_nonicu_i", "onedose_i", "twodose_i"), .(total = sum(value)), by = .(run, population, group, compartment, t)] %>% mutate(scenario_run="run7")
    rm(run7)
    gc()
    
    # get all scenarios
    
    results_df = do.call("rbind", list(results_run1, results_run2, results_run3, results_run4, results_run5, results_run6, results_run7)) %>% dplyr::filter(compartment=="cases")
    
    # same for qalys 
    
    results_qalys = do.call("rbind", list(results_run1, results_run2, results_run3, results_run4, results_run5, results_run6, results_run7)) %>%
      group_by(scenario_run, run, group, population, compartment) %>% summarise(total=sum(total))
    
    results_qalys <- results_qalys %>%
      group_by(scenario_run, population, group, compartment) %>% 
      summarise(total=mean(total)) %>% ungroup()
    
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
    
    ### END OF VACC SHORT ###
    
    results_df <- results_df %>% dplyr::filter(compartment=="cases") %>% 
      mutate(run=k)
    results_psa <- rbind(results_psa, results_df)
    results_qalys <- results_qalys %>% mutate(run=k)
    qaly_psa <- rbind(results_qalys, qaly_psa)
    results_deaths <- results_deaths %>% mutate(run=k)
    deaths_psa <- rbind(deaths_psa, results_deaths)
    results_dose <- results_dose %>% mutate(run=k)
    dose_psa <- rbind(dose_psa, results_dose)
    rm(results_df)
    rm(results_qalys)
    rm(results_deaths)
    rm(results_dose)
  }
  toc()
  psa <- psa %>% mutate(run=1:n_psa)
  cases_prcc <- results_psa %>% group_by(scenario_run, run) %>% summarise(total=sum(total)) %>% left_join(psa, by=c("run"))
  case.total <- results_psa %>% group_by(scenario_run, run) %>% summarise(total=sum(total))
  case.time <- results_psa %>% group_by(scenario_run, t, run) %>% summarise(total=sum(total))
  rm(results_psa)
  qaly.total <- qaly_psa %>% group_by(scenario_run, run) %>% summarise(total=sum(qaly.loss))
  qaly_prcc <- qaly_psa %>% group_by(scenario_run, run) %>% summarise(total=sum(qaly.loss)) %>% left_join(psa, by=c("run"))
  rm(qaly_psa)
  dose.total <- dose_psa %>% filter(compartment=="twodose_i") %>% group_by(scenario_run, run) %>% summarise(vacc.count=sum(total))
  rm(dose_psa)
  death.total <- deaths_psa %>% group_by(scenario_run, run) %>% summarise(total=sum(total))
  rm(deaths_psa)
  return(list(case.total, case.time, qaly.total, death.total, cases_prcc, qaly_prcc, dose.total))}
