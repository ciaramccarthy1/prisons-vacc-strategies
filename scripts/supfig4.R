#################################
###### SUPPLEMENTARY FIG 4 ######
#################################

# Waning vs. no waning

rm(list = ls())
source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")
source("~/Documents/prisons-vacc-strategies/scripts/psa5yr.R")
timehorizon <- 5
params$time1 <- "2026-02-01"
n_psa <- 100
seeds <- c(123, 456, 789, 101, 999)

for(z in seeds){
  source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")
  timehorizon <- 5
  params$time1 <- "2026-02-01"
  set.seed(z)
  n_parms <- 14
  
  sampleLSH <-randomLHS(n=n_psa, k=n_parms)
  
  results_psa <- data.frame()
  psa.wane <- data.frame(
    imm_nat=c(rep(0, n_psa)),
    imm_vac=c(rep(0, n_psa)),
    staff_to=qbeta(sampleLSH[,1], shape1=staff.parms$alpha, shape2=staff.parms$beta), # /365.25 to get daily rate
    b_rate=qbeta(sampleLSH[,2], shape1=pris.parms$alpha, shape2=pris.parms$beta),
    qaly.sym=qbeta(sampleLSH[,3], shape1=sym.parms$alpha, shape2=sym.parms$beta),
    qaly.nonicu=qbeta(sampleLSH[,4], shape1=nonicu.parms$alpha, shape2=nonicu.parms$beta),
    qaly.icu=qbeta(sampleLSH[,5], shape1=icu.parms$alpha, shape2=icu.parms$beta),
    eff_inf1=qbeta(sampleLSH[,6], shape1=infeff1.parms$alpha, shape2=infeff1.parms$beta),
    eff_inf2=qbeta(sampleLSH[,7], shape1=infeff2.parms$alpha, shape2=infeff2.parms$beta),
    eff_dis1=qbeta(sampleLSH[,8], shape1=diseff1.parms$alpha, shape2=diseff1.parms$beta),
    eff_dis2=qbeta(sampleLSH[,9], shape1=diseff2.parms$alpha, shape2=diseff2.parms$beta),
    target_R0=qlnorm(sampleLSH[,10], meanlog=R0.parms$meanlog, sdlog=R0.parms$sdlog),
    prev=qunif(sampleLSH[,11], min=0.0003, max=0.0206),
    LFD.uptake=qbeta(sampleLSH[,12], shape1=lfd.uptake.parms$alpha, shape2=lfd.uptake.parms$beta),
    LFD.sens=qbeta(sampleLSH[,13], shape1=lfd.parms$alpha, shape2=lfd.parms$beta),
    vac.uptake=qbeta(sampleLSH[,14], shape1=vac.parms$alpha, shape2=vac.parms$beta)
  )
  
  # No waning
  params_list <- list()
  results_nowane <- psafive(n_psa, psa.wane)
  results_nowane[[1]] <- results_nowane[[1]] %>% mutate(wane="No waning")
  write.csv(results_nowane[[1]], paste0(save_path,"nowane", z, ".csv"))
  
  psa.wane <- psa.wane %>% mutate(imm_nat=c(rep(1/(365.25*0.5), n_psa)), imm_vac=c(rep(1/(365.25*0.5), n_psa)))
  results_wane <- psafive(n_psa, psa.wane)
  results_wane[[1]] <- results_wane[[1]] %>% mutate(wane="Waning")
  write.csv(results_wane[[1]], paste0(save_path,"wane", z, ".csv"))
  rm(list=setdiff(ls(), c("n_psa", "n_parms", "seeds", "psa.wane", "psafive")))
}


# results_psa <- results_psa %>% group_by(scenario_run, run) %>% summarise(total=sum(total)) %>% mutate(wane="No waning")
# write.csv(results_psa, paste0("nowane", z, ".csv"))

source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")

cases_fig11 <- data.frame()
for(i in 1:5){
  case.nowane <- read.csv(paste0(save_path,"nowane", seeds[i], ".csv")) %>% select(scenario_run, run, total, wane)
  case.nowane <- case.nowane %>% mutate(run=run + (i-1)*n_psa)
  case.wane <- read.csv(paste0(save_path,"wane", seeds[i], ".csv")) %>% select(scenario_run, run, total, wane)
  case.wane <- case.wane %>% mutate(run=run + (i-1)*n_psa)
  cases_fig11 <- rbind(cases_fig11, case.nowane) %>% rbind(case.wane)
}

cases_fig11 <- labelforplots(cases_fig11)

# cases_fig11 <- cases_fig11 %>% group_by(scenario_run, scenario, scenario_nr, wane) %>% summarise(median=median(total), lower=quantile(total, probs=0.25), upper=quantile(total, probs=0.75), min=min(total), max=max(total))

fig11a <- ggplot(cases_fig11) +
  geom_boxplot(aes(x=wane, y=total, fill=scenario), outlier.shape=NA, alpha=0.7) + facet_wrap(~scenario) + theme_bw() +
  scale_fill_viridis(discrete=TRUE, name="Scenario") + labs(x="Waning", y="Sum of clinical cases over five years")

fig11b <- ggplot(cases_fig11) +
  geom_boxplot(aes(x=scenario_nr, y=total, fill=scenario), outlier.shape=NA, alpha=0.7) + facet_wrap(~wane, ncol=1) + theme_bw() +
  scale_fill_viridis(discrete=TRUE, name="Scenario") + labs(x="Scenario", y="Sum of clinical cases over five years")

fig11c <- ggplot(cases_fig11) +
  geom_boxplot(aes(x=scenario_nr, y=total, fill=scenario), outlier.shape=NA, alpha=0.7) + facet_wrap(~wane) + theme_bw() +
  scale_fill_viridis(discrete=TRUE, name="Scenario") + labs(x="Scenario", y="Sum of clinical cases over five years") + theme(text=element_text(size=7))

ggsave("fig11a.png", fig11a, width=300, height=210, units="mm")
ggsave("fig11b.png", fig11b, width=180, height=120, units="mm")
ggsave(paste0(save_path,"supfig4new.png"),fig11c, width=190, height=55, units="mm")

# Back to base case
target_R0 <- 5
timehorizon <- 1
params$time1 <- "2022-02-01"
source(paste0(pris_path, "/scripts/sensitivity.R"))
params$pop[[1]]$seed_times <- round(c(seq(from=365.25, to=365.25*(timehorizon+1), by=seed_freq.NPF)), 0)
params$pop[[2]]$seed_times <- round(c(seq(from=365.25, to=365.25*(timehorizon+1), by=seed_freq.PF)), 0)
params$pop[[3]]$seed_times <- 365.25
immune1 <- ((params$time0)+28)
immune2 <- ((params$time0+(7*12)+14))
time0 <- 0

