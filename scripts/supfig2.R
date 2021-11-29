##################################
####### SUPPLEMENTARY FIG 2 ######
##################################

# PRCC - QALYs

rm(list = ls())

n_psa <- 100
n_parms <- 16
seeds <- c(123, 456, 789, 101, 999)
source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")

qaly_prcc <- data.frame()
for(i in 1:5){
  qaly.prcc.part <- read.csv(paste0(save_path,"qalys-prcc", seeds[i], ".csv")) %>% select(2:(n_parms+4))
  qaly.prcc.part <- qaly.prcc.part %>% mutate(run=run + (i-1)*n_psa)
  qaly_prcc <- rbind(qaly_prcc, qaly.prcc.part)
}

prcc_list <- lapply(unique(qaly_prcc$scenario_run), function(x) qaly_prcc[qaly_prcc$scenario_run == x,])
for(i in 1:7){
  df <- prcc_list[[i]]
  df <- df %>% ungroup %>% select(imm_vac:vac.uptake, total)
  
  assign(paste0("run", i),  df)
}

for(i in 1:7){
  x <- get(paste0("run", i))
  prcc <- epi.prcc(x, sided.test=2, conf.level=0.95)
  prcc <- prcc %>% mutate(scenario_run=paste0("run", i), parameter=colnames(qaly_prcc)[4:(n_parms+3)])
  assign(paste0("run", i, ".prcc"), prcc)
}
qaly_prcc <- do.call("rbind", list(run1.prcc, run2.prcc, run3.prcc, run4.prcc, run5.prcc, run6.prcc, run7.prcc))
qaly_prcc <- labelforplots(qaly_prcc)

parm.labs <- c("Duration of vaccine immunity", "Duration of natural immunity", "Staff turnover rate",
               "Prisoner turnover rate", "QALY loss - symptomatic",
               "QALY loss - non-ICU hospitalisation",
               "QALY loss - ICU",
               "VE against infection - 1st dose",
               "VE against infection - 2nd dose",
               "VE against disease - 1st dose",
               "VE against disease - 2nd dose",
               "R0", 
               "Community prevalence",
               "LFD uptake",
               "LFD sensitivity",
               "Vaccine uptake")


qaly_prcc$parameter <- ordered(qaly_prcc$parameter, levels=c("eff_inf1", "eff_inf2", "eff_dis1", "eff_dis2", 
                                                             "imm_vac", "imm_nat", "staff_to", "b_rate", 
                                                             "vac.uptake", "LFD.uptake", "LFD.sens", "prev",
                                                             "target_R0", "qaly.sym", "qaly.nonicu", "qaly.icu"))
param.qaly <- unique(qaly_prcc$parameter)
names(parm.labs) <- param.qaly

fig_prcc_qaly <- ggplot(qaly_prcc, aes(x=scenario_nr, y=est, fill=scenario)) + geom_col() + 
  scale_fill_viridis(discrete=TRUE) +
  theme_bw() + theme(text=element_text(size=7.5)) +
  labs(x="Scenario", y="PRCC", fill="Scenario") + 
  geom_errorbar(mapping=aes(ymin=lower, ymax=upper), width=0.2) + 
  facet_wrap(~parameter, labeller=labeller(parameter=parm.labs))

ggsave(paste0(save_path,"supfig2new.png"), fig_prcc_qaly, width=190, height=170, units="mm")
