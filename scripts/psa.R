##############
#### PSA #####
##############

n_psa <- 100
n_parms <- 16
seeds <- c(123, 456, 789, 101, 999)
tic()
save_path <- "C:/Users/CiaraMcCarthy/Documents/prisons-vacc-strategies/outputs/"
for(z in seeds){
  source("~/Documents/prisons-vacc-strategies/scripts/set-up-upgrading.R")
  psa <- createpsa(z, n_psa)
  psa_temp <- psafunc(n_psa, psa)
  write.csv(psa_temp[[1]], paste0(save_path,"case-total", z, ".csv"))
  write.csv(psa_temp[[2]], paste0(save_path,"case-time", z, ".csv"))
  write.csv(psa_temp[[3]], paste0(save_path,"qaly-total", z, ".csv"))
  write.csv(psa_temp[[4]], paste0(save_path,"death-total", z, ".csv"))
  write.csv(psa_temp[[5]], paste0(save_path,"cases-prcc", z, ".csv"))
  write.csv(psa_temp[[6]], paste0(save_path,"qalys-prcc", z, ".csv"))
  write.csv(psa_temp[[7]], paste0(save_path,"dose-total", z, ".csv"))
  rm(list=setdiff(ls(), c("n_psa", "n_parms", "seeds", "psa")))
}
toc()

source("~/Documents/prisons-vacc-strategies/scripts/figure1.R")
source("~/Documents/prisons-vacc-strategies/scripts/figure2.R")
source("~/Documents/prisons-vacc-strategies/scripts/figure3.R")
source("~/Documents/prisons-vacc-strategies/scripts/figure4.R")
source("~/Documents/prisons-vacc-strategies/scripts/figure5.R")
source("~/Documents/prisons-vacc-strategies/scripts/figure6.R")
source("~/Documents/prisons-vacc-strategies/scripts/figure7.R")

source("~/Documents/prisons-vacc-strategies/scripts/supfig1.R")
source("~/Documents/prisons-vacc-strategies/scripts/supfig2.R")
source("~/Documents/prisons-vacc-strategies/scripts/supfig3.R")
source("~/Documents/prisons-vacc-strategies/scripts/supfig4.R")
source("~/Documents/prisons-vacc-strategies/scripts/supfig5-7.R")

