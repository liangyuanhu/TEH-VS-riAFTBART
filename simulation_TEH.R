#============================================================================#
# Simulations                                                                #                  
#============================================================================#
library(survival)
library(tidyverse)
library(riAFTBART)
library(coxme)
library(lme4)
library(pammtools)
library(mgcv)
library(timereg)
library(WeightIt)
library(RISCA)
# Source the functions to calculate median survival times
source("code/data_gen_TEH.R")
source("code/riAFT-BART-TEH.R")
source("code/DR-riAH-TEH.R")
source("code/IPW-riCox-TEH.R")
source("code/PEAMM-TEH.R")

# Simulation scenario HS(i): PH-------------------------------------------------------- 
for (i in 1:250){
  set.seed(1111+i)
  sim_data <- data_gen(nk=10, K=200, PH = TRUE, overlap = "moderate to strong", HS = "i")
  res_ATEs_riAFT_BART <- ri_AFTBART(M.burnin = 2000, M.keep = 3500, status = sim_data$delta, y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl, verbose = FALSE)
  res_ATEs_IPW_ri_Cox <- IPW_ri_Cox(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  res_PEAMM <- PEAMM(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  res_DR_riAH <- DR_riAH(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  save(res_ATEs_riAFT_BART, res_ATEs_IPW_ri_Cox, res_PEAMM, res_DR_riAH, file = paste0("result_scenario_HSi_PH_",i,".Rdata"))
}

# Simulation scenario HS(ii): PH-------------------------------------------------------- 
for (i in 1:250){
  set.seed(1111+i)
  sim_data <- data_gen(nk=10, K=200, PH = TRUE, overlap = "moderate to strong", HS = "ii")
  res_ATEs_riAFT_BART <- ri_AFTBART(M.burnin = 2000, M.keep = 3500, status = sim_data$delta, y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl, verbose = FALSE)
  res_ATEs_IPW_ri_Cox <- IPW_ri_Cox(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  res_PEAMM <- PEAMM(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  res_DR_riAH <- DR_riAH(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  save(res_ATEs_riAFT_BART, res_ATEs_IPW_ri_Cox, res_PEAMM, res_DR_riAH, file = paste0("result_scenario_HSii_PH_",i,".Rdata"))
}
# Simulation scenario HS(iii): PH-------------------------------------------------------- 
for (i in 1:250){
  set.seed(1111+i)
  sim_data <- data_gen(nk=10, K=200, PH = TRUE, overlap = "moderate to strong", HS = "iii")
  res_ATEs_riAFT_BART <- ri_AFTBART(M.burnin = 2000, M.keep = 3500, status = sim_data$delta, y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl, verbose = FALSE)
  res_ATEs_IPW_ri_Cox <- IPW_ri_Cox(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  res_PEAMM <- PEAMM(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  res_DR_riAH <- DR_riAH(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  save(res_ATEs_riAFT_BART, res_ATEs_IPW_ri_Cox, res_PEAMM, res_DR_riAH, file = paste0("result_scenario_HSiii_PH_",i,".Rdata"))
}
# Simulation scenario HS(i): nPH-------------------------------------------------------- 
for (i in 1:250){
  set.seed(1111+i)
  sim_data <- data_gen(nk=10, K=200, PH = FALSE, overlap = "moderate to strong", HS = "i")
  res_ATEs_riAFT_BART <- ri_AFTBART(M.burnin = 2000, M.keep = 3500, status = sim_data$delta, y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl, verbose = FALSE)
  res_ATEs_IPW_ri_Cox <- IPW_ri_Cox(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  res_PEAMM <- PEAMM(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  res_DR_riAH <- DR_riAH(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  save(res_ATEs_riAFT_BART, res_ATEs_IPW_ri_Cox, res_PEAMM, res_DR_riAH, file = paste0("result_scenario_HSi_nPH_",i,".Rdata"))
}

# Simulation scenario HS(ii): nPH-------------------------------------------------------- 
for (i in 1:250){
  set.seed(1111+i)
  sim_data <- data_gen(nk=10, K=200, PH = FALSE, overlap = "moderate to strong", HS = "ii")
  res_ATEs_riAFT_BART <- ri_AFTBART(M.burnin = 2000, M.keep = 3500, status = sim_data$delta, y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl, verbose = FALSE)
  res_ATEs_IPW_ri_Cox <- IPW_ri_Cox(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  res_PEAMM <- PEAMM(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  res_DR_riAH <- DR_riAH(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  save(res_ATEs_riAFT_BART, res_ATEs_IPW_ri_Cox, res_PEAMM, res_DR_riAH, file = paste0("result_scenario_HSii_nPH_",i,".Rdata"))
}
# Simulation scenario HS(iii): nPH-------------------------------------------------------- 
for (i in 1:250){
  set.seed(1111+i)
  sim_data <- data_gen(nk=10, K=200, PH = FALSE, overlap = "moderate to strong", HS = "iii")
  res_ATEs_riAFT_BART <- ri_AFTBART(M.burnin = 2000, M.keep = 3500, status = sim_data$delta, y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl, verbose = FALSE)
  res_ATEs_IPW_ri_Cox <- IPW_ri_Cox(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  res_PEAMM <- PEAMM(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  res_DR_riAH <- DR_riAH(y = sim_data$Tobs_C, x = sim_data$X, trt = sim_data$A, cluster.id = sim_data$cl,status = sim_data$delta)
  save(res_ATEs_riAFT_BART, res_ATEs_IPW_ri_Cox, res_PEAMM, res_DR_riAH, file = paste0("result_scenario_HSiii_nPH_",i,".Rdata"))
}