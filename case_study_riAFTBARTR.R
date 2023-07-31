library(tidyverse)
library(riAFTBART)
library(mice)
data_post_ICU_varselect <- read_csv("code/riAFTBART VS/covid_post_icu.csv") %>% 
  mutate_if(is.character, as.factor) %>% 
  select(-treatment) %>% 
  as.data.frame()
# Fit riAFTBART using selected variables---------------------------------------------------------------
library(riAFTBART)
data_post_icu <- read_csv("code/riAFTBART VS/covid_post_icu.csv") %>% 
  mutate_if(is.character, as.factor) %>% 
  select(-treatment) %>% 
  as.data.frame()
data_post_icu %>% count(treatment)
X_covariates <- data_post_icu %>% 
  select(d_dimer, Creatinine, Oxygen_ordinal, age,wbc, Race, sofa_score, ldh) %>%
  mutate_if(is.character, as.factor) %>% 
  as.data.frame()
for (i in 1:n_boot_impute){
  # Bootstarp
  id_boot<-sample(c(1:nrow(X_covariates)), nrow(X_covariates), replace=TRUE)
  X_covariates_b<-X_covariates[c(id_boot),]
  mice_b<-mice(X_covariates_b, m=1, maxit = maxit_impute)
  mice_b_list<-mice_b 
  mice_b_dat<-complete(mice_b, include=FALSE)
  mice_b_dat_list[[i]] <- mice_b_dat
}
impute_dat_list <- mice_b_dat_list

for (i in 1:n_boot_impute){
  riAFTBART_fit_trt_1 <- riAFTBART_fit(M.burnin = 1000,
                                       M.keep = 3500,
                                       status = data_post_icu$death,
                                       y.train = data_post_icu$time_icu_to_outcome,
                                       x.train = X_covariates,
                                       trt.train = data_post_icu$treatment,
                                       x.test = X_covariates,
                                       trt.test = 1,
                                       cluster.id = data_post_icu$LocationName,
                                       verbose = T)
  save(riAFTBART_fit_trt_1, file = paste0("riAFTBART_fit_trt_1_",n_boot_impute,"._Rdata"))
  riAFTBART_fit_trt_2 <- riAFTBART_fit(M.burnin = 1000,
                                       M.keep = 3500,
                                       status = data_post_icu$death,
                                       y.train = data_post_icu$time_icu_to_outcome,
                                       x.train = X_covariates,
                                       trt.train = data_post_icu$treatment,
                                       x.test = X_covariates,
                                       trt.test = 2,
                                       cluster.id = data_post_icu$LocationName,
                                       verbose = T)
  save(riAFTBART_fit_trt_2, file = paste0("riAFTBART_fit_trt_2_",n_boot_impute,"._Rdata"))
  riAFTBART_fit_trt_3 <- riAFTBART_fit(M.burnin = 1000,
                                       M.keep = 3500,
                                       status = data_post_icu$death,
                                       y.train = data_post_icu$time_icu_to_outcome,
                                       x.train = X_covariates,
                                       trt.train = data_post_icu$treatment,
                                       x.test = X_covariates,
                                       trt.test = 3,
                                       cluster.id = data_post_icu$LocationName,
                                       verbose = T)
  save(riAFTBART_fit_trt_3, file = paste0("riAFTBART_fit_trt_3_",n_boot_impute,"._Rdata"))
}


