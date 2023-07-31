library(tidyverse)
library(riAFTBART)
library(mice)
data_post_ICU_varselect <- read_csv("code/riAFTBART VS/covid_post_icu.csv") %>% 
  mutate_if(is.character, as.factor) %>% 
  select(-treatment) %>% 
  as.data.frame()
# Fit the fit using Intree ---------------------------------------------------------------------
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

library(inTrees)
library(randomForest)



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
  
  ite23 <- apply(riAFTBART_fit_trt_2$tree.pred,1,mean)-apply(riAFTBART_fit_trt_3$tree.pred,1,mean)
  X_covariates_ite23 <- X_covariates %>% 
    as_tibble() %>% 
    mutate(ite = apply(riAFTBART_fit_trt_2$tree.pred,1,mean)-apply(riAFTBART_fit_trt_3$tree.pred,1,mean))
  ite12 <- apply(riAFTBART_fit_trt_1$tree.pred,1,mean)-apply(riAFTBART_fit_trt_2$tree.pred,1,mean)
  X_covariates_ite12 <- X_covariates %>% 
    as_tibble() %>% 
    mutate(ite = apply(riAFTBART_fit_trt_1$tree.pred,1,mean)-apply(riAFTBART_fit_trt_2$tree.pred,1,mean))
  ite13 <- apply(riAFTBART_fit_trt_1$tree.pred,1,mean)-apply(riAFTBART_fit_trt_3$tree.pred,1,mean)
  X_covariates_ite13 <- X_covariates %>% 
    as_tibble() %>% 
    mutate(ite = apply(riAFTBART_fit_trt_1$tree.pred,3,mean)-apply(riAFTBART_fit_trt_2$tree.pred,1,mean))
  rf_ite12 <- randomForest(x = X_covariates, y = ite12) 
  tree_list_ite12 <- RF2List(rf_ite12)
  rule_exec_ite12 <- unique(extractRules(tree_list_ite12,X_covariates,digits=3))
  rule_metric_ite12 <- getRuleMetric(rule_exec_ite12,X_covariates,ite12)
  rule_metric_ite12 <- pruneRule(rule_metric_ite12,X_covariates,ite12,typeDecay = 1)
  rule_select_ite12 <- selectRuleRRF(rule_metric_ite12, X_covariates,ite12)
  save(rule_select_ite12, file = paste0("Intree_ite12_",n_boot_impute,"._Rdata"))
  rf_ite23 <- randomForest(x = X_covariates, y = ite23) 
  tree_list_ite23 <- RF2List(rf_ite23)
  rule_exec_ite23 <- unique(extractRules(tree_list_ite23,X_covariates,digits=3))
  rule_metric_ite23 <- getRuleMetric(rule_exec_ite23,X_covariates,ite23)
  rule_metric_ite23 <- pruneRule(rule_metric_ite23,X_covariates,ite23,typeDecay = 1)
  rule_select_ite23 <- selectRuleRRF(rule_metric_ite23, X_covariates,ite23)
  save(rule_select_ite23, file = paste0("Intree_ite23_",n_boot_impute,"._Rdata"))
  rf_ite13 <- randomForest(x = X_covariates, y = ite13) 
  tree_list_ite13 <- RF2List(rf_ite13)
  rule_exec_ite13 <- unique(extractRules(tree_list_ite13,X_covariates,digits=3))
  rule_metric_ite13 <- getRuleMetric(rule_exec_ite13,X_covariates,ite13)
  rule_metric_ite13 <- pruneRule(rule_metric_ite13,X_covariates,ite13,typeDecay = 1)
  rule_select_ite13 <- selectRuleRRF(rule_metric_ite13, X_covariates,ite13)
  save(rule_select_ite13, file = paste0("Intree_ite13_",n_boot_impute,"._Rdata"))
}


