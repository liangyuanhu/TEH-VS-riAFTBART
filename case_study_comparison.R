library(tidyverse)
library(riAFTBART)
library(mice)
data_post_ICU_varselect <- read_csv("code/riAFTBART VS/covid_post_icu.csv") %>% 
  mutate_if(is.character, as.factor) %>% 
  select(-treatment) %>% 
  as.data.frame()
# Variable selection using riAFTBART ------------------------------------------------------
n_boot_impute <- 100
maxit_impute<-50
set.seed(123)
mice_b_dat_list <- vector("list",n_boot_impute)
for (i in 1:n_boot_impute){
  # Bootstarp
  id_boot<-sample(c(1:nrow(data_post_ICU_varselect)), nrow(data_post_ICU_varselect), replace=TRUE)
  data_post_ICU_varselect_b<-data_post_ICU_varselect[c(id_boot),]
  mice_b<-mice(data_post_ICU_varselect_b, m=1, maxit = maxit_impute)
  mice_b_list<-mice_b 
  mice_b_dat<-complete(mice_b, include=FALSE)
  mice_b_dat_list[[i]] <- mice_b_dat
}
impute_dat_list <- mice_b_dat_list
## Output list
selected_local_var_i<-list(NA)
selected_globalmax_var_i<-list(NA)
selected_globalse_var_i<-list(NA)

selected_var_i<-list(NA)
cores=detectCores()
cl <- makeCluster(cores-2) #not to overload your computer


registerDoParallel(cl)

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

finalMatrix<-NA
registerDoParallel(cl)
finalMatrix <- foreach(i=1:n_boot_impute, .combine="comb", .multicombine = TRUE,
                       .init=list(list(), list(), list(), list(), list())) %dopar% {
                         library(caret)
                         library(pROC)
                         library(riAFTBART)
                         
                         mice_b_dat <- impute_dat_list[[i]]
                         
                         cov_mat<-mice_b_dat[,-which(colnames(mice_b_dat)%in%c("death","time_icu_to_outcome","LocationName"))]
                         
                         
                         # Sys.time() # 10:25
                         riAFTBART_result <- riAFTBART_fit(M.burnin = 1000,
                                                           M.keep = 3500,
                                                           status = mice_b_dat$death,
                                                           y.train = mice_b_dat$time_icu_to_outcome,
                                                           x.train = cov_mat,
                                                           trt.train = NULL,
                                                           x.test = cov_mat,
                                                           trt.test = NULL,
                                                           cluster.id = as.numeric(as.factor(mice_b_dat$LocationName)),
                                                           verbose = T)
                         
                         var_select_result<-var_select(M.burnin = 1000,
                                                       M.keep = 3500,
                                                       status = mice_b_dat$death,
                                                       y.train = mice_b_dat$time_icu_to_outcome,
                                                       x.train = cov_mat,
                                                       trt.train = NULL,
                                                       x.test = cov_mat,
                                                       trt.test = NULL,
                                                       cluster.id = as.numeric(as.factor(mice_b_dat$LocationName)),
                                                       verbose = T,
                                                       n_permuate = 100 ) 
                         colnames(var_select_result$vip_perm) <- colnames(riAFTBART_result$vip)
                         # H0: combine the probabilities for each factor with different levels, need to be extremely careful here
                         ## Replace () and space in colnames
                         colnames(var_select_result$vip_perm)<-gsub("\\(","_",colnames(var_select_result$vip_perm))
                         colnames(var_select_result$vip_perm)<-gsub("\\)","_",colnames(var_select_result$vip_perm))
                         colnames(var_select_result$vip_perm)<-gsub("\\ ","_",colnames(var_select_result$vip_perm))
                         colnames(var_select_result$vip_perm)<-gsub("\\$","",colnames(var_select_result$vip_perm))
                         colnames(var_select_result$vip_perm)<-gsub("\\,","",colnames(var_select_result$vip_perm))
                         colnames(var_select_result$vip_perm)<-gsub("<","",colnames(var_select_result$vip_perm))
                         colnames(var_select_result$vip_perm)<-gsub("\\=","",colnames(var_select_result$vip_perm))
                         colnames(var_select_result$vip_perm)<-gsub(">","",colnames(var_select_result$vip_perm))
                         colnames(var_select_result$vip_perm)<-gsub("\\/","",colnames(var_select_result$vip_perm))
                         colnames(var_select_result$vip_perm)<-gsub("-","",colnames(var_select_result$vip_perm))
                         colnames(var_select_result$vip_perm)<-gsub("\\,","",colnames(var_select_result$vip_perm))
                         colnames(var_select_result$vip_perm)<-gsub("'","",colnames(var_select_result$vip_perm))
                         colnames(var_select_result$vip_perm)<-gsub("\\+","",colnames(var_select_result$vip_perm))
                         colnames(var_select_result$vip_perm)<-gsub("&","",colnames(var_select_result$vip_perm))
                         var_count<-var_select_result$vip_perm
                         sum_tree<-apply(var_select_result$vip_perm, 1, sum)
                         p_var_selected1<-var_count
                         for(j in 1:ncol(var_count)){
                           var_select_result$vip_perm[,j]<-var_count[,j]/sum_tree
                         }
                         
                         library(dplyr)
                         if(TRUE){
                           
                           vs_permutate_factor<-var_select_result$vip_perm %>%  
                             as_tibble() %>% 
                             rename(
                               gender = gender.Male,
                             )
                           
                         }
                         
                         # Final check: if probabilities sum to 1 for all 100 permutations, then let's celebrate! Otherwise, we need to check whether we've missed any levels or factors....
                         ave_check<-mean(vs_permutate_factor %>% apply(1, sum))
                         
                         
                         # Prop included: combine the probabilities for each factor with different levels, need to be extremely careful here
                         ## Replace () and space in colnames
                         var_count<-riAFTBART_result$vip
                         sum_tree<-apply(var_count, 1, sum)
                         p_var_selected1<-var_count
                         for(j in 1:ncol(var_count)){
                           p_var_selected1[,j]<-var_count[,j]/sum_tree
                         }
                         
                         colnames(p_var_selected1)<-gsub("\\(","_",colnames(p_var_selected1))
                         colnames(p_var_selected1)<-gsub("\\)","_",colnames(p_var_selected1))
                         colnames(p_var_selected1)<-gsub("\\ ","_",colnames(p_var_selected1))
                         colnames(p_var_selected1)<-gsub("\\$","",colnames(p_var_selected1))
                         colnames(p_var_selected1)<-gsub("\\,","",colnames(p_var_selected1))
                         colnames(p_var_selected1)<-gsub("<","",colnames(p_var_selected1))
                         colnames(p_var_selected1)<-gsub("\\=","",colnames(p_var_selected1))
                         colnames(p_var_selected1)<-gsub(">","",colnames(p_var_selected1))
                         colnames(p_var_selected1)<-gsub("\\/","",colnames(p_var_selected1))
                         colnames(p_var_selected1)<-gsub("-","",colnames(p_var_selected1))
                         colnames(p_var_selected1)<-gsub("\\,","",colnames(p_var_selected1))
                         colnames(p_var_selected1)<-gsub("'","",colnames(p_var_selected1))
                         colnames(p_var_selected1)<-gsub("\\+","",colnames(p_var_selected1))
                         colnames(p_var_selected1)<-gsub("&","",colnames(p_var_selected1))
                         
                         library(dplyr)
                         library(tibble)
                         vs_prop_factor0<-p_var_selected1 %>%  
                           as_tibble() %>% 
                           rename(
                             gender = gender.Male,
                           )
                         # Final check: if probabilities sum to 1 for all 100 permutations, then let's celebrate! Otherwise, we need to check whether we've missed any levels or factors....
                         ave_check2<-mean(vs_prop_factor0 %>% apply(1, sum))
                         #if(ave_check2!=1){
                         #   break
                         #}
                         detach("package:dplyr")
                         vs_prop_factor <- vs_prop_factor0 %>% apply(2,mean)
                         # Var seletion (local): if prop > 1-alpha(0.05) quantils in permutation dist
                         ## Caltulate cutoff
                         localcutoff<-apply(vs_permutate_factor, 2, quantile, probs=0.90)
                         
                         ## Merge cutoff and vs_prop_factor by var names 
                         cut_mat<-data.frame("localcutoff"=localcutoff, "varnames"=names(localcutoff))
                         prop_mat<-data.frame("prop"=vs_prop_factor, "varnames"=names(vs_prop_factor))
                         vs_mat<-merge(cut_mat, prop_mat, by="varnames", all=TRUE)
                         
                         vs_mat$local_results<-ifelse(vs_mat$prop>vs_mat$localcutoff,1,0)
                         #table(vs_mat$local_results) # 24 selected
                         vs_mat$varnames[vs_mat$local_results==1]
                         
                         
                         selected_local_var_i<-vs_mat2$varnames[vs_mat2$local_results==1]
                         
                         
                         
                         list(selected_local_var_i, 
                              vs_permutate_factor, vs_prop_factor)
                       }

stopCluster(cl)
save(finalMatrix,  file = "result/BIBART.Rdata")
final_var0<-list(NA)

selected_var_i<-finalMatrix[[1]]
select_freq_i<-table(unlist(selected_var_i))
final_var_i<-boot_finalvar(cutpt_boot_varselect=cutpt_boot_varselect, 
                           varlist=selected_var_i, n_boot=n_boot_impute)

final_var_i_dat<-data.frame(matrix(NA,nrow=length(cutpt_boot_varselect), ncol = 3))
names(final_var_i_dat)<-c("pi", "N_selected", "var_selected")
for(i in 1:length(cutpt_boot_varselect)){
  cut_i<-cutpt_boot_varselect[i]
  final_var_i_dat[i,1]<-cut_i
  final_var_i_dat[i,2]<-length(final_var_i[[i]])
  final_var_i_dat[i,3]<-paste(unlist(final_var_i[[i]]), collapse=", ")
}

# Variable selection using riCox ------------------------------------------------------
n_boot_impute <- 100
maxit_impute<-50
set.seed(123)
mice_b_dat_list <- vector("list",n_boot_impute)
for (i in 1:n_boot_impute){
  # Bootstarp
  id_boot<-sample(c(1:nrow(data_post_ICU_varselect)), nrow(data_post_ICU_varselect), replace=TRUE)
  data_post_ICU_varselect_b<-data_post_ICU_varselect[c(id_boot),]
  mice_b<-mice(data_post_ICU_varselect_b, m=1, maxit = maxit_impute)
  mice_b_list<-mice_b 
  mice_b_dat<-complete(mice_b, include=FALSE)
  mice_b_dat_list[[i]] <- mice_b_dat
}
impute_dat_list <- mice_b_dat_list
## Output list
selected_local_var_i<-list(NA)
selected_globalmax_var_i<-list(NA)
selected_globalse_var_i<-list(NA)

selected_var_i<-list(NA)
cores=detectCores()
cl <- makeCluster(cores-2) #not to overload your computer


registerDoParallel(cl)

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

finalMatrix<-NA
registerDoParallel(cl)
finalMatrix <- foreach(i=1:n_boot_impute, .combine="comb", .multicombine = TRUE,
                       .init=list(list(), list(), list(), list(), list())) %dopar% {
                         library(caret)
                         library(pROC)
                         library(riAFTBART)
                         
                         mice_b_dat <- impute_dat_list[[i]]
                         
                         m1<- try(stepcoxme(data = mice_b_dat, event = "death",time = "time_icu_to_outcome",cl = "LocationName", direction = directionx))
                         if(inherits(m1, "try-error"))
                         {
                           m1 <- data.frame()
                           saveRDS(mice_b_dat, paste("result/cox_i",i,"b",b,"mice_b_dat_debug.rds",sep=""))
                         }
                         
                         
                         if (nrow(m1) == 0){
                           selected_var_i <- NA
                         } else{
                           selected_var_i<-strsplit(m1$Var_in_model[nrow(m1)], split="\\+")[[1]]                          
                         } 
                         
                         # selected_var_b<-m1$var_in_model[[1]]
                         list("selected_var_i"=selected_var_i)
                       }

stopCluster(cl)
save(finalMatrix,  file = "result/riCox.Rdata")
final_var0<-list(NA)

selected_var_i<-finalMatrix[[1]]
select_freq_i<-table(unlist(selected_var_i))
final_var_i<-boot_finalvar(cutpt_boot_varselect=cutpt_boot_varselect, 
                           varlist=selected_var_i, n_boot=n_boot_impute)

final_var_i_dat<-data.frame(matrix(NA,nrow=length(cutpt_boot_varselect), ncol = 3))
names(final_var_i_dat)<-c("pi", "N_selected", "var_selected")
for(i in 1:length(cutpt_boot_varselect)){
  cut_i<-cutpt_boot_varselect[i]
  final_var_i_dat[i,1]<-cut_i
  final_var_i_dat[i,2]<-length(final_var_i[[i]])
  final_var_i_dat[i,3]<-paste(unlist(final_var_i[[i]]), collapse=", ")
}

# Variable selection using FrailtyHL ------------------------------------------------------
n_boot_impute <- 100
maxit_impute<-50
set.seed(123)
mice_b_dat_list <- vector("list",n_boot_impute)
for (i in 1:n_boot_impute){
  # Bootstarp
  id_boot<-sample(c(1:nrow(data_post_ICU_varselect)), nrow(data_post_ICU_varselect), replace=TRUE)
  data_post_ICU_varselect_b<-data_post_ICU_varselect[c(id_boot),]
  mice_b<-mice(data_post_ICU_varselect_b, m=1, maxit = maxit_impute)
  mice_b_list<-mice_b 
  mice_b_dat<-complete(mice_b, include=FALSE)
  mice_b_dat_list[[i]] <- mice_b_dat
}
impute_dat_list <- mice_b_dat_list
## Output list
selected_local_var_i<-list(NA)
selected_globalmax_var_i<-list(NA)
selected_globalse_var_i<-list(NA)

selected_var_i<-list(NA)
cores=detectCores()
cl <- makeCluster(cores-2) #not to overload your computer


registerDoParallel(cl)

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

finalMatrix<-NA
registerDoParallel(cl)
finalMatrix <- foreach(i=1:n_boot_impute, .combine="comb", .multicombine = TRUE,
                       .init=list(list(), list(), list(), list(), list())) %dopar% {
                         library(caret)
                         library(pROC)
                         library(riAFTBART)
                         
                         mice_b_dat <- impute_dat_list[[i]]
                         
                         lasso_result <- frailty.vs(Surv(time_icu_to_outcome,death) ~ age+sex+race+ethnicity+smoking_status+hypertension+coronary_artery_disease+cancer+diabetes+asthma+copd+temperature+sys_blood+dia_blood+oxygen_saturation+heart_rate+fraction_oxygen+bmi+Oxygen_ordinal+sofa_score+glasgow_score+pressure_oxygen+ddimer+ldh+ferritin+crp+creatinine+wbc +(1|LocationName),data = mice_b_dat, model = "lognorm", penalty = "lasso",tun1 = seq(0, 0.1, 0.001), B = lognormal_result$coefficients[-1])
                         scad_result <- frailty.vs(Surv(time_icu_to_outcome,death) ~ age+sex+race+ethnicity+smoking_status+hypertension+coronary_artery_disease+cancer+diabetes+asthma+copd+temperature+sys_blood+dia_blood+oxygen_saturation+heart_rate+fraction_oxygen+bmi+Oxygen_ordinal+sofa_score+glasgow_score+pressure_oxygen+ddimer+ldh+ferritin+crp+creatinine+wbc +(1|LocationName),data = mice_b_dat, model = "lognorm", penalty = "scad",tun1 = seq(0, 0.1, 0.001), B = lognormal_result$coefficients[-1])
                         hl_result <- frailty.vs(Surv(time_icu_to_outcome,death) ~ age+sex+race+ethnicity+smoking_status+hypertension+coronary_artery_disease+cancer+diabetes+asthma+copd+temperature+sys_blood+dia_blood+oxygen_saturation+heart_rate+fraction_oxygen+bmi+Oxygen_ordinal+sofa_score+glasgow_score+pressure_oxygen+ddimer+ldh+ferritin+crp+creatinine+wbc+(1|LocationName),data = mice_b_dat, model = "lognorm", penalty = "hl",tun1 = seq(0, 0.1, 0.001), tun2=seq(0.001,0.25, 0.001),B = lognormal_result$coefficients[-1])
                         
                         min_BIC_index <- which.min(c(lasso_result$BIC, scad_result$BIC, hl_result$BIC))
                         if (min_BIC_index == 1)  {
                           list("selected_var_i"= names(mice_b_dat)[4:31][lasso_result$beta!=0])
                         }
                         if (min_BIC_index == 2)  {
                           list("selected_var_i"= names(mice_b_dat)[4:31][scad_result$beta!=0])
                         }
                         if (min_BIC_index == 3)  {
                           list("selected_var_i"= names(mice_b_dat)[4:31][hl_result$beta!=0])
                         }
                       }

stopCluster(cl)
save(finalMatrix,  file = "result/FrailtyHL.Rdata")
final_var0<-list(NA)

selected_var_i<-finalMatrix[[1]]
select_freq_i<-table(unlist(selected_var_i))
final_var_i<-boot_finalvar(cutpt_boot_varselect=cutpt_boot_varselect, 
                           varlist=selected_var_i, n_boot=n_boot_impute)

final_var_i_dat<-data.frame(matrix(NA,nrow=length(cutpt_boot_varselect), ncol = 3))
names(final_var_i_dat)<-c("pi", "N_selected", "var_selected")
for(i in 1:length(cutpt_boot_varselect)){
  cut_i<-cutpt_boot_varselect[i]
  final_var_i_dat[i,1]<-cut_i
  final_var_i_dat[i,2]<-length(final_var_i[[i]])
  final_var_i_dat[i,3]<-paste(unlist(final_var_i[[i]]), collapse=", ")
}
# Variable selection using Frailtypack ------------------------------------------------------
n_boot_impute <- 100
maxit_impute<-50
set.seed(123)
mice_b_dat_list <- vector("list",n_boot_impute)
for (i in 1:n_boot_impute){
  # Bootstarp
  id_boot<-sample(c(1:nrow(data_post_ICU_varselect)), nrow(data_post_ICU_varselect), replace=TRUE)
  data_post_ICU_varselect_b<-data_post_ICU_varselect[c(id_boot),]
  mice_b<-mice(data_post_ICU_varselect_b, m=1, maxit = maxit_impute)
  mice_b_list<-mice_b 
  mice_b_dat<-complete(mice_b, include=FALSE)
  mice_b_dat_list[[i]] <- mice_b_dat
}
impute_dat_list <- mice_b_dat_list
## Output list
selected_local_var_i<-list(NA)
selected_globalmax_var_i<-list(NA)
selected_globalse_var_i<-list(NA)

selected_var_i<-list(NA)
cores=detectCores()
cl <- makeCluster(cores-2) #not to overload your computer


registerDoParallel(cl)

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

finalMatrix<-NA
registerDoParallel(cl)
finalMatrix <- foreach(i=1:n_boot_impute, .combine="comb", .multicombine = TRUE,
                       .init=list(list(), list(), list(), list(), list())) %dopar% {
                         library(caret)
                         library(pROC)
                         library(riAFTBART)
                         
                         mice_b_dat <- impute_dat_list[[i]]
                         
                         m1 <- frailtyPenal(Surv(time_icu_to_outcome,death) ~ cluster(LocationName) + age+sex+race+ethnicity+smoking_status+hypertension+coronary_artery_disease+cancer+diabetes+asthma+copd+temperature+sys_blood+dia_blood+oxygen_saturation+heart_rate+fraction_oxygen+bmi+Oxygen_ordinal+sofa_score+glasgow_score+pressure_oxygen+ddimer+ldh+ferritin+crp+creatinine+wbc, n.knots = 12, kappa = 10, data=  mice_b_dat, cross.validation=TRUE)
                         selected_var_i<-names(abs(m1$coef) > 0.001) 
                       }

stopCluster(cl)
save(finalMatrix,  file = "result/frailtypack.Rdata")
final_var0<-list(NA)

selected_var_i<-finalMatrix[[1]]
select_freq_i<-table(unlist(selected_var_i))
final_var_i<-boot_finalvar(cutpt_boot_varselect=cutpt_boot_varselect, 
                           varlist=selected_var_i, n_boot=n_boot_impute)

final_var_i_dat<-data.frame(matrix(NA,nrow=length(cutpt_boot_varselect), ncol = 3))
names(final_var_i_dat)<-c("pi", "N_selected", "var_selected")
for(i in 1:length(cutpt_boot_varselect)){
  cut_i<-cutpt_boot_varselect[i]
  final_var_i_dat[i,1]<-cut_i
  final_var_i_dat[i,2]<-length(final_var_i[[i]])
  final_var_i_dat[i,3]<-paste(unlist(final_var_i[[i]]), collapse=", ")
}
# Variable selection using PEAMM ------------------------------------------------------
n_boot_impute <- 100
maxit_impute<-50
set.seed(123)
mice_b_dat_list <- vector("list",n_boot_impute)
for (i in 1:n_boot_impute){
  # Bootstarp
  id_boot<-sample(c(1:nrow(data_post_ICU_varselect)), nrow(data_post_ICU_varselect), replace=TRUE)
  data_post_ICU_varselect_b<-data_post_ICU_varselect[c(id_boot),]
  mice_b<-mice(data_post_ICU_varselect_b, m=1, maxit = maxit_impute)
  mice_b_list<-mice_b 
  mice_b_dat<-complete(mice_b, include=FALSE)
  mice_b_dat_list[[i]] <- mice_b_dat
}
impute_dat_list <- mice_b_dat_list
## Output list
selected_local_var_i<-list(NA)
selected_globalmax_var_i<-list(NA)
selected_globalse_var_i<-list(NA)

selected_var_i<-list(NA)
cores=detectCores()
cl <- makeCluster(cores-2) #not to overload your computer


registerDoParallel(cl)

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

finalMatrix<-NA
registerDoParallel(cl)
finalMatrix <- foreach(i=1:n_boot_impute, .combine="comb", .multicombine = TRUE,
                       .init=list(list(), list(), list(), list(), list())) %dopar% {
                         library(caret)
                         library(pROC)
                         library(riAFTBART)
                         
                         mice_b_dat <- impute_dat_list[[i]]
                         
                         peamm_penalize(data = mice_b_dat, event = "death",time = "time_icu_to_outcome",cl = "cl") 
                       }

stopCluster(cl)
save(finalMatrix,  file = "result/PEAMM.Rdata")
final_var0<-list(NA)

selected_var_i<-finalMatrix[[1]]
select_freq_i<-table(unlist(selected_var_i))
final_var_i<-boot_finalvar(cutpt_boot_varselect=cutpt_boot_varselect, 
                           varlist=selected_var_i, n_boot=n_boot_impute)

final_var_i_dat<-data.frame(matrix(NA,nrow=length(cutpt_boot_varselect), ncol = 3))
names(final_var_i_dat)<-c("pi", "N_selected", "var_selected")
for(i in 1:length(cutpt_boot_varselect)){
  cut_i<-cutpt_boot_varselect[i]
  final_var_i_dat[i,1]<-cut_i
  final_var_i_dat[i,2]<-length(final_var_i[[i]])
  final_var_i_dat[i,3]<-paste(unlist(final_var_i[[i]]), collapse=", ")
}
# Concordence index comparison -------------------------------------------------------------------
library(riAFTBART)
library(SurvMetrics)
library(survival)
library(coxme)
data_post_icu <- read_csv("code/riAFTBART VS/covid_post_icu.csv") %>% 
  mutate_if(is.character, as.factor) %>% 
  select(-treatment) %>% 
  as.data.frame()
n_boot_impute <- 250
maxit_impute<-50
set.seed(123)
mice_b_dat_list <- vector("list",n_boot_impute)
impute_dat_list <- mice_b_dat_list
metrics_riAFTBART <- rep(NA, n_boot_impute)
metrics_riCox <- rep(NA, n_boot_impute)
metrics_FrailtyHL <- rep(NA, n_boot_impute)
metrics_Frailtypack <- rep(NA, n_boot_impute)
metrics_PEAMM <- rep(NA, n_boot_impute)
for (i in 1:n_boot_impute){
  # Bootstarp
  index_data = createFolds(1:nrow(mice_b_dat_list_once), 2)
  train_data = mice_b_dat_list_once[index_data[[1]],]
  test_data = mice_b_dat_list_once[index_data[[2]],] 
  id_boot<-sample(c(1:nrow(train_data)), nrow(train_data), replace=TRUE)
  data_post_ICU_varselect_b<-train_data[c(id_boot),]
  mice_b<-mice(data_post_ICU_varselect_b, m=1, maxit = maxit_impute)
  mice_b_list<-mice_b 
  mice_b_dat<-complete(mice_b, include=FALSE)
  mice_b_dat_list[[i]] <- mice_b_dat
  mice_b_dat <- impute_dat_list[[i]]
  
  riAFTBART_result <- riAFTBART_fit(M.burnin = 1000,
                                    M.keep = 3500,
                                    status = mice_b_dat$death,
                                    y.train = mice_b_dat$time_icu_to_outcome,
                                    x.train = cov_mat,
                                    trt.train = NULL,
                                    x.test = cov_mat,
                                    trt.test = NULL,
                                    cluster.id = as.numeric(as.factor(mice_b_dat$LocationName)),
                                    verbose = T)
  
  var_select_result<-var_select(M.burnin = 1000,
                                M.keep = 3500,
                                status = mice_b_dat$death,
                                y.train = mice_b_dat$time_icu_to_outcome,
                                x.train = cov_mat,
                                trt.train = NULL,
                                x.test = cov_mat,
                                trt.test = NULL,
                                cluster.id = as.numeric(as.factor(mice_b_dat$LocationName)),
                                verbose = T,
                                n_permuate = 100 ) 
  colnames(var_select_result$vip_perm) <- colnames(riAFTBART_result$vip)
  # H0: combine the probabilities for each factor with different levels, need to be extremely careful here
  ## Replace () and space in colnames
  colnames(var_select_result$vip_perm)<-gsub("\\(","_",colnames(var_select_result$vip_perm))
  colnames(var_select_result$vip_perm)<-gsub("\\)","_",colnames(var_select_result$vip_perm))
  colnames(var_select_result$vip_perm)<-gsub("\\ ","_",colnames(var_select_result$vip_perm))
  colnames(var_select_result$vip_perm)<-gsub("\\$","",colnames(var_select_result$vip_perm))
  colnames(var_select_result$vip_perm)<-gsub("\\,","",colnames(var_select_result$vip_perm))
  colnames(var_select_result$vip_perm)<-gsub("<","",colnames(var_select_result$vip_perm))
  colnames(var_select_result$vip_perm)<-gsub("\\=","",colnames(var_select_result$vip_perm))
  colnames(var_select_result$vip_perm)<-gsub(">","",colnames(var_select_result$vip_perm))
  colnames(var_select_result$vip_perm)<-gsub("\\/","",colnames(var_select_result$vip_perm))
  colnames(var_select_result$vip_perm)<-gsub("-","",colnames(var_select_result$vip_perm))
  colnames(var_select_result$vip_perm)<-gsub("\\,","",colnames(var_select_result$vip_perm))
  colnames(var_select_result$vip_perm)<-gsub("'","",colnames(var_select_result$vip_perm))
  colnames(var_select_result$vip_perm)<-gsub("\\+","",colnames(var_select_result$vip_perm))
  colnames(var_select_result$vip_perm)<-gsub("&","",colnames(var_select_result$vip_perm))
  var_count<-var_select_result$vip_perm
  sum_tree<-apply(var_select_result$vip_perm, 1, sum)
  p_var_selected1<-var_count
  for(j in 1:ncol(var_count)){
    var_select_result$vip_perm[,j]<-var_count[,j]/sum_tree
  }
  
  library(dplyr)
  if(TRUE){
    
    vs_permutate_factor<-var_select_result$vip_perm %>%  
      as_tibble() %>% 
      rename(
        gender = gender.Male,
      )
    
  }
  
  # Final check: if probabilities sum to 1 for all 100 permutations, then let's celebrate! Otherwise, we need to check whether we've missed any levels or factors....
  ave_check<-mean(vs_permutate_factor %>% apply(1, sum))
  
  
  # Prop included: combine the probabilities for each factor with different levels, need to be extremely careful here
  ## Replace () and space in colnames
  var_count<-riAFTBART_result$vip
  sum_tree<-apply(var_count, 1, sum)
  p_var_selected1<-var_count
  for(j in 1:ncol(var_count)){
    p_var_selected1[,j]<-var_count[,j]/sum_tree
  }
  
  colnames(p_var_selected1)<-gsub("\\(","_",colnames(p_var_selected1))
  colnames(p_var_selected1)<-gsub("\\)","_",colnames(p_var_selected1))
  colnames(p_var_selected1)<-gsub("\\ ","_",colnames(p_var_selected1))
  colnames(p_var_selected1)<-gsub("\\$","",colnames(p_var_selected1))
  colnames(p_var_selected1)<-gsub("\\,","",colnames(p_var_selected1))
  colnames(p_var_selected1)<-gsub("<","",colnames(p_var_selected1))
  colnames(p_var_selected1)<-gsub("\\=","",colnames(p_var_selected1))
  colnames(p_var_selected1)<-gsub(">","",colnames(p_var_selected1))
  colnames(p_var_selected1)<-gsub("\\/","",colnames(p_var_selected1))
  colnames(p_var_selected1)<-gsub("-","",colnames(p_var_selected1))
  colnames(p_var_selected1)<-gsub("\\,","",colnames(p_var_selected1))
  colnames(p_var_selected1)<-gsub("'","",colnames(p_var_selected1))
  colnames(p_var_selected1)<-gsub("\\+","",colnames(p_var_selected1))
  colnames(p_var_selected1)<-gsub("&","",colnames(p_var_selected1))
  
  library(dplyr)
  library(tibble)
  vs_prop_factor0<-p_var_selected1 %>%  
    as_tibble() %>% 
    rename(
      gender = gender.Male,
    )
  # Final check: if probabilities sum to 1 for all 100 permutations, then let's celebrate! Otherwise, we need to check whether we've missed any levels or factors....
  ave_check2<-mean(vs_prop_factor0 %>% apply(1, sum))
  #if(ave_check2!=1){
  #   break
  #}
  detach("package:dplyr")
  vs_prop_factor <- vs_prop_factor0 %>% apply(2,mean)
  # Var seletion (local): if prop > 1-alpha(0.05) quantils in permutation dist
  ## Caltulate cutoff
  localcutoff<-apply(vs_permutate_factor, 2, quantile, probs=0.90)
  
  ## Merge cutoff and vs_prop_factor by var names 
  cut_mat<-data.frame("localcutoff"=localcutoff, "varnames"=names(localcutoff))
  prop_mat<-data.frame("prop"=vs_prop_factor, "varnames"=names(vs_prop_factor))
  vs_mat<-merge(cut_mat, prop_mat, by="varnames", all=TRUE)
  
  vs_mat$local_results<-ifelse(vs_mat$prop>vs_mat$localcutoff,1,0)
  #table(vs_mat$local_results) # 24 selected
  vs_mat$varnames[vs_mat$local_results==1]
  
  
  selected_local_var_riaftbart<-vs_mat2$varnames[vs_mat2$local_results==1]
  
  stepcoxme(data = mice_b_dat, event = "death",time = "time_icu_to_outcome",cl = "LocationName", direction = directionx)
  
  selected_var_i_ricox<-strsplit(m1$Var_in_model[nrow(m1)], split="\\+")[[1]]   
  lasso_result <- frailty.vs(Surv(time_icu_to_outcome,death) ~ age+sex+race+ethnicity+smoking_status+hypertension+coronary_artery_disease+cancer+diabetes+asthma+copd+temperature+sys_blood+dia_blood+oxygen_saturation+heart_rate+fraction_oxygen+bmi+Oxygen_ordinal+sofa_score+glasgow_score+pressure_oxygen+ddimer+ldh+ferritin+crp+creatinine+wbc +(1|LocationName),data = mice_b_dat, model = "lognorm", penalty = "lasso",tun1 = seq(0, 0.1, 0.001), B = lognormal_result$coefficients[-1])
  scad_result <- frailty.vs(Surv(time_icu_to_outcome,death) ~ age+sex+race+ethnicity+smoking_status+hypertension+coronary_artery_disease+cancer+diabetes+asthma+copd+temperature+sys_blood+dia_blood+oxygen_saturation+heart_rate+fraction_oxygen+bmi+Oxygen_ordinal+sofa_score+glasgow_score+pressure_oxygen+ddimer+ldh+ferritin+crp+creatinine+wbc +(1|LocationName),data = mice_b_dat, model = "lognorm", penalty = "scad",tun1 = seq(0, 0.1, 0.001), B = lognormal_result$coefficients[-1])
  hl_result <- frailty.vs(Surv(time_icu_to_outcome,death) ~ age+sex+race+ethnicity+smoking_status+hypertension+coronary_artery_disease+cancer+diabetes+asthma+copd+temperature+sys_blood+dia_blood+oxygen_saturation+heart_rate+fraction_oxygen+bmi+Oxygen_ordinal+sofa_score+glasgow_score+pressure_oxygen+ddimer+ldh+ferritin+crp+creatinine+wbc+(1|LocationName),data = mice_b_dat, model = "lognorm", penalty = "hl",tun1 = seq(0, 0.1, 0.001), tun2=seq(0.001,0.25, 0.001),B = lognormal_result$coefficients[-1])
  
  min_BIC_index <- which.min(c(lasso_result$BIC, scad_result$BIC, hl_result$BIC))
  if (min_BIC_index == 1)  {
    selected_var_i_frailtyhl= names(mice_b_dat)[4:31][lasso_result$beta!=0]
  }
  if (min_BIC_index == 2)  {
    selected_var_i_frailtyhl= names(mice_b_dat)[4:31][scad_result$beta!=0]
  }
  if (min_BIC_index == 3)  {
    selected_var_i_frailtyhl= names(mice_b_dat)[4:31][hl_result$beta!=0]
  }
  m1 <- frailtyPenal(Surv(time_icu_to_outcome,death) ~ cluster(LocationName) + age+sex+race+ethnicity+smoking_status+hypertension+coronary_artery_disease+cancer+diabetes+asthma+copd+temperature+sys_blood+dia_blood+oxygen_saturation+heart_rate+fraction_oxygen+bmi+Oxygen_ordinal+sofa_score+glasgow_score+pressure_oxygen+ddimer+ldh+ferritin+crp+creatinine+wbc, n.knots = 12, kappa = 10, data=  mice_b_dat, cross.validation=TRUE)
  selected_var_i_frailtypack<-names(abs(m1$coef) > 0.001) 
  selected_var_i_peamm <- peamm_penalize(data = mice_b_dat, event = "death",time = "time_icu_to_outcome",cl = "cl") 
  riAFTBART_model <- riAFTBART_fit(M.burnin = 1000,
                                   M.keep = 3500,
                                   status = test_data$death,
                                   y.train = test_data$time_icu_to_outcome,
                                   x.train = test_data %>% 
                                     select(d_dimer, Creatinine, Oxygen_ordinal, age,wbc, Race, sofa_score, ldh) %>%
                                     mutate_if(is.character, as.factor) %>% 
                                     as.data.frame(),
                                   trt.train = NULL,
                                   x.test = test_data %>% 
                                     select(d_dimer, Creatinine, Oxygen_ordinal, age,wbc, Race, sofa_score, ldh) %>%
                                     mutate_if(is.character, as.factor) %>% 
                                     as.data.frame(),
                                   trt.test = NULL,
                                   cluster.id = test_data$LocationName,
                                   verbose = T)
  
  riCox <- coxme(Surv(time_icu_to_outcome, death)~ age + Oxygen_ordinal+creatinine+wbc+sofa_score+race+ddimer+ldh + (1|LocationName), test_data)
  FrailtyHL <- coxme(Surv(time_icu_to_outcome, death)~ age + Oxygen_ordinal+creatinine+wbc+sofa_score+race+ddimer+ldh + (1|LocationName), test_data)
  Frailtypack <- coxme(Surv(time_icu_to_outcome, death)~ age + Oxygen_ordinal+creatinine+wbc+sofa_score+race+ddimer+ldh + (1|LocationName), test_data)
  PEAMM_pred <- test_data %>% as_ped(Surv(time_icu_to_outcome, death)~ age + Oxygen_ordinal+creatinine+wbc+sofa_score+race+ddimer+ldh + LocationName)
  PEAMM <- gam(ped_status ~ s(tend) + age + Oxygen_ordinal+creatinine+wbc+sofa_score+race+ddimer+ldh, data = PEAMM_pred,family = poisson(), offset = offset)
  mat_riAFTBART = predictSurvProb(riAFTBART_model, test_data, test_data$time_icu_to_outcome)
  mat_riCox = predictSurvProb(riCox, test_data, test_data$time_icu_to_outcome)
  mat_FrailtyHL = predictSurvProb(FrailtyHL, test_data, test_data$time_icu_to_outcome)
  mat_Frailtypack = predictSurvProb(Frailtypack, test_data, test_data$time_icu_to_outcome)
  mat_PEAMM = predictSurvProb(PEAMM, test_data, test_data$time_icu_to_outcome)
  surv_obj = Surv(test_data$time_icu_to_outcome, test_data$death)
  metrics_riAFTBART_once = Cindex(surv_obj, sp_matrix = mat_riAFTBART, test_data$time_icu_to_outcome)
  metrics_riCox_once = Cindex(surv_obj, sp_matrix = mat_riCox, test_data$time_icu_to_outcome)
  metrics_FrailtyHL_once = Cindex(surv_obj, sp_matrix = mat_FrailtyHL, test_data$time_icu_to_outcome)
  metrics_Frailtypack_once = Cindex(surv_obj, sp_matrix = mat_Frailtypack, test_data$time_icu_to_outcome)
  metrics_PEAMM_once = Cindex(surv_obj, sp_matrix = mat_PEAMM, test_data$time_icu_to_outcome)
  metrics_riAFTBART[i] <- metrics_riAFTBART_once
  metrics_riCox[i] <- metrics_riCox_once
  metrics_FrailtyHL[i] <- metrics_FrailtyHL_once
  metrics_Frailtypack[i] <- metrics_Frailtypack_once
  metrics_PEAMM[i] <- metrics_PEAMM_once
}

