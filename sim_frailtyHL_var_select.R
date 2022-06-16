library(foreach)
library(doParallel)
library(doSNOW)
library(Metrics)
library(survival)
library(tidyverse)
library(frailtyHL)
source("data_gen_var_select.R")

n_rep<-10 # Num of replications (250 to be)
n_boot_impute<-10 # Num of bootstraps to be imputed (100 to be)
n_impute<-1 # Num of imputation (1 to be)
maxit_impute<-5 # (5 to be)

cutpt_boot_varselect<-c(1/n_boot_impute, seq(0.1,1.0,0.1))

# Output
comp_dat_list<-list(NA)
var_select_model<-list(NA)
final_var<-list(NA)

#setup parallel backend to use many processors
cores=detectCores()

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}
library(mice)
comp_dat_list<-list(NA)
na_dat_list<-list(NA)
set.seed(12)
for(i in 1:n_rep) {
  dat_simulation <- dat_gen(n = 2000, K = 10, nK = 200, PH=TRUE)
  comp_dat_list[[i]] <- dat_simulation$dat_comp
  na_dat_list[[i]] <- dat_simulation$dat_na
  print(i)
}

cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
impute_boot<-foreach(i=1:n_rep, .combine="comb", .multicombine = TRUE,  .init=list(list()), 
                     .packages = c("doParallel", "doSNOW")) %:%
  
  foreach(b=1:n_boot_impute, .combine="comb", .multicombine = TRUE, .init=list(list()),
          .packages = c("doParallel")) %dopar%{
            library(caret)
            library(pROC)
            library(timereg)
            library(mice)
            library(missForest)
            library(car)
            dat_na<-na_dat_list[[i]]
            # Bootstarp
            id_boot<-sample(c(1:nrow(dat_na)), nrow(dat_na), replace=TRUE)
            dat_b<-dat_na[c(id_boot),]
            # Impute
            missfrst_b<-missForest(dat_b)
            mice_b_dat<-missfrst_b$ximp
            
            # mice_b<-mice(dat_b, m=1, maxit = 5)
            # mice_b_dat<-complete(mice_b, include=FALSE)
            names(mice_b_dat)[1:2] <- c("time", "event")
            mice_b_dat
            # Model
            lognormal_result <- survival::survreg(survival::Surv(time,event) ~ x1+x2+x3+x4+x5+x6+x7+x8+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12+z13+z14+z15+z16+z17+z18+z19+z20, dist="lognormal", data = mice_b_dat)
            
            lasso_result <- frailty.vs(Surv(time,event) ~ x1+x2+x3+x4+x5+x6+x7+x8+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12+z13+z14+z15+z16+z17+z18+z19+z20 +(1|cl),data = mice_b_dat, model = "lognorm", penalty = "lasso",tun1 = seq(0, 0.1, 0.001), B = lognormal_result$coefficients[-1])
            scad_result <- frailty.vs(Surv(time,event) ~ x1+x2+x3+x4+x5+x6+x7+x8+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12+z13+z14+z15+z16+z17+z18+z19+z20 +(1|cl),data = mice_b_dat, model = "lognorm", penalty = "scad",tun1 = seq(0, 0.1, 0.001), B = lognormal_result$coefficients[-1])
            hl_result <- frailty.vs(Surv(time,event) ~ x1+x2+x3+x4+x5+x6+x7+x8+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12+z13+z14+z15+z16+z17+z18+z19+z20 +(1|cl),data = mice_b_dat, model = "lognorm", penalty = "hl",tun1 = seq(0, 0.1, 0.001), tun2=seq(0.001,0.25, 0.001),B = lognormal_result$coefficients[-1])
            
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
timestamp()

## Metric
final_var0<-list(NA)
prec_recall_f1_met<-list(NA)
metrics_mat0<-list(NA)
final_number_fitted_model <- list(NA)
boot_finalvar<-function(cutpt_boot_varselect, varlist, n_boot){
  final_var<-list(NA)
  met_max<-list(NA)
  
  for(j in 1:length(cutpt_boot_varselect)){
    boot_cut_j<-cutpt_boot_varselect[j]
    varlist_i<-table(unlist(varlist))
    # Identify final var by cut_j
    final_var[[j]]<-names(varlist_i[varlist_i>=boot_cut_j*n_boot])
  }
  return(final_var)
}

prec_recall_f1<-function(gold_pos,gold_neg,var_selected){
  tp<-intersect(gold_pos,var_selected); tp_n<-length(tp)
  fp<-intersect(gold_neg,var_selected); fp_n<-length(fp)
  fn<-setdiff(gold_pos,var_selected); fn_n<-length(fn)
  precision<-tp_n/(tp_n+fp_n)
  recall<-tp_n/(tp_n+fn_n)
  f1<-2*precision*recall/(precision+recall)
  return(list("tp"=tp,"fp"=fp,"fn"=fn,"precision"=precision,"recall"=recall,"f1"=f1))
}

for(i in 1:n_rep) {
  selected_var_i<-impute_boot[[1]][[i]]
  number_fitted_model <- (1-sum(is.na(unlist(selected_var_i)))/n_boot_impute)
  select_freq_i<-table(unlist(selected_var_i))
  final_var_i<-boot_finalvar(cutpt_boot_varselect=cutpt_boot_varselect, 
                             varlist=selected_var_i, n_boot=n_boot_impute)
  
  # 5. Metrics
  prec_recall_f1_met_i<-list(NA)
  metrics_mat_i<-data.frame(matrix(NA, nrow=length(cutpt_boot_varselect), ncol=4))
  names(metrics_mat_i)<-c("Cutpt","Precision","Recall","F1")
  
  for(j in 1:length(cutpt_boot_varselect)){
    final_var_j<-final_var_i[[j]]
    ## Precision, Recall, F1
    prec_recall_f1_met_i[[j]]<-prec_recall_f1(gold_pos=c("x1","x2","x3","x4","x5","x6","x7","x8"),
                                              gold_neg=c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                                         "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20"),
                                              var_selected=final_var_j)
    
    ## Precision, recall f1
    metrics_mat_i[j,]<-c(cutpt_boot_varselect[j], prec_recall_f1_met_i[[j]]$precision, 
                         prec_recall_f1_met_i[[j]]$recall, prec_recall_f1_met_i[[j]]$f1)
    
  }
  final_var0[[i]]<-final_var_i
  metrics_mat0[[i]]<-metrics_mat_i
  prec_recall_f1_met[[i]]<-prec_recall_f1_met_i
  final_number_fitted_model[[i]] <- number_fitted_model
}
finalMatrix<-list("final_var0"=final_var0, "prec_recall_f1_met"=prec_recall_f1_met, "metrics_mat0"=metrics_mat0,
                  "final_number_fitted_model" = final_number_fitted_model)
final_results<-list("dat"=list(comp_dat_list, na_dat_list), "impute_boot"=impute_boot, "met"=finalMatrix)

saveRDS(final_results,  paste("final_results_frailtyHL.rds"))