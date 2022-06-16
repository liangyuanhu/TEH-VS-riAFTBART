library(foreach)
library(doParallel)
library(doSNOW)
library(pammtools)
library(mgcv)
library(Metrics)
library(survival)
library(tidyverse)
source("code/data_gen_var_select.R")
n_rep<-250 # Num of replications (250 to be)
n_boot_impute<-100 # Num of bootstraps to be imputed (100 to be)
n_impute<-1 # Num of imputation (1 to be)
maxit_impute<-5 # (5 to be)

directionx<-"back_for"

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

peamm_penalize <- function(data, event,time,cl){
  varname <- names(mice_b_dat)[-which(names(mice_b_dat)%in%c("time","event","cl"))]
  mice_b_dat_ped <- data %>% as_ped(Surv(time, event)~ ., id = "id")
  var_list<-varname
  fm_i <- paste("ped_status~tend+", paste0(var_list, collapse = "+"),  sep = "")
  riGAPH_mydata <- mgcv::gam(as.formula(fm_i), data = mice_b_dat_ped, family = poisson(), offset = offset, select = T)
  Var_in_model <- names(riGAPH_mydata$coefficients)[-1][abs(riGAPH_mydata$coefficients)[-1] > 0]
}

cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
impute_boot<-foreach(i=1:n_rep, .combine="comb", .multicombine = TRUE,  .init=list(list()), 
                     .packages = c("doParallel", "doSNOW")) %:%
  
  foreach(b=1:n_boot_impute, .combine="comb", .multicombine = TRUE, .init=list(list()),
          .packages = c("doParallel")) %dopar%{
            library(caret)
            library(pROC)
            library(glmnet)
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
            mice_b_dat %>% names
            # Model
            m1<- try(peamm_penalize(data = mice_b_dat, event = "event",time = "time",cl = "cl"))
            if(inherits(m1, "try-error"))
            {
              m1 <- data.frame()
              saveRDS(mice_b_dat, paste("result/gaph_i",i,"b",b,"mice_b_dat_debug.rds",sep=""))
            }
            
            
            # if (nrow(m1) == 0){
            #   selected_var_i <- NA
            # } else{
            #   selected_var_i<-strsplit(m1$Var_in_model[nrow(m1)], split="\\+")[[1]]                          
            # } 
            
            # selected_var_b<-m1$var_in_model[[1]]
            
            selected_var_i <- m1$Var_in_model
            list("selected_var_i"=selected_var_i)
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

saveRDS(final_results,  paste("result/final_results_gaph.rds",sep=""))


## Organize results
final_var0<-finalMatrix[[1]]
prec_recall_f1_i<-finalMatrix[[2]]
metrics_mat0<-finalMatrix[[3]]


metrics_mat_list<-list(NA)
for(i in 1:length(cutpt_boot_varselect)){
  # Precision, recall, f1
  metrics_mat_list_i0<- lapply(metrics_mat0, function(x){x[i,]})
  metrics_mat_list_i<-data.frame(matrix(unlist(metrics_mat_list_i0), nrow=n_rep, byrow=T))
  names(metrics_mat_list_i)<-c("Cutpt","Precision","Recall","F1")
  
  ## Power & type I error
  final_var<-lapply(final_var0, function(x){x[[i]]})
  final_var_factor<-factor(unlist(final_var), 
                           levels = c("x1","x2","x3","x4","x5","x6","x7","x8",
                                      "z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                      "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20"))
  metrics_mat_list_i[,c("power_x1","power_x2","power_x3","power_x4","power_x5",
                        "power_x6","power_x7","power_x8",
                        "type1err_z1","type1err_z2","type1err_z3","type1err_z4","type1err_z5","type1err_z6","type1err_z7","type1err_z8","type1err_z9","type1err_z10",
                        "type1err_z11","type1err_z12","type1err_z13","type1err_z14","type1err_z15","type1err_z16","type1err_z17","type1err_z18","type1err_z19","type1err_z20")]<-
    matrix(rep(table(final_var_factor)/n_rep, n_rep), ncol=50, byrow = T)
  metrics_mat_list_i$power_overall<-apply(metrics_mat_list_i[,c("power_x1","power_x2","power_x3","power_x4","power_x5",
                                                                "power_x6","power_x7","power_x8")],1,mean)
  metrics_mat_list_i$type1err_overall<-apply(metrics_mat_list_i[,c("type1err_z1","type1err_z2","type1err_z3","type1err_z4","type1err_z5","type1err_z6",
                                                                   "type1err_z7","type1err_z8","type1err_z9","type1err_z10",
                                                                   "type1err_z11","type1err_z12","type1err_z13","type1err_z14","type1err_z15","type1err_z16","type1err_z17","type1err_z18","type1err_z19","type1err_z20")],1,mean)
  
  
  
  
  metrics_mat_list[[i]]<-metrics_mat_list_i
}
saveRDS(metrics_mat_list, paste("result/metrics_mat_gaph_boot.rds",sep=""))
