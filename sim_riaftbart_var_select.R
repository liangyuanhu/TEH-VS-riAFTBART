library(foreach)
library(doParallel)
# library(doSNOW) 
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
# code from bartmachine R package for variable selection using global SE method
bisectK <- function(tol, coverage, permute_mat, x_left, x_right, countLimit, perm_mean, perm_se){
  count = 0
  guess = mean(c(x_left, x_right))
  while ((x_right - x_left) / 2 >= tol & count < countLimit){
    empirical_coverage = mean(sapply(1 : nrow(permute_mat),
                                     function(s){all(permute_mat[s,] - perm_mean <= guess * perm_se)}))
    if (empirical_coverage - coverage == 0){
      break
    } else if (empirical_coverage - coverage < 0){
      x_left = guess
    } else {
      x_right = guess
    }
    guess = mean(c(x_left, x_right))
    count = count + 1
  }
  guess
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
            library(riAFTBART)
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
            cov_mat <- mice_b_dat[,-which(colnames(mice_b_dat)%in%c("event","time","cl"))]
            riAFTBART_result <- riAFTBART_fit(M.burnin = 2000,
                                              M.keep = 3500,
                                              status = mice_b_dat$event,
                                              y.train = mice_b_dat$time,
                                              x.train = cov_mat,
                                              trt.train = NULL,
                                              x.test = cov_mat,
                                              trt.test = NULL,
                                              cluster.id = as.numeric(as.factor(mice_b_dat$cl)),
                                              verbose = T)
            
            var_select_result<-var_select(M.burnin = 2000,
                                          M.keep = 3500,
                                          status = mice_b_dat$event,
                                          y.train = mice_b_dat$time,
                                          x.train = cov_mat,
                                          trt.train = NULL,
                                          x.test = cov_mat,
                                          trt.test = NULL,
                                          cluster.id = as.numeric(as.factor(mice_b_dat$cl)),
                                          verbose = T,
                                          n_permuate = 2 ) 
            colnames(var_select_result$vip_perm) <- colnames(riAFTBART_result$vip)
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
                as_tibble() 
              
            }
            
            # Final check: if probabilities sum to 1 for all 100 permutations, then let's celebrate! Otherwise, we need to check whether we've missed any levels or factors....
            ave_check<-mean(vs_permutate_factor %>% apply(1, sum))
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
              as_tibble()
            # Final check: if probabilities sum to 1 for all 100 permutations, then let's celebrate! Otherwise, we need to check whether we've missed any levels or factors....
            ave_check2<-mean(vs_prop_factor0 %>% apply(1, sum))
            detach("package:dplyr")
            vs_prop_factor <- vs_prop_factor0 %>% apply(2,mean)
            ## Caltulate cutoff
            localcutoff<-apply(vs_permutate_factor, 2, quantile, probs=0.95)
            
            ## Merge cutoff and vs_prop_factor by var names 
            cut_mat<-data.frame("localcutoff"=localcutoff, "varnames"=names(localcutoff))
            prop_mat<-data.frame("prop"=vs_prop_factor, "varnames"=names(vs_prop_factor))
            vs_mat<-merge(cut_mat, prop_mat, by="varnames", all=TRUE)
            
            vs_mat$local_results<-ifelse(vs_mat$prop>vs_mat$localcutoff,1,0)
            table(vs_mat$local_results) # 11 selected
            vs_mat$varnames[vs_mat$local_results==1]
            
            # Global max
            vs_mat$globalmaxcut<-quantile(apply(vs_permutate_factor,1,max), probs=0.95)
            vs_mat$globalmax_results<-ifelse(vs_mat$prop>vs_mat$globalmaxcut,1,0)
            
            # Global SE
            mk<-apply(vs_permutate_factor,2,mean) # k is number of cov
            sk<-apply(vs_permutate_factor,2,sd)
            cover_constant<-bisectK(tol = 0.01, coverage = 1 - 0.05, 
                                    permute_mat = vs_permutate_factor, x_left = 1, x_right = 20, 
                                    countLimit = 100, perm_mean = mk, perm_se = sk)
            globalse_cut<-mk+cover_constant*sk
            cut_mat2<-data.frame("globalse_cut"=globalse_cut, "varnames"=names(globalse_cut))
            vs_mat2<-merge(cut_mat2, vs_mat, by="varnames", all=TRUE)
            vs_mat2$globalse_results<-ifelse(vs_mat2$prop>vs_mat2$globalse_cut,1,0)
            selected_local_var_i<-vs_mat2$varnames[vs_mat2$local_results==1]
            selected_globalmax_var_i<-vs_mat2$varnames[vs_mat2$globalmax_results==1]
            selected_globalse_var_i<-vs_mat2$varnames[vs_mat2$globalse_results==1]
            selected_var_b<-list("local_var_boot"=selected_local_var_i, 
                                 "globalmax_var_boot"=selected_globalmax_var_i, 
                                 "globalse_var_boot"=selected_globalse_var_i)
            
            list(selected_var_b)
          }
stopCluster(cl)   
timestamp()

# Metrics
cl <- makeCluster(25)
registerDoParallel(cl)
finalMatrix <- foreach(i=1:n_rep, .combine="comb", .multicombine = TRUE,  
                       .init=list(list(), list(), list(), list(), list())) %dopar% {
                         selected_var_i<-impute_boot[[1]][[i]]
                         
                         selected_local_var_i<-lapply(selected_var_i, function(x){x$local_var_boot})
                         selected_globalmax_var_i<-lapply(selected_var_i, function(x){x$globalmax_var_boot})
                         selected_globalse_var_i<-lapply(selected_var_i, function(x){x$globalse_var_boot})
                         
                         select_local_freq_i<-table(unlist(selected_local_var_i))
                         select_globalmax_freq_i<-table(unlist(selected_globalmax_var_i))
                         select_globalse_freq_i<-table(unlist(selected_globalse_var_i))
                         
                         final_local_var_i<-boot_finalvar(cutpt_boot_varselect=cutpt_boot_varselect, 
                                                          varlist=selected_local_var_i, n_boot=n_boot_impute)
                         final_globalmax_var_i<-boot_finalvar(cutpt_boot_varselect=cutpt_boot_varselect, 
                                                              varlist=selected_globalmax_var_i, n_boot=n_boot_impute)
                         final_globalse_var_i<-boot_finalvar(cutpt_boot_varselect=cutpt_boot_varselect, 
                                                             varlist=selected_globalse_var_i, n_boot=n_boot_impute)
                         final_var_i<-list("local_finalvar"=final_local_var_i,
                                           "globalmax_finalvar"=final_globalmax_var_i,
                                           "globalse_finalvar"=final_globalse_var_i)
                         
                         # 5. Metrics
                         prec_recall_f1_met_local<-list(NA)
                         prec_recall_f1_met_globalmax<-list(NA)
                         prec_recall_f1_met_globalse<-list(NA)
                         
                         metrics_mat_local<-data.frame(matrix(NA, nrow=length(cutpt_boot_varselect), ncol=4))
                         names(metrics_mat_local)<-c("Cutpt","Precision","Recall","F1")
                         metrics_mat_globalmax<-data.frame(matrix(NA, nrow=length(cutpt_boot_varselect), ncol=4))
                         names(metrics_mat_globalmax)<-c("Cutpt","Precision","Recall","F1")
                         metrics_mat_globalse<-data.frame(matrix(NA, nrow=length(cutpt_boot_varselect), ncol=4))
                         names(metrics_mat_globalse)<-c("Cutpt","Precision","Recall","F1")
                         
                         for(j in 1:length(cutpt_boot_varselect)){
                           #final_var_j<-final_var_i[[j]]
                           ## Precision, Recall, F1
                           if(length(final_var_i$local_finalvar)>0){
                             prec_recall_f1_met_local[[j]]<-prec_recall_f1(gold_pos=c("x1","x2","x3","x4","x5","x6","x7","x8"),
                                                                           gold_neg=c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                                                                      "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20"),
                                                                           var_selected=final_var_i$local_finalvar[[j]])
                             
                           }else{
                             prec_recall_f1_met_globalmax[[j]]<-list("tp"=NA, "fp"=NA, "fn"=NA, 
                                                                     "precision"=NA, "recall"=NA, "f1"=NA)
                           }
                           if(length(final_var_i$globalmax_finalvar)>0){
                             prec_recall_f1_met_globalmax[[j]]<-prec_recall_f1(gold_pos=c("x1","x2","x3","x4","x5","x6","x7","x8"),
                                                                               gold_neg=c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                                                                          "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20"),
                                                                               var_selected=final_var_i$globalmax_finalvar[[j]])
                           }else{
                             prec_recall_f1_met_globalmax[[j]]<-list("tp"=NA, "fp"=NA, "fn"=NA, 
                                                                     "precision"=NA, "recall"=NA, "f1"=NA)
                           }
                           if(length(final_var_i$globalse_finalvar)>0){
                             prec_recall_f1_met_globalse[[j]]<-prec_recall_f1(gold_pos=c("x1","x2","x3","x4","x5","x6","x7","x8"),
                                                                              gold_neg=c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                                                                         "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20"),
                                                                              var_selected=final_var_i$globalse_finalvar[[j]])
                             
                           }else{
                             prec_recall_f1_met_globalmax[[j]]<-list("tp"=NA, "fp"=NA, "fn"=NA, 
                                                                     "precision"=NA, "recall"=NA, "f1"=NA)
                           }
                           ## Precision, recall f1
                           metrics_mat_local[j,]<-c(cutpt_boot_varselect[j], prec_recall_f1_met_local[[j]]$precision, 
                                                    prec_recall_f1_met_local[[j]]$recall, prec_recall_f1_met_local[[j]]$f1)
                           metrics_mat_globalmax[j,]<-c(cutpt_boot_varselect[j], prec_recall_f1_met_globalmax[[j]]$precision, 
                                                        prec_recall_f1_met_globalmax[[j]]$recall, prec_recall_f1_met_globalmax[[j]]$f1)
                           metrics_mat_globalse[j,]<-c(cutpt_boot_varselect[j], prec_recall_f1_met_globalse[[j]]$precision, 
                                                       prec_recall_f1_met_globalse[[j]]$recall, prec_recall_f1_met_globalse[[j]]$f1)
                           
                           
                         }
                         prec_recall_f1_met<-list("local_met"=prec_recall_f1_met_local, 
                                                  "globalmax_met"=prec_recall_f1_met_globalmax, 
                                                  "globalse_met"=prec_recall_f1_met_globalse)
                         
                         list(final_var_i, prec_recall_f1_met, metrics_mat_local, metrics_mat_globalmax, metrics_mat_globalse)
                         
                       }

final_results<-list("dat"=dat_parallel, "impute_boot"=impute_boot, "met"=finalMatrix)
saveRDS(final_results,  paste(path,"x",seedx,"final_results_riaftbart.rds",sep=""))
final_results<-readRDS(paste(path,"x",seedx,"final_results_riaftbart.rds",sep=""))
finalMatrix<-final_results[[3]]

#stop cluster
stopCluster(cl)
Sys.time()

final_var0<-finalMatrix[[1]]
prec_recall_f1_i<-finalMatrix[[2]]
metrics_mat_local0<-finalMatrix[[3]]
metrics_mat_globalmax0<-finalMatrix[[4]]
metrics_mat_globalse0<-finalMatrix[[5]]



# Metrics
metrics_mat_local_list<-list(NA)
metrics_mat_globalmax_list<-list(NA)
metrics_mat_globalse_list<-list(NA)

for(i in 1:length(cutpt_boot_varselect)){
  # Precision, recall, f1
  metrics_mat_local_list_i0<- lapply(metrics_mat_local0, function(x){x[i,]})
  metrics_mat_local_list_i<-data.frame(matrix(unlist(metrics_mat_local_list_i0), nrow=n_rep, byrow=T))
  names(metrics_mat_local_list_i)<-c("Cutpt","Precision","Recall","F1")
  
  metrics_mat_globalmax_list_i0<- lapply(metrics_mat_globalmax0, function(x){x[i,]})
  metrics_mat_globalmax_list_i<-data.frame(matrix(unlist(metrics_mat_globalmax_list_i0), 
                                                  nrow=n_rep, byrow=T))
  names(metrics_mat_globalmax_list_i)<-c("Cutpt","Precision","Recall","F1")
  
  metrics_mat_globalse_list_i0<- lapply(metrics_mat_globalse0, function(x){x[i,]})
  metrics_mat_globalse_list_i<-data.frame(matrix(unlist(metrics_mat_globalse_list_i0), nrow=n_rep, 
                                                 byrow=T))
  names(metrics_mat_globalse_list_i)<-c("Cutpt","Precision","Recall","F1")
  
  ## Power & type I error
  final_var_local0<-lapply(final_var0, function(x){x$local_finalvar})
  final_var_local<-list(NA)       
  for(j in 1:n_rep){
    final_var_local_j<-final_var_local0[[j]]
    if(length(final_var_local_j)>0){
      final_var_local[[j]]<-final_var_local_j[[i]]
    }
  }
  final_var_local_factor<-factor(unlist(final_var_local), 
                                 levels = c("x1","x2","x3","x4","x5","x6","x7","x8",
                                            "z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                            "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20"))
  metrics_mat_local_list_i[,c("power_x1","power_x2","power_x3","power_x4","power_x5",
                              "power_x6","power_x7","power_x8",
                              "type1err_z1","type1err_z2","type1err_z3","type1err_z4","type1err_z5","type1err_z6","type1err_z7","type1err_z8","type1err_z9","type1err_z10",
                              "type1err_z11","type1err_z12","type1err_z13","type1err_z14","type1err_z15","type1err_z16","type1err_z17","type1err_z18","type1err_z19","type1err_z20")]<-
    matrix(rep(table(final_var_local_factor)/n_rep, n_rep), ncol=50, byrow = T)
  metrics_mat_local_list_i$power_overall<-apply(metrics_mat_local_list_i[,c("power_x1","power_x2","power_x3","power_x4","power_x5",
                                                                            "power_x6","power_x7","power_x8")],1,mean)
  metrics_mat_local_list_i$type1err_overall<-apply(metrics_mat_local_list_i[,c("type1err_z1","type1err_z2","type1err_z3","type1err_z4","type1err_z5","type1err_z6",
                                                                               "type1err_z7","type1err_z8","type1err_z9","type1err_z10",
                                                                               "type1err_z11","type1err_z12","type1err_z13","type1err_z14","type1err_z15","type1err_z16","type1err_z17","type1err_z18","type1err_z19","type1err_z20")],1,mean)
  
  ## Max
  final_var_globalmax0<-lapply(final_var0, function(x){x$globalmax_finalvar})
  final_var_globalmax<-list(NA)       
  for(j in 1:n_rep){
    final_var_globalmax_j<-final_var_globalmax0[[j]]
    if(length(final_var_globalmax_j)>0){
      final_var_globalmax[[j]]<-final_var_globalmax_j[[i]]
    }
  }
  
  final_var_globalmax_factor<-factor(unlist(final_var_globalmax), 
                                     levels = c("x1","x2","x3","x4","x5","x6","x7","x8",
                                                "z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                                "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20"))
  metrics_mat_globalmax_list_i[,c("power_x1","power_x2","power_x3","power_x4","power_x5",
                                  "power_x6","power_x7","power_x8",
                                  "type1err_z1","type1err_z2","type1err_z3","type1err_z4","type1err_z5","type1err_z6","type1err_z7","type1err_z8","type1err_z9","type1err_z10",
                                  "type1err_z11","type1err_z12","type1err_z13","type1err_z14","type1err_z15","type1err_z16","type1err_z17","type1err_z18","type1err_z19","type1err_z20")]<-
    matrix(rep(table(final_var_globalmax_factor)/n_rep, n_rep), ncol=50, byrow = T)
  metrics_mat_globalmax_list_i$power_overall<-apply(metrics_mat_globalmax_list_i[,c("power_x1","power_x2","power_x3","power_x4","power_x5",
                                                                                    "power_x6","power_x7","power_x8")],1,mean)
  metrics_mat_globalmax_list_i$type1err_overall<-apply(metrics_mat_globalmax_list_i[,c("type1err_z1","type1err_z2","type1err_z3","type1err_z4","type1err_z5","type1err_z6",
                                                                                       "type1err_z7","type1err_z8","type1err_z9","type1err_z10",
                                                                                       "type1err_z11","type1err_z12","type1err_z13","type1err_z14","type1err_z15","type1err_z16","type1err_z17","type1err_z18","type1err_z19","type1err_z20")],1,mean)
  
  ## SE
  final_var_globalse0<-lapply(final_var0, function(x){x$globalse_finalvar})
  final_var_globalse<-list(NA)       
  for(j in 1:n_rep){
    final_var_globalse_j<-final_var_globalse0[[j]]
    if(length(final_var_globalse_j)>0){
      final_var_globalse[[j]]<-final_var_globalse_j[[i]]
    }
  }
  
  final_var_globalse<-lapply(final_var0, function(x){x$globalse_finalvar[[i]]})
  final_var_globalse_factor<-factor(unlist(final_var_globalse), 
                                    levels = c("x1","x2","x3","x4","x5","x6","x7","x8",
                                               "z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                                               "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20"))
  metrics_mat_globalse_list_i[,c("power_x1","power_x2","power_x3","power_x4","power_x5",
                                 "power_x6","power_x7","power_x8",
                                 "type1err_z1","type1err_z2","type1err_z3","type1err_z4","type1err_z5","type1err_z6","type1err_z7","type1err_z8","type1err_z9","type1err_z10",
                                 "type1err_z11","type1err_z12","type1err_z13","type1err_z14","type1err_z15","type1err_z16","type1err_z17","type1err_z18","type1err_z19","type1err_z20")]<-
    matrix(rep(table(final_var_globalse_factor)/n_rep, n_rep), ncol=50, byrow = T)
  metrics_mat_globalse_list_i$power_overall<-apply(metrics_mat_globalse_list_i[,c("power_x1","power_x2","power_x3","power_x4","power_x5",
                                                                                  "power_x6","power_x7","power_x8")],1,mean)
  metrics_mat_globalse_list_i$type1err_overall<-apply(metrics_mat_globalse_list_i[,c("type1err_z1","type1err_z2","type1err_z3","type1err_z4","type1err_z5","type1err_z6",
                                                                                     "type1err_z7","type1err_z8","type1err_z9","type1err_z10",
                                                                                     "type1err_z11","type1err_z12","type1err_z13","type1err_z14","type1err_z15","type1err_z16","type1err_z17","type1err_z18","type1err_z19","type1err_z20")],1,mean)
  
  
  metrics_mat_local_list[[i]]<-metrics_mat_local_list_i
  metrics_mat_globalmax_list[[i]]<-metrics_mat_globalmax_list_i  
  metrics_mat_globalse_list[[i]]<-metrics_mat_globalse_list_i  
}

metrics_mat_list<-list("met_local"=metrics_mat_local_list, "met_globalmax"=metrics_mat_globalmax_list, "met_globalse"=metrics_mat_globalse_list)
# Export
saveRDS(metrics_mat_list, paste(path,"x",seedx, "metrics_mat_riaftbart.rds",sep=""))