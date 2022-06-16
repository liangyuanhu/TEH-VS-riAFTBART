library(foreach)
library(doParallel)
# library(doSNOW) 
library(Metrics)
library(survival)
library(tidyverse)
source("code/data_gen_var_select.R")
source("code/stepcoxme.R")
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

stepcoxme <- function(data, event,time,cl, direction, alpha_in=0.05, alpha_out=0.05){
  # Exclude variables not use in the selection
  varname <- names(mice_b_dat)[-which(names(mice_b_dat)%in%c("time","event","cl"))]
  
  
  
  if(direction%in%c("for_back","forward")){
    FORWARD=TRUE
    BACKWARD=FALSE
    var_in_mod<-1
  }
  if(direction%in%c("back_for","backward")){
    FORWARD=FALSE
    BACKWARD=TRUE
    var_in_mod<-varname
  }
  # Output matrix
  outmat<-data.frame(matrix(NA, ncol = 5, nrow=1))
  names(outmat)<-c("nstep","Var_in_model","VarIn","VarOut", "p")
  nstep<-1
  var_in_mod0<-var_in_mod
  while(TRUE){
    # Compare var_in_model before and after last step. If same, then break
    if(nstep!=1 & 
       length(setdiff(var_in_mod0, var_in_mod))==0 &
       length(setdiff(var_in_mod, var_in_mod0))==0){
      break
    }else{
      var_in_mod0<-var_in_mod
    }
    
    if(FORWARD==TRUE){
      var_list<-varname[varname%in%var_in_mod==0]
      var_p<-data.frame(matrix(NA,ncol=2,nrow=length(var_list)))
      names(var_p)<-c("VarName","p")
      
      # Fit modl w/ each of addition var and select the most sig
      for(i in 1:length(var_list)){
        var_i<-var_list[i]
        fm_i <- paste("Surv(time,event)~", paste0(var_in_mod, collapse = "+"), "+",var_i,"+ (1|cl)",  sep = "")
        coxme_fit <- coxme(as.formula(fm_i), data = mice_b_dat)
        
        anova_coxme <-car::Anova(coxme_fit, test="Wald", type="III")
        var_p$VarName[i]<-rownames(anova_coxme)[length(var_in_mod)+1]
        var_p$p[i]<-anova_coxme$`Pr(>Chisq)`[length(var_in_mod)+1]
      } # end for loop
      p_in<-var_p$p[var_p$p==min(var_p$p)]
      var_in<-var_p$VarName[var_p$p==min(var_p$p)&var_p$VarName!="(Intercept)"]
      if(p_in<alpha_in){
        out_p<-ifelse(p_in<0.001, "<0.001",
                      as.character(round(p_in,3)))
        var_in_mod<-c(var_in_mod, var_in)
        outvec<-list(nstep, paste0(var_in_mod,collapse="+"), var_in, NA, out_p)
        outmat<-rbind(outmat, unlist(outvec))
        
        FORWARD=TRUE
        
        if(length(var_in_mod)==(ncol(mice_b_dat)-1)){
          if(direction=="backward"){
            break
          }else if(direction%in%c("for_back","back_for")){
            # Change to backward
            FORWARD=TRUE
            BACKWARD=FALSE
          }}
      }else {
        if(direction=="forward"){
          break
        }else if(direction%in%c("for_back","back_for")){
          # Change to backward
          FORWARD=FALSE
          BACKWARD=TRUE
        }
      }
    } # end if(FORWARD=T)
    if(BACKWARD==TRUE){
      var_list<-var_in_mod
      var_p<-data.frame(matrix(NA,ncol=2,nrow=length(var_list)))
      names(var_p)<-c("VarName","p")
      fm_i <- paste("Surv(time,event)~", paste0(var_list, collapse = "+"),"+ (1|cl)",  sep = "")
      coxme_fit <- coxme(as.formula(fm_i), data = mice_b_dat)
      anova_coxme <- car::Anova(coxme_fit, test="Wald", type="III")
      if(min(anova_coxme$`Pr(>Chisq)`[-1])==1){ # if all se inflated but still converge, then break
        break
      }
      p_out<-anova_coxme$`Pr(>Chisq)`[anova_coxme$`Pr(>Chisq)`==max(anova_coxme$`Pr(>Chisq)`)]
      var_out<-rownames(anova_coxme)[anova_coxme$`Pr(>Chisq)`==max(anova_coxme$`Pr(>Chisq)`)]
      
      if(p_out>alpha_out){
        out_p<-ifelse(p_out<0.001, "<0.001",
                      as.character(round(p_out,3)))
        var_in_mod<-var_in_mod[-which(var_in_mod==var_out)]
        
        outvec<-list(nstep, paste0(var_in_mod,collapse="+"), NA, var_out, out_p)
        outmat<-rbind(outmat, unlist(outvec))
        
        if(length(var_in_mod)==0){
          if(direction=="backward"){
            break
          }else if(direction%in%c("for_back","back_for")){
            # Change to backward
            FORWARD=TRUE
            BACKWARD=FALSE
          }
        }
        
      }else {
        if(direction=="backward"){
          break
        }else if(direction%in%c("for_back","back_for")){
          # Change to backward
          FORWARD=TRUE
          BACKWARD=FALSE
        }
      }
      
    } # end if(BACKWARD=T)
    
    nstep<-nstep+1
    
  } # end while
  outmat<-outmat[-1,]
  return(outmat)
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
            library(coxme)
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
            
            # Model
            m1<- try(stepcoxme(data = mice_b_dat, event = "event",time = "time",cl = "cl", direction = directionx))
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

saveRDS(final_results,  paste("result/final_results_cox.rds",sep=""))


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
saveRDS(metrics_mat_list, paste("result/metrics_mat_cox_boot.rds",sep=""))
