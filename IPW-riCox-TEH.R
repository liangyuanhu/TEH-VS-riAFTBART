#=================================================================================#
# Function to calcualte ATE in terms of median survival using IPW-riCox           #
#                                                                                 #
#=================================================================================#

library(timereg)
IPW_ri_Cox <- function(y, x, trt,status, cluster.id){
  mydata_df <- tibble(Time = y, Event= status, cluster.id=cluster.id,treatment = trt) %>% 
    bind_cols(as.data.frame(x)) %>% 
    mutate(trt = as.factor(trt),
           cluster.id = as.factor(cluster.id)) 
  # IPTW with weight estimated by superleaner 
  weight_sl <- WeightIt::weightit(trt ~ V1+V2+V3+V4+V5+x6+x7+x8+x9+x10 + cluster.id, data = mydata_df, method = "super", SL.library = c("SL.xgboost", "SL.bartMachine","SL.gbm"))
  cox_model <- coxph(Surv(Time, Event) ~ ., data = mydata_df, weight = weight_sl$weights)
  # Calculate the predicted 5 year survival and RMST for treatment = 1
  survfit_cox_trt_1 <- survfit(cox_model, newdata = mydata_df %>% mutate(trt = 1) %>% mutate(trt = as.factor(trt)))
  predict_trt_1_survival <- rep(NA, dim(mydata_df)[1])
  predict_trt_1_rmst <- rep(NA, dim(mydata_df)[1])
  for (i in 1:dim(mydata_df)[1]){
    predict_trt_1_survival[i] <- survfit_cox_trt_1$surv[,i][which.min(abs(survfit_cox_trt_1$time - 5 *12))]
    predict_trt_1_rmst[i] <- rmst(survfit_cox_trt_1$time, survfit_cox_trt_1$surv[,i], max.time = 60)
  }
  
  # Calculate the predicted 5 year survival and RMST for treatment = 2
  survfit_cox_trt_2 <- survfit(cox_model, newdata = mydata_df %>% mutate(trt = 2) %>% mutate(trt = as.factor(trt)))
  predict_trt_2_survival <- rep(NA, dim(mydata_df)[1])
  predict_trt_2_rmst <- rep(NA, dim(mydata_df)[1])
  for (i in 1:dim(mydata_df)[1]){
    predict_trt_2_survival[i] <- survfit_cox_trt_2$surv[,i][which.min(abs(survfit_cox_trt_2$time - 5 *12))]
    predict_trt_2_rmst[i] <- rmst(survfit_cox_trt_2$time, survfit_cox_trt_2$surv[,i], max.time = 60)
  }
  
  # Calculate the predicted 5 year survival and RMST for treatment = 3
  survfit_cox_trt_3 <- survfit(cox_model, newdata = mydata_df %>% mutate(trt = 3) %>% mutate(trt = as.factor(trt)))
  predict_trt_3_survival <- rep(NA, dim(mydata_df)[1])
  predict_trt_3_rmst <- rep(NA, dim(mydata_df)[1])
  for (i in 1:dim(mydata_df)[1]){
    predict_trt_3_survival[i] <- survfit_cox_trt_3$surv[,i][which.min(abs(survfit_cox_trt_3$time - 5 *12))]
    predict_trt_3_rmst[i] <- rmst(survfit_cox_trt_3$time, survfit_cox_trt_3$surv[,i], max.time = 60)
  }
  
  pred_survival_rmst <- tibble(predict_trt_1_survival = predict_trt_1_survival,
                               predict_trt_2_survival = predict_trt_2_survival,
                               predict_trt_3_survival = predict_trt_3_survival,
                               predict_trt_1_rmst = predict_trt_1_rmst,
                               predict_trt_2_rmst = predict_trt_2_rmst,
                               predict_trt_3_rmst = predict_trt_3_rmst)
  return(pred_survival_rmst)
}