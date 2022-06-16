#===============================================================================#
# Function to calcualte median survival using DR-riAH #                         #                  
#===============================================================================#
# y = sim_data$Tobs_C
# x = sim_data$X
# trt = sim_data$A
# cluster.id = sim_data$cl
# status = sim_data$delta

library(timereg)
library(RISCA)
DR_riAH <- function(y, x, trt,status, cluster.id){
  mydata_df <- tibble(Time = y, Event= status, cluster.id=cluster.id,treatment = trt) %>% 
    bind_cols(as.data.frame(x)) %>% 
    mutate(trt = as.factor(trt),
           cluster.id = as.factor(cluster.id))
  # Fit the treatment model to estimate ps using Superlearner
  weight_sl <- WeightIt::weightit(trt ~ V1+V2+V3+V4+V5+x6+x7+x8+x9+x10 + cluster.id, data = mydata_df, method = "super", SL.library = c("SL.xgboost", "SL.bartMachine","SL.gbm"))
  ps_sl <- 1/weight_sl$weights
  
  # Calculate the predicted 5 year surviavl probabilty for treatment = 1
  ah_model_trt_1 <- aalen(Surv(Time, Event)~ V1+V2+V3+V4+V5+x6+x7+x8+x9+x10, data = mydata_df %>% as.data.frame() %>% filter(trt == 1), clusters = cluster.id, resample.iid=1)
  predict_ah_trt_1 <- predict(ah_model_trt_1, newdata = mydata_df %>% as.data.frame())
  predict_ah_trt_1_survival <- predict_ah_trt_1$S0[,which.min(abs(predict_ah_trt_1$time - 5 * 12))]
  
  # Calculate the predicted 5 year surviavl probabilty for treatment = 2
  ah_model_trt_2 <- aalen(Surv(Time, Event)~ V1+V2+V3+V4+V5+x6+x7+x8+x9+x10, data = mydata_df %>% as.data.frame() %>% filter(trt == 2), clusters = cluster.id, resample.iid=1)
  predict_ah_trt_2 <- predict(ah_model_trt_2, newdata = mydata_df %>% as.data.frame())
  predict_ah_trt_2_survival <- predict_ah_trt_2$S0[,which.min(abs(predict_ah_trt_2$time - 5 * 12))]
  
  # Calculate the predicted 5 year surviavl probabilty for treatment = 3
  ah_model_trt_3 <- aalen(Surv(Time, Event)~ V1+V2+V3+V4+V5+x6+x7+x8+x9+x10, data = mydata_df %>% as.data.frame() %>% filter(trt == 3), clusters = cluster.id, resample.iid=1)
  predict_ah_trt_3 <- predict(ah_model_trt_3, newdata = mydata_df %>% as.data.frame())
  predict_ah_trt_3_survival <- predict_ah_trt_3$S0[,which.min(abs(predict_ah_trt_3$time - 5 * 12))]
  # Calculate the predicted 5 year RMST for treatment = 1, 2, 3 respectively
  predict_ah_trt_3_rmst <- predict_ah_trt_2_rmst <- predict_ah_trt_1_rmst <- rep(NA, nrow = dim(mydata_df)[1])
  for (i in 1:dim(mydata_df)[1]){
    predict_ah_trt_1_rmst[i] <- rmst(predict_ah_trt_1$time,
                                     predict_ah_trt_1$S0[i,],
                                     max.time = 60)
    predict_ah_trt_2_rmst[i] <- rmst(predict_ah_trt_2$time,
                                     predict_ah_trt_2$S0[i,],
                                     max.time = 60)
    predict_ah_trt_3_rmst[i] <- rmst(predict_ah_trt_3$time,
                                     predict_ah_trt_3$S0[i,],
                                     max.time = 60)
    # print(i)
  }
  
  pred_survival_rmst <- tibble(predict_trt_1_survival = predict_ah_trt_1_survival,
                               predict_trt_2_survival = predict_ah_trt_2_survival,
                               predict_trt_3_survival = predict_ah_trt_3_survival,
                               predict_trt_1_rmst = predict_ah_trt_1_rmst,
                               predict_trt_2_rmst = predict_ah_trt_2_rmst,
                               predict_trt_3_rmst = predict_ah_trt_3_rmst)
  return(pred_survival_rmst)
}