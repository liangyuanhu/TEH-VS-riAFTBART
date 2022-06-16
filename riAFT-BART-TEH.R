#==================================================================================#
# Function to calcualte median survival using riAFT-BART                           #
#==================================================================================#

library(riAFTBART)
ri_AFTBART <- function(M.burnin, M.keep, status, y, x, trt, cluster.id, verbose = TRUE){
  # Counterfactual survival prob for treatment = 1
    riAFTBART_ATE_1 <- riAFTBART_fit(M.burnin = M.burnin,
                                                          M.keep = M.keep,
                                                          status = status,
                                                          y.train = y,
                                                          x.train = x,
                                                          trt.train = trt,
                                                          x.test = x,
                                                          trt.test = 1,
                                                          cluster.id = cluster.id,
                                                          verbose = verbose)
   
    riAFTBART_ATE_surv_prob_trt_1 <- cal_surv_prob(object = riAFTBART_ATE_1, time.points = seq(1,100,0.01), test.only = T,cluster.id = cluster.id)
    # Calculate the counterfactual 5 year survival probability and RMST for treatment = 1
    predict_trt_1_survival <- riAFTBART_ATE_surv_prob_trt_1$Surv.test[,riAFTBART_ATE_surv_prob_trt_1$time.points == 5 * 12] 
    predict_trt_1_rmst <- rep(NA, dim(x)[1])
    for (i in 1:dim(x)[1]){
      predict_trt_1_rmst[i] <- rmst(riAFTBART_ATE_surv_prob_trt_1$time.points, riAFTBART_ATE_surv_prob_trt_1$Surv.test[i,], max.time = 60)
    }
    
    # Counterfactual survival prob for treatment = 2
    riAFTBART_ATE_2 <- riAFTBART_fit(M.burnin = M.burnin,
                                          M.keep = M.keep,
                                          status = status,
                                          y.train = y,
                                          x.train = x,
                                          trt.train = trt,
                                          x.test = x,
                                          trt.test = 2,
                                          cluster.id = cluster.id,
                                          verbose = verbose)
    
    riAFTBART_ATE_surv_prob_trt_2 <- cal_surv_prob(object = riAFTBART_ATE_2, time.points = seq(1,80,0.01), test.only = T,cluster.id = cluster.id)
    # Calculate the counterfactual 5 year survival probability and RMST for treatment = 2
    predict_trt_2_survival <- riAFTBART_ATE_surv_prob_trt_2$Surv.test[,riAFTBART_ATE_surv_prob_trt_2$time.points == 5 * 12] 
    predict_trt_2_rmst <- rep(NA, dim(x)[1])
    for (i in 1:dim(x)[1]){
      predict_trt_2_rmst[i] <- rmst(riAFTBART_ATE_surv_prob_trt_2$time.points, riAFTBART_ATE_surv_prob_trt_2$Surv.test[i,], max.time = 60)
    }
    
    # Counterfactual survival prob for treatment = 3
    riAFTBART_ATE_3 <- riAFTBART_fit(M.burnin = M.burnin,
                                          M.keep = M.keep,
                                          status = status,
                                          y.train = y,
                                          x.train = x,
                                          trt.train = trt,
                                          x.test = x,
                                          trt.test = 3,
                                          cluster.id = cluster.id,
                                          verbose = verbose)
    
    riAFTBART_ATE_surv_prob_trt_3 <- cal_surv_prob(object = riAFTBART_ATE_3, time.points = seq(1,80,0.01), test.only = T,cluster.id = cluster.id)
    # Calculate the counterfactual 5 year survival probability and RMST for treatment = 3
    predict_trt_3_survival <- riAFTBART_ATE_surv_prob_trt_3$Surv.test[,riAFTBART_ATE_surv_prob_trt_3$time.points == 5 * 12] 
    predict_trt_3_rmst <- rep(NA, dim(x)[1])
    for (i in 1:dim(x)[1]){
      predict_trt_3_rmst[i] <- rmst(riAFTBART_ATE_surv_prob_trt_3$time.points, riAFTBART_ATE_surv_prob_trt_3$Surv.test[i,], max.time = 60)
    }
   
    pred_survival_rmst <- tibble(predict_trt_1_survival = predict_trt_1_survival,
                                   predict_trt_2_survival = predict_trt_2_survival,
                                   predict_trt_3_survival = predict_trt_3_survival,
                                   predict_trt_1_rmst = predict_trt_1_rmst,
                                   predict_trt_2_rmst = predict_trt_2_rmst,
                                   predict_trt_3_rmst = predict_trt_3_rmst)
    return(pred_survival_rmst)
}