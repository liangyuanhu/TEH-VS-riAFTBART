#==============================================================================#
# Function to calcualte ATE in terms of median survival using PEAMM           #
#==============================================================================#

library(pammtools)
library(mgcv)
PEAMM <- function(y, x, trt,status, cluster.id){
  mydata_df <- tibble(Time = y, Event= status, cluster.id=cluster.id,treatment = trt) %>% 
    bind_cols(as.data.frame(x)) %>% 
    mutate(trt = as.factor(trt)) 
  # Transform the data into the format for ri-GAPH
  mydata_df_ped <- mydata_df %>% as_ped(Surv(Time, Event)~ V1+V2+V3+V4+V5+x6+x7+x8+x9+x10 + trt + cluster.id, id = "id")
  # Fit the ri-GAPH model
  riGAPH_mydata <- gam(ped_status ~ s(tend) + V1+V2+V3+V4+V5+x6+x7+x8+x9+x10 + trt, data = mydata_df_ped,family = poisson(), offset = offset)
  # Counterfactual survival curve for treatment = 1
  riGAPH_trt_1 <- mydata_df_ped  %>% mutate(intlen = tend - tstart, trt = 1)
  riGAPH_trt_1_surv <- riGAPH_trt_1 %>% add_surv_prob(riGAPH_mydata)
  predict_trt_1_survival <- riGAPH_trt_1_surv %>% 
    select(id, tstart, surv_prob) %>% 
    unique %>% 
    as_tibble() %>% 
    mutate(tstart_close = abs(tstart - 5 *12)) %>% 
    group_by(id) %>% 
    filter(tstart_close == min(tstart_close)) %>% 
    ungroup %>% 
    arrange(id) %>% 
    pull(surv_prob)
  predict_trt_1_rmst <- rep(NA, dim(sim_data$X)[1])
  for (i in 1:dim(sim_data$X)[1]){
    riGAPH_trt_1_surv_once <- riGAPH_trt_1_surv %>% 
      select(id, tstart, surv_prob) %>% 
      unique %>% 
      as_tibble() %>% 
      filter(id == i)
    predict_trt_1_rmst[i] <- rmst(riGAPH_trt_1_surv_once$tstart, riGAPH_trt_1_surv_once$surv_prob,
         max.time = 5*12)
  }
    
  # Counterfactual survival curve for treatment = 2
  riGAPH_trt_2 <- mydata_df_ped  %>% mutate(intlen = tend - tstart, trt = 2)
  riGAPH_trt_2_surv <- riGAPH_trt_2 %>% add_surv_prob(riGAPH_mydata)
  predict_trt_2_survival <- riGAPH_trt_2_surv %>% 
    select(id, tstart, surv_prob) %>% 
    unique %>% 
    as_tibble() %>% 
    mutate(tstart_close = abs(tstart - 5 *12)) %>% 
    group_by(id) %>% 
    filter(tstart_close == min(tstart_close)) %>% 
    ungroup %>% 
    arrange(id) %>% 
    pull(surv_prob)
  predict_trt_2_rmst <- rep(NA, dim(sim_data$X)[1])
  for (i in 1:dim(sim_data$X)[1]){
    riGAPH_trt_2_surv_once <- riGAPH_trt_2_surv %>% 
      select(id, tstart, surv_prob) %>% 
      unique %>% 
      as_tibble() %>% 
      filter(id == i)
    predict_trt_2_rmst[i] <- rmst(riGAPH_trt_2_surv_once$tstart, riGAPH_trt_2_surv_once$surv_prob,
                                  max.time = 5*12)
  }
  
  # Counterfactual survival curve for treatment = 3
  riGAPH_trt_3 <- mydata_df_ped  %>% mutate(intlen = tend - tstart, trt = 3)
  riGAPH_trt_3_surv <- riGAPH_trt_3 %>% add_surv_prob(riGAPH_mydata)
  predict_trt_3_survival <- riGAPH_trt_3_surv %>% 
    select(id, tstart, surv_prob) %>% 
    unique %>% 
    as_tibble() %>% 
    mutate(tstart_close = abs(tstart - 5 *12)) %>% 
    group_by(id) %>% 
    filter(tstart_close == min(tstart_close)) %>% 
    ungroup %>% 
    arrange(id) %>% 
    pull(surv_prob)
  predict_trt_3_rmst <- rep(NA, dim(sim_data$X)[1])
  for (i in 1:dim(sim_data$X)[1]){
    riGAPH_trt_3_surv_once <- riGAPH_trt_3_surv %>% 
      select(id, tstart, surv_prob) %>% 
      unique %>% 
      as_tibble() %>% 
      filter(id == i)
    predict_trt_3_rmst[i] <- rmst(riGAPH_trt_3_surv_once$tstart, riGAPH_trt_3_surv_once$surv_prob,
                                  max.time = 5*12)
  }
  
  pred_survival_rmst <- tibble(predict_trt_1_survival = predict_trt_1_survival,
                               predict_trt_2_survival = predict_trt_2_survival,
                               predict_trt_3_survival = predict_trt_3_survival,
                               predict_trt_1_rmst = predict_trt_1_rmst,
                               predict_trt_2_rmst = predict_trt_2_rmst,
                               predict_trt_3_rmst = predict_trt_3_rmst)
  return(pred_survival_rmst)
  
}