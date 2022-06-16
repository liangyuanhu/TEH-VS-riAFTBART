dat_gen <- function(n = 2000, K = 10, nK = 200, PH=TRUE){
  x1<-rbinom(n,1,0.5)
  x2<-rbinom(n,1,0.5)
  x3<-rnorm(n,0,1)
  x4<-rnorm(n,0,1)
  x5<-rnorm(n, 0.3*x2+0.2*x3,1)
  x6<-rnorm(n,-0.4*x3+0.4*x4+0.3*x3*x4,1)
  x7<-rnorm(n, 0.1*x4*(x5-2)^2-0.1*x6^2, 1)
  x8<-rnorm(n, 0.5*x6+0.3*x7-0.3*x5^2+0.2*x6*x7, 1)
  
  # Z: 10 N(0,1) and 10 bin(1,0.5)
  z1<-rnorm(n,0,1); z2<-rnorm(n,0,1); z3<-rnorm(n,0,1); z4<-rnorm(n,0,1); z5<-rnorm(n,0,1)
  z6<-rnorm(n,0,1); z7<-rnorm(n,0,1); z8<-rnorm(n,0,1); z9<-rnorm(n,0,1); z10<-rnorm(n,0,1)
  
  z11<-rbinom(n,1,0.5); z12<-rbinom(n,1,0.5); z13<-rbinom(n,1,0.5); z14<-rbinom(n,1,0.5); z15<-rbinom(n,1,0.5)
  z16<-rbinom(n,1,0.5); z17<-rbinom(n,1,0.5); z18<-rbinom(n,1,0.5); z19<-rbinom(n,1,0.5); z20<-rbinom(n,1,0.5)
  
  
  if (PH==FALSE) {eta = exp(0.7+ 0.5*x1)} #shape parameter>0}
  if (PH==TRUE) {eta = 2}
  bk = rnorm (K,0,4)
  cl = rep(1:K, each=nK)
  bkT = rep(bk, each = nK)
  #generate U ~ unif(0,1)
  U = runif(n,0,1)
  LP = 1.8*x1+0.5*x2+1.1*x3-0.4*exp(x5)+0.4*(x6-1.5)^2+0.1*(x7-0.1)^3-5*sin(0.1*pi*x4*x8)-0.4*x5*x7
  summary(LP)
  time = (3000*(-log(U))/exp(LP/3))^(1/eta)
  summary(time)
  time_censor <- rexp(n, rate = 0.02)
  summary(time_censor)
  time_observed = pmin(time,time_censor)
  #censoring rate
  sum(time>time_censor)/n 
  
  #censoring indicator
  delta = as.numeric(time>time_censor)
  table(delta)
  dat_comp<-data.frame(time_observed,delta,cl,x1,x2,x3,x4,x5,x6,x7,x8,
                       z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,
                       z11,z12,z13,z14,z15,z16,z17,z18,z19,z20) 
  dim(dat_comp)
  # 3. Ampute
  # Create interaction and polynomial for amputation
  dat_comp$x3x4<-dat_comp$x3*dat_comp$x4 # used in x5 x6 and x6 * x8
  dat_comp$x5x5<-dat_comp$x5*dat_comp$x5 # used in x6
  dat_comp$x4x5<-dat_comp$x4*dat_comp$x5 # used in x7 and x7 * x8
  dat_comp$x6x6<-dat_comp$x6*dat_comp$x6 # used in x7
  dat_comp$x6x7<-dat_comp$x6*dat_comp$x7 # used in x8
  
  na_pattern<-matrix(c( 
    c(1, 1, 1, rep(1,4), 0, 1, 1, 1, rep(1,20), rep(1,dim(dat_comp)[2]-31)), # x5
    c(1, 1, 1, rep(1,4), 1, 0, 1, 1, rep(1,20), rep(1,dim(dat_comp)[2]-31)), # x6
    c(1, 1, 1, rep(1,4), 1, 1, 0, 1, rep(1,20), rep(1,dim(dat_comp)[2]-31)), # x7
    c(1, 1, 1, rep(1,4), 1, 1, 1, 0, rep(1,20), rep(1,dim(dat_comp)[2]-31)), # x8
    c(1, 1, 1, rep(1,4), 0, 0, 1, 1, rep(1,20), rep(1,dim(dat_comp)[2]-31)), # x5 x6
    c(1, 1, 1, rep(1,4), 1, 0, 0, 1, rep(1,20), rep(1,dim(dat_comp)[2]-31)), # x6 x7
    c(1, 1, 1, rep(1,4), 1, 1, 0, 0, rep(1,20), rep(1,dim(dat_comp)[2]-31)), # x7 x8
    c(1, 1, 1, rep(1,4), 1, 0, 1, 0, rep(1,20), rep(1,dim(dat_comp)[2]-31)) # x6 x8
  ), ncol=dim(dat_comp)[2], byrow=T)
  
  na_wt<-matrix(c(
    c(0,0,0, 0,0,1,1,0,0, 0, 0, rep(0,20), 1, rep(0,4)), # x5
    c(0,0,0, 0,0,1,1,1,0, 0, 0, rep(0,20), 1, 1, rep(0,3)), # x6
    c(0,0,0, 0,0,0,1,1,0, 0, 0, rep(0,20), 0, 0, 1, 1, 0), # x7
    c(0,0,0, 0,0,0,0,1,1, 1, 0, rep(0,20), rep(0,4),1), # x8
    c(0,0,0, 0,0,1,1,0,0, 0, 0, rep(0,20), rep(0,5)), # x5 x6
    c(0,0,0, 0,0,0,0,1,0, 0, 0, rep(0,20), rep(0,5)), # x6 x7
    c(0,0,0, 0,0,0,1,1,0, 0, 0, rep(0,20), 0,0,0.5,0,0), # x7 x8
    c(0,0,0, 0,0,1,1,0,0, 0, 0, rep(0,20), 1,rep(0,4)) # x6 x8
  ), ncol=dim(dat_comp)[2], byrow=T)
  colnames(na_pattern)<-colnames(na_wt)<-names(dat_comp)
  dim(na_wt)
  dim(na_pattern)
  dat_na0<-ampute(dat_comp, prop=0.40, mech="MAR", patterns = na_pattern, 
                  freq = c(0.30,0.09,0.09,0.08,0.08,0.16,0.10,0.10), weights=na_wt,
                  cont = TRUE, 
                  type=c("RIGHT","RIGHT", "RIGHT","RIGHT", "RIGHT","TAIL", "TAIL","TAIL"))$amp
  
  dat_na<-dat_na0[,c("time_observed","delta", "cl","x1","x2","x3","x4","x5","x6","x7","x8",
                     "z1","z2","z3","z4","z5","z6","z7","z8","z9","z10",
                     "z11","z12","z13","z14","z15","z16","z17","z18","z19","z20")]
  return(list(dat_na = dat_na, 
         dat_comp = dat_comp))
}
