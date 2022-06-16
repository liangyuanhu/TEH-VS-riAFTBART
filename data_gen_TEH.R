data_gen = function(nk=10, K=200, PH = TRUE, overlap = "moderate to strong", HS = "iii") {
  # nk=20; K=850; PH = TRUE
  n <- nk * K
  if (overlap == "moderate to strong"){
    phi = 1
  } else if (overlap == "weak"){
    phi = 2
  }
  n <- nk * K
  alpha1=1.5
  alpha2=0.7
  
  X = matrix(rnorm(5*n), nrow=n, ncol=5)
  x1 = X[,1]; x2 = X[,2]; x3 = X[,3]
  x4 = X[,4]; x5 = X[,5]
  x6 <- sample(size = n, c(0,1,2), prob = c(.3,.3,.4), replace = TRUE)
  x7 <- sample(size = n, c(0,1,2), prob = c(.3,.3,.4), replace = TRUE)
  X <- cbind(X, x6, x7)
  # random intercepts for Amodel
  xi = rnorm(K, 0,1)
  # cluster label
  cl = rep(1:K, each=nk)
  xiX  = rep(xi, each=nk)
  # A-Model, nonlinear treatment assignment
  ex1 = exp(alpha1 + phi*(0.1*x1 + 0.1*x2+ 0.1*x3+ 0.5*x4+0.4*x5+0.2*x6+0.3*x7+0.4*x2^2+0.4*x2^2*x5+xiX))
  ex2 = exp(alpha2 + phi*(0.1*x1 + 0.3*x2+ 0.2*x3+ 0.2*x4+0.1*x5+0.4*x6+0.5*x7-0.3*x2*x4+0.7*x2^2*x4+xiX))
  Ap1 = ex1 / (1 + ex1 + ex2)
  Ap2 = ex2 / (1 + ex1 + ex2)
  Ap3 = 1 - Ap1 - Ap2
  
  Ap1 = ifelse(Ap1 < 10^(-10), 0, Ap1)
  Ap2 = ifelse(Ap2 < 10^(-10), 0, Ap2)
  Ap3 = ifelse(Ap3 < 10^(-10), 0, Ap3)
  
  
  A = NULL
  for (i in 1:n) A[i] = sample(c(0,2,1), size=1, replace=TRUE, prob=c(Ap1[i],Ap2[i],Ap3[i]))
  table(A)
  if (PH==FALSE) {eta = exp(0.7+ 0.5*x1)} #shape parameter>0}
  if (PH==TRUE) {eta = 2}
  #partially Aligned Non-parallel T-Model with nPH  
  # cluster-specific random intercepts
  bk = rnorm (K,0,4)
  bkT = rep(bk, each = nk)
  #generate U ~ unif(0,1)
  U = runif(n,0,1)
  #scale parameter lambda>0, exp(Linear predictor)
  lambda1 = 5000
  lambda2 = 800
  lambda3 = 1200
  #linear predictors 
  if (HS == "iii") {
    LP1 = 1 * x1 + 0.3*x2+ sin(pi*x3)+ 0.6*x4+0.5*x5+1.2*x6 +0.3*x2^2+0.5*x4*x5 +bkT-1
    
    LP2 = 0.4*x1 + 1.2*sin(pi*x3)+0.4*x4+0.3*x5+1.0*x6+0.8*x7+0.7*x1^2+0.4*x1*x4+bkT
    
    LP3 = 0.4*sin(pi*x2) + 0.9*x3+ 0.9*x4+0.4*x5+0.4*x6+0.9*x7+0.4*x4^2-0.3*x2*x3 +bkT-3
  } else if (HS == "ii"){
    LP1 = 1 * x1 + 0.3*x2+ sin(pi*x3)+ 0.6*x4+0.5*x5+1.2*x6+ 0.3 * x7 +0.3*x2^2+0.5*x4*x5 +bkT -1
    
    LP2 = 0.4*x1 + 1.2*sin(pi*x2)+0.4*x3+0.3*x4+1.0*x5+0.8*x6 + 0.1 * x7 +0.7*x1^2+0.4*x1*x4+bkT
    
    LP3 = 0.4*sin(pi*x1) + 0.9*x2+ 0.9*x3+0.4*x4+0.4*x5+0.9*x6 + 0.3 * x7 +0.4*x4^2-0.3*x2*x3 +bkT -3
    
  } else if (HS == "i"){
    LP1 = 1 * x1 + 0.3*x2+ sin(pi*x3)+ 0.6*x4+0.5*x5+1.2*x6+0.4 *x7 +0.3*x2^2+0.5*x4*x5 +bkT -1
    
    LP2 = 0.4*x1 + 1.2*sin(pi*x2)+0.4*x3+0.3*x4+1.0*x5+0.8*x6+0.2 *x7+0.7*x1^2+0.4*x1*x4+bkT 
    
    LP3 = 0.4*x1 + 0.9*x2 + 0.4*x3 + 0.9*x4+0.4*x5+0.4*x6 +0.3 *x7+bkT-2
  }
  
  
  #potential survival times 
  T1 = (lambda1*(-log(U))/exp(LP1))^(1/eta)
  summary(T1)
  T2 = (lambda2*(-log(U))/exp(LP2))^(1/eta)
  summary(T2)
  T3 = (lambda3*(-log(U))/exp(LP3))^(1/eta)
  summary(T3)
  # observed outcomes
  T = cbind(T1,T2,T3)
  
  TA = cbind(T,A)
  Tobs = apply(TA, 1, function(x) x[1:3][x[4]+1]) #observed when trt is received 
  
  
  C <- rexp(n, rate = 0.02)
  
  #C = rexp(n,rate=0.00025) #for PH model
  Tobs_C = pmin(Tobs,C)
  summary(Tobs_C)
  #censoring rate
  sum(Tobs>C)/n 
  #censoring indicator
  delta = as.numeric(Tobs>C)
  return(list(T=T, Tobs=Tobs,Tobs_C=Tobs_C,delta=delta, A=A+1, X=X, nk=nk, K=K,cl=cl, C = C, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, LP1 = LP1, LP2 = LP2, LP3 = LP3, eta = eta, Ap1 = Ap1, Ap2 = Ap2, Ap3 = Ap3))
}
