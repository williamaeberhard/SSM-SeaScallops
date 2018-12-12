SS2v4simul <- function(gt,gRt,Nt,
                       B1,R1,m1,C1,
                       sigmatau,sigmaphi,sigmam,
                       sigmaeps,sigmaups,sigmakap,sigmaC,
                       qI,qR,S,a,chi,
                       seedvalue=NULL,warn=TRUE,tol.depletion=1e-6){
  # Simulate data according to the alterative SSM introduced in:
  #   Yin, Y., Aeberhard, W. H., Smith, S. J., and Mills Flemming, J. (2018).
  #   Identifiable state-space models: A case study of the Bay of Fundy sea
  #   scallop fishery. Canadian Journal of Statistics, In press.
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Setup ####
  #/////////////////////////////////////////////////////////////////////////////
  
  NY <- length(gt) # time series length, number of years
  if (length(gRt)!=NY || length(Nt)!=NY){
    stop('gt, gRt and Nt must be vectors of the same length.')
  }
  Bt <- double(NY) # biomass commercial size, randeff
  Rt <- Bt # number of recruits, randeff
  mt <- Bt # natural mortality rate, randeff
  It <- Bt # survey index biomass commercial size, response
  Jt <- Bt # survey index number of recruits, response
  Lt <- Bt # survey estimate of clappers, response
  Ct <- Bt # survey estimate of clappers, response
  if (!is.null(seedvalue)){set.seed(seedvalue)}
  
  #/////////////////////////////////////////////////////////////////////////////
  ##### Generate iid error, randeff and response ####
  #/////////////////////////////////////////////////////////////////////////////
  
  tau <- rlnorm(n=NY,meanlog=-sigmatau^2/2,sdlog=sigmatau) # proc error Bt
  phi <- rlnorm(n=NY,meanlog=-sigmaphi^2/2,sdlog=sigmaphi) # proc error Rt
  mmm <- rlnorm(n=NY,meanlog=-sigmam^2/2,sdlog=sigmam) # proc error mt
  
  eps <- rlnorm(n=NY,meanlog=-sigmaeps^2/2,sdlog=sigmaeps) # obs error It
  ups <- rlnorm(n=NY,meanlog=-sigmaups^2/2,sdlog=sigmaups) # obs error Jt
  kap <- rlnorm(n=NY,meanlog=-sigmakap^2/2,sdlog=sigmakap) # obs error Lt
  ccc <- rlnorm(n=NY,meanlog=-sigmaC^2/2,sdlog=sigmaC) # obs error Ct
  
  ### t==1: fixed randeff (B1, R1 and m1), fixed C1
  mt[1] <- m1 # no error for t=1, fixed
  Rt[1] <- R1 # no error for t=1, fixed
  Bt[1] <- B1 # no error for t=1, fixed
  mean.It <- qI*Bt[1]
  It[1] <- mean.It*eps[1] # lognormal error, E=1
  mean.Jt <- qR*Rt[1]
  Jt[1] <- mean.Jt*ups[1] # lognormal error, E=1
  mean.Lt <- mt[1]*S*Nt[1]
  Lt[1] <- mean.Lt*kap[1] # lognormal error, E=1
  # mean.Ct <- Bt[1]/10 # arbitrary...
  mean.Ct <- C1 # used to default to B1/10
  Ct[1]  <- mean.Ct # deterministic
  ### 2<=t<=NY
  t <- 2
  stockdepleted <- F
  while (t<=NY && !stockdepleted){
    # proc
    Rt[t] <- Rt[t-1]*phi[t] # log-RW, lognormal error, E=1
    mt[t] <- mt[t-1]*mmm[t] # log-RW, lognormal error, E=1
    mean.Bt <- exp(-mt[t])*gt[t-1]*(Bt[t-1]-Ct[t-1])+exp(-mt[t])*gRt[t-1]*Rt[t-1]
    Bt[t] <- mean.Bt*tau[t] # lognormal error, E=1
    # obs
    mean.It <- qI*Bt[t]
    It[t] <- mean.It*eps[t] # lognormal error, E=1
    mean.Jt <- qR*Rt[t]
    Jt[t] <- mean.Jt*ups[t] # lognormal error, E=1
    mean.Lt <- mt[t]*S*(S/2*Nt[t-1]+(1-S/2)*Nt[t])
    Lt[t] <- mean.Lt*kap[t] # lognormal error, E=1
    mean.Ct <- Ct[t-1]/Bt[t-1]*Bt[t]*(Bt[t-1]/(a*Bt[1]/2))^chi
    Ct[t] <- mean.Ct*ccc[t] # lognormal error, E=1
    if ((Bt[t]-Ct[t]) <= tol.depletion){
      if (warn){warning(paste0('Stock got depleted at t = ',t,'.'))}
      stockdepleted <- T
    }
    t <- t+1
  }
  
  #/////////////////////////////////////////////////////////////////////////////
  ##### Outputs ####
  #/////////////////////////////////////////////////////////////////////////////
  
  return(list('mt'=mt,'Rt'=Rt,'Bt'=Bt, # randeff
              'It'=It,'Jt'=Jt,'Lt'=Lt,'Ct'=Ct, # responses
              'gt'=gt,'gRt'=gRt,'Nt'=Nt)) # covariates
}
# END SS2v4simul
