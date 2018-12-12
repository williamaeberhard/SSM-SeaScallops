#///////////////////////////////////////////////////////////////
#### SSM-SeaScallops: Simulate data and fit alternative SSM ####
#///////////////////////////////////////////////////////////////

# Reference for alternative SSM for Bay of Fundy SPA4 sea scallops fishery:
#   Yin, Y., Aeberhard, W. H., Smith, S. J., and Mills Flemming, J. (2018).
#   Identifiable state-space models: A case study of the Bay of Fundy sea
#   scallop fishery. Canadian Journal of Statistics, In press.


rm(list=ls())

### load libraries and functions
library(TMB)
compile('SS2v4.cpp') # compile TMB C++ template, to do only once
dyn.load(dynlib("SS2v4")) # load the compiled dynamic library

source('SS2v4simul.r') # creates SS2v4simul function, simulates data


### create design following simulation study in paper
NY <- 23 # number of time steps (years)
yearsvec <- 1:NY

set.seed(1234) # fix design for reproducibility
gt <- rep(1.3,NY) # growth rate commercial size
gRt <- rep(1.8,NY) # growth rate recruitment size
Nt <- exp(16.8+as.numeric(arima.sim(model=list(ar=c(1.2,-0.5)),n=NY,sd=0.3)))
# ^ survey estimates of live scallops, ARMA mimics well true data

sigmatau <- 0.1 # sd of log-normal error in process eq for Bt
sigmaphi <- 0.1 # sd of log-normal error in random walk for Rt
sigmam <- 0.1 # sd of log-normal error in random walk for mt
sigmaepsilon <- 0.1 # sd of log-normal error in obs eq for It
sigmaupsilon <- 0.1 # sd of log-normal error in obs eq for I^R_t
sigmakappa <- 0.1 # sd of log-normal error in obs eq for Lt
sigmaC <- 0.1 # sd of log-normal error in obs eq for Ct

qI <- 0.4 # catchability coeff in obs eq for It, commercial size
qR <- 0.2 # catchability coeff in obs eq for I^R_t, recruitment size
S <- 0.5 # dissolution rate, average hinge separation time of a clapper (years)
a <- 1.6 # ratio of biomass at bioeconomic equilibrium to biomass at MSY
chi <- 0.1 # rate at which effort enters or exits the fishery

theta.true <- c(sigmatau,sigmaphi,sigmam,
                sigmaepsilon,sigmaupsilon,sigmakappa,sigmaC,
                qI,qR,S,a,chi)
length.theta <- length(theta.true) # 12

B1 <- 1000 # initial condition for Bt, rough value
R1 <- 0.2*B1 # initial condition for Rt, rough value
m1 <- 0.2 # initial condition for mt, rough value
C1 <- 0.2*B1 # initial condition for Ct, rough value


### simul data
dat <- SS2v4simul(gt=gt,gRt=gRt,Nt=Nt,
                  B1=B1,R1=R1,m1=m1,C1=C1,
                  sigmatau=sigmatau,sigmaphi=sigmaphi,sigmam=sigmam,
                  sigmaeps=sigmaepsilon,sigmaups=sigmaupsilon,
                  sigmakap=sigmakappa,sigmaC=sigmaC,
                  qI=qI,qR=qR,S=S,a=a,chi=chi,seedvalue=NULL)
str(dat) # T=23 observations, all variables in there, see eq. (6)-(12) in paper


### fit alternative SSM to simulated data
datalist <- list('I'=dat$It,'IR'=dat$Jt,'L'=dat$Lt,'C'=dat$Ct,
                 'g'=dat$gt,'gR'=dat$gRt,'N'=dat$Nt)
# ^ all "fixed" variables, both response variable and covariates

parlist <- list('log_sigma_tau'=0,'log_sigma_phi'=0,'log_sigma_m'=0,
                'log_sigma_epsilon'=0,'log_sigma_upsilon'=0,
                'log_sigma_kappa'=0,'log_sigma_C'=0,
                'log_q_I'=-1,'log_q_R'=-1,
                'log_S'=0,'log_a'=0,'log_chi'=-1,
                'log_B'=rep(log(max(dat$It)*10),NY),
                'log_R'=rep(log(max(dat$Jt)*10),NY),
                'log_m'=rep(log(0.5),NY))
# ^ starting values for all model parameters, incl randeff

system.time(obj <- MakeADFun(data=datalist,parameters=parlist,
                             random=c('log_B','log_R','log_m'),
                             DLL="SS2v4",silent=T))
# ^ create TMB object that contains Laplace-approximated marginal negative log-
#   likelihood and its gradient wrt theta as functions

obj$fn() # Laplace-approximated marginal neg log-lik
obj$gr() # corresponding gradient wrt to theta (not incl randeff)

system.time(opt <- try(nlminb(start=obj$par,obj=obj$fn,gr=obj$gr,
                              control=list(eval.max=1000,iter.max=1000)),T))
# ^ minimize wrt theta

opt # apparently converged fine

rep <- try(sdreport(obj,bias.correct=F,bias.correct.control=list(sd=F)),T)
# ^ reports all estimates+se, incl quantities declared in ADREPORT in SS2v4.cpp
summary.rep <- summary(rep)

cbind(theta.true,summary.rep[(length.theta+3*NY+1):(2*length.theta+3*NY),])
# ^ estimates+se of all model parameters (theta), original scale


### plot predicted randeff against true, incl envelope from +-1.96se
logmt <- summary.rep[dimnames(summary.rep)[[1]]=='log_m',1]
sd.logmt <- summary.rep[dimnames(summary.rep)[[1]]=='log_m',2]
logBt <- summary.rep[dimnames(summary.rep)[[1]]=='log_B',1]
sd.logBt <- summary.rep[dimnames(summary.rep)[[1]]=='log_B',2]
logRt <- summary.rep[dimnames(summary.rep)[[1]]=='log_R',1]
sd.logRt <- summary.rep[dimnames(summary.rep)[[1]]=='log_R',2]
# ^ extract pred and se on log scale

pred.m <- exp(logmt)
lb.m <-  exp(logmt-1.96*sd.logmt)
ub.m <-  exp(logmt+1.96*sd.logmt)
pred.B <- exp(logBt)
lb.B <-  exp(logBt-1.96*sd.logBt)
ub.B <-  exp(logBt+1.96*sd.logBt)
pred.R <- exp(logRt)
lb.R <-  exp(logRt-1.96*sd.logRt)
ub.R <-  exp(logRt+1.96*sd.logRt)
# ^ rough 95% CI on original scale

colmed <- '#2b05ff'
colenv <- paste0(colmed,'30')

par(mfrow=c(3,1))
# Bt
plot(yearsvec,dat$Bt,type='l',lwd=2,
     xlab='Years',ylab=expression(italic(B[t])~'(log scale)'),
     main=expression('Predicted commercial biomass'~italic(B[t])),
     ylim=c(100,6000),log='y')
grid(nx=NA,ny=NULL,equilogs=F)
polygon(c(yearsvec,yearsvec[NY:1]),
        c(lb.B,ub.B[NY:1]),col=colenv,border=NA,xpd=F)
lines(yearsvec,pred.B,type='o',lty=2,col=colmed)
# Rt
plot(yearsvec,dat$Rt,type='l',lwd=2,
     xlab='Years',ylab=expression(italic(R[t])~'(log scale)'),
     main=expression('Predicted recruitment biomass'~italic(R[t])),
     ylim=c(100,1000),log='y')
grid(nx=NA,ny=NULL,equilogs=F)
polygon(c(yearsvec,yearsvec[NY:1]),
        c(lb.R,ub.R[NY:1]),col=colenv,border=NA,xpd=F)
lines(yearsvec,pred.R,type='o',lty=2,col=colmed)
# mt
plot(yearsvec,dat$mt,type='l',lwd=2,
     xlab='Years',ylab=expression(italic(m[t])),
     main=expression('Predicted natural mortality rate'~italic(m[t])),
     ylim=c(0,1))
grid(nx=NA,ny=NULL,equilogs=F)
polygon(c(yearsvec,yearsvec[NY:1]),
        c(lb.m,ub.m[NY:1]),col=colenv,border=NA,xpd=F)
lines(yearsvec,pred.m,type='o',lty=2,col=colmed)
par(mfrow=c(1,1))
# ^ for that particular sample: recruitment biomass R_t generally overestimated,
#   (stock less productive than predicted), but compensated by overestimated
#   natural mortality m_t, so that commercial biomass B_t is well predicted.

### END SS2v4_main
