
# PBA 1-14-2019
# Run a series of simulations to find the ESS germination
# fraction and ESS germination variance for one environment.

rm(list=ls())
setwd("C:/Repos/climate-strategies/code")

# read in functions
library(boot)
source("model_functions.R")

# set parameters for variable, unpredictable environment
myFec=50
myFecSigma = 30   # this makes it variable
myRho = 0       # 0 =  unpredictable, 1 = perfectly predictable
myAlpha = 1
mySeedSurv = 0.9
tot_time = 5000 # length of each simulation
burn_in = tot_time/5 + 1

# find ESS g fraction given a series of g variances
gVars = c(0,0.1,0.5,1,2)
ESS_gFrac_given_gVar = rep(NA,length(gVars))
for(iGvar in 1:length(gVars)){
  mySigmaG = gVars[iGvar]
  print(paste0("Doing ",iGvar," in ",length(gVars), " g vars"))
  source("find_gFrac.R")
  ESS_gFrac_given_gVar[iGvar] = ESS
}

# find ESS g var given a series of g fractions
gFracs = logit(seq(0.1,0.9,0.1))
ESS_gVar_given_gFrac = rep(NA,length(gFracs))
for(iGfrac in 1:length(gFracs)){
  myGfrac = gFracs[iGfrac]
  print(paste0("Doing ",iGfrac," in ",length(gFracs), " g vars"))
  source("find_gVar.R")
  ESS_gVar_given_gFrac[iGfrac] = ESS
}

# plot all partial ESS's
plot(gVars,ESS_gFrac_given_gVar,type="o",ylim=c(0,1),pch=16,col="blue",
     xlab="Var(g)",ylab="g",main=paste("Var(F)=",myFecSigma,", rho=",myRho))
points(ESS_gVar_given_gFrac,inv.logit(gFracs),pch=16,col="red")
lines(ESS_gVar_given_gFrac,inv.logit(gFracs),col="red")
