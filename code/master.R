
# PBA 1-14-2019
# Run a series of simulations to find the ESS germination
# fraction and ESS germination variance for one environment.

rm(list=ls())
setwd("C:/Repos/climate-strategies/code")

# read in functions
library(boot)
source("model_functions.R")

### set parameters for variable, unpredictable environment ---------------------------------
myFec=2
myFecSigma = 0.5*myFec   # this makes it variable
myRho = 0          # 0 =  unpredictable, 1 = perfectly predictable
myAlpha = 0.04
mySeedSurv = 0.9
tot_time = 5000 # length of each simulation
burn_in = tot_time/5 + 1

### find ESS g fraction given a series of g variances
gVars = c(0,0.1,0.5,1,2,4)
ESS_gFrac_given_gVar = rep(NA,length(gVars))
for(iGvar in 1:length(gVars)){
  mySigmaG = gVars[iGvar]
  print(paste0("Doing ",iGvar," in ",length(gVars), " g vars"))
  source("find_gFrac.R")
  print(paste0("ESS - ESS2 = ", ESS-ESS2)) # compare two methods for determining ESS
  ESS_gFrac_given_gVar[iGvar] = ESS
}

### find ESS g var given a series of g fractions
gFracs = logit(seq(0.3,0.9,0.1))
ESS_gVar_given_gFrac = rep(NA,length(gFracs))
for(iGfrac in 1:length(gFracs)){
  myGfrac = gFracs[iGfrac]
  print(paste0("Doing ",iGfrac," in ",length(gFracs), " g vars"))
  source("find_gVar.R")
  ESS_gVar_given_gFrac[iGfrac] = ESS
}

### plot all partial ESS's
plot(gVars,ESS_gFrac_given_gVar,type="o",ylim=c(0,1),pch=16,col="blue",
     xlab="Var(g)",ylab="g",main=paste("Var(F)=",myFecSigma,", rho=",myRho))
points(ESS_gVar_given_gFrac,inv.logit(gFracs),pch=16,col="red")
lines(ESS_gVar_given_gFrac,inv.logit(gFracs),col="red")

### simulate and plot the precip - production relationship
gFrac = 0.68 # eyeballed from previous figure
gVar = 0   # eyeballed from previous figure
# simulate rates
rates = get_F_G(tot_time, mu=c(myFec,gFrac,gFrac),sigma=c(myFecSigma,gVar,gVar), rho=myRho)
rates=data.frame(rates)
#simulate resident population
production = seeds = rep(NA, tot_time)
seeds[1] = 2   # initial population
for(i in 2:tot_time){
  out = grow_res(seeds_res=seeds[i-1],Fec=rates$Fec[i],alpha=myAlpha,seedSurv=mySeedSurv,G_res=rates$G1[i],G_inv=rates$G2[i])
  seeds[i] = out$seeds
  production[i] = out$yield
}

print(paste("Resource mean=",round(mean(rates$Fec[burn_in:tot_time]),2)))
print(paste("Resource SD=",round(sd(rates$Fec[burn_in:tot_time]),2)))
print(paste("Production mean=",round(mean(production[burn_in:tot_time]),2)))
print(paste("Production SD=",round(sd(production[burn_in:tot_time]),2)))

sensitivity = lm(production[burn_in:tot_time]~rates$Fec[burn_in:tot_time])

plot(production[burn_in:tot_time]~rates$Fec[burn_in:tot_time],xlab="Resources",ylab="Production",
     xlim=c(0,3*myFec),ylim=c(0,100),
    main=paste("Var(F)=",myFecSigma,", rho=",myRho))
abline(sensitivity, col="red",lwd=2)
legend("topleft",paste("Slope=",round(coef(sensitivity)[2],2)),bty="n")

# look for lag effects
#plot(seeds[burn_in:(tot_time-1)],production[(burn_in+1):tot_time])
plot(production[burn_in:(tot_time-1)],production[(burn_in+1):tot_time],
     xlab="Production last year",ylab="Production this year")


### set parameters for variable, predictable environment ---------------------------------
myFec=2
myFecSigma = 0.5*myFec   # this makes it variable
myRho = 0.8          # 0 =  unpredictable, 1 = perfectly predictable
myAlpha = 0.04
mySeedSurv = 0.9
tot_time = 5000 # length of each simulation
burn_in = tot_time/5 + 1

### find ESS g fraction given a series of g variances
gVars = c(0,0.5,1,2,4,6)
ESS_gFrac_given_gVar = rep(NA,length(gVars))
for(iGvar in 1:length(gVars)){
  mySigmaG = gVars[iGvar]
  print(paste0("Doing ",iGvar," in ",length(gVars), " g vars"))
  source("find_gFrac.R")
  ESS_gFrac_given_gVar[iGvar] = ESS
}

### find ESS g var given a series of g fractions
tot_time = 10000 # length of each simulation
gFracs = logit(seq(0.3,0.9,0.1))
ESS_gVar_given_gFrac = rep(NA,length(gFracs))
for(iGfrac in 1:length(gFracs)){
  myGfrac = gFracs[iGfrac]
  print(paste0("Doing ",iGfrac," in ",length(gFracs), " g vars"))
  source("find_gVar.R")
  ESS_gVar_given_gFrac[iGfrac] = ESS
}

### plot all partial ESS's
plot(gVars,ESS_gFrac_given_gVar,type="o",ylim=c(0,1),pch=16,col="blue",
     xlab="Var(g)",ylab="g",main=paste("Var(F)=",myFecSigma,", rho=",myRho))
points(ESS_gVar_given_gFrac,inv.logit(gFracs),pch=16,col="red")
lines(ESS_gVar_given_gFrac,inv.logit(gFracs),col="red")

### simulate and plot the precip - production relationship
gFrac = 0.8 # eyeballed from previous figure
gVar = 3.9   # eyeballed from previous figure
# simulate rates
rates = get_F_G(tot_time, mu=c(myFec,gFrac,gFrac),sigma=c(myFecSigma,gVar,gVar), rho=myRho)
rates=data.frame(rates)
#simulate resident population
production = seeds = rep(NA, tot_time)
seeds[1] = 2   # initial population
for(i in 2:tot_time){
  out = grow_res(seeds_res=seeds[i-1],Fec=rates$Fec[i],alpha=myAlpha,seedSurv=mySeedSurv,G_res=rates$G1[i],G_inv=rates$G2[i])
  seeds[i] = out$seeds
  production[i] = out$yield
}

print(paste("Resource mean=",round(mean(rates$Fec[burn_in:tot_time]),2)))
print(paste("Resource SD=",round(sd(rates$Fec[burn_in:tot_time]),2)))
print(paste("Production mean=",round(mean(production[burn_in:tot_time]),2)))
print(paste("Production SD=",round(sd(production[burn_in:tot_time]),2)))

sensitivity = lm(production[burn_in:tot_time]~rates$Fec[burn_in:tot_time])

plot(production[burn_in:tot_time]~rates$Fec[burn_in:tot_time],xlab="Resources",ylab="Production",
     xlim=c(0,3*myFec),ylim=c(0,100),
    main=paste("Var(F)=",myFecSigma,", rho=",myRho))
abline(sensitivity, col="red",lwd=2)
legend("topleft",paste("Slope=",round(coef(sensitivity)[2],2)),bty="n")

# look for lag effects
#plot(seeds[burn_in:(tot_time-1)],production[(burn_in+1):tot_time])
plot(production[burn_in:(tot_time-1)],production[(burn_in+1):tot_time],
     xlab="Production last year",ylab="Production this year")

### set parameters for low variability environment ---------------------------------
myFec=2
myFecSigma = 0.05*myFec   # this makes it variable
myRho = 0.4       # 0 =  unpredictable, 1 = perfectly predictable
myAlpha = 0.04
mySeedSurv = 0.9
tot_time = 5000 # length of each simulation
burn_in = tot_time/5 + 1

### find ESS g fraction given a series of g variances
gVars = c(0,0.1,0.5,1,2)
ESS_gFrac_given_gVar = rep(NA,length(gVars))
for(iGvar in 1:length(gVars)){
  mySigmaG = gVars[iGvar]
  print(paste0("Doing ",iGvar," in ",length(gVars), " g vars"))
  source("find_gFrac.R")
  ESS_gFrac_given_gVar[iGvar] = ESS
}

### find ESS g var given a series of g fractions
gFracs = logit(seq(0.1,0.9,0.1))
ESS_gVar_given_gFrac = rep(NA,length(gFracs))
for(iGfrac in 1:length(gFracs)){
  myGfrac = gFracs[iGfrac]
  print(paste0("Doing ",iGfrac," in ",length(gFracs), " g vars"))
  source("find_gVar.R")
  ESS_gVar_given_gFrac[iGfrac] = ESS
}

### plot all partial ESS's
plot(gVars,ESS_gFrac_given_gVar,type="o",ylim=c(0,1),pch=16,col="blue",
     xlab="Var(g)",ylab="g",main=paste("Var(F)=",myFecSigma,", rho=",myRho))
points(ESS_gVar_given_gFrac,inv.logit(gFracs),pch=16,col="red")
lines(ESS_gVar_given_gFrac,inv.logit(gFracs),col="red")

### simulate and plot the precip - production relationship
gFrac = 1 # eyeballed from previous figure
gVar = 0   # eyeballed from previous figure
# simulate rates
rates = get_F_G(tot_time, mu=c(myFec,gFrac,gFrac),sigma=c(myFecSigma,gVar,gVar), rho=myRho)
rates=data.frame(rates)
#simulate resident population
production = seeds = rep(NA, tot_time)
seeds[1] = 2   # initial population
for(i in 2:tot_time){
  out = grow_res(seeds_res=seeds[i-1],Fec=rates$Fec[i],alpha=myAlpha,seedSurv=mySeedSurv,G_res=rates$G1[i],G_inv=rates$G2[i])
  seeds[i] = out$seeds
  production[i] = out$yield
}
sensitivity = lm(production[burn_in:tot_time]~rates$Fec[burn_in:tot_time])
plot(production[burn_in:tot_time]~rates$Fec[burn_in:tot_time],xlab="Resources",ylab="Production",
     xlim=c(0,3*myFec),ylim=c(0,100),
    main=paste("Var(F)=",myFecSigma,", rho=",myRho))
abline(sensitivity, col="red",lwd=2)
legend("topleft",paste("Slope=",round(coef(sensitivity)[2],2)),bty="n")

print(paste("Resource mean=",round(mean(rates$Fec[burn_in:tot_time]),2)))
print(paste("Resource SD=",round(sd(rates$Fec[burn_in:tot_time]),2)))
print(paste("Production mean=",round(mean(production[burn_in:tot_time]),2)))
print(paste("Production SD=",round(sd(production[burn_in:tot_time]),2)))
