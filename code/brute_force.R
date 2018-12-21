# PBA 12/20/2018
# Seeing how long it will take to do this by brute force.
# Pseudo code:
# 1. Create the environment (set variance in fecundity and rho)
# 2. Holding germination variance constant, find ESS germination fraction
#     a. create a coarse grid showing all combinations of invader and resident germ. fractions
#     b. for each combination, calculate invasion growth rate
#     c. find the combination where invasion growth rate is zero, and any change makes it worse
#     d. repeat a-c on finer grid
# If this works:
# 3. Repeat step for different values of germination variance
# 4. Repeat steps 2-3, holding germ fraction constant each time, changing germ variance
# 5. Hope that the lines cross

# bits of code to test the functoins

rm(list=ls())

# read in functions
library(boot)
source("model_functions.R")

# set constant parameters
myFec=50
myFecSigma = 20
myRho = 0
myAlpha = 1
mySeedSurv = 0.8
tot_time = 5000 # length of each simulation
burn_in = tot_time/5 + 1

# set germination variance
mySigmaG = 0

# loop through germination fraction grid
gRes = gInv = logit(seq(0.02,0.98,0.04))
res_grid = rbar_grid = matrix(NA,length(gRes),length(gInv))
counter=0
for(iR in 1:length(gRes)){
  for(iI in 1:length(gInv)){
    
   # simulate rates
   rates = get_F_G(tot_time, mu=c(myFec,gRes[iR],gInv[iI]),sigma=c(myFecSigma,mySigmaG,mySigmaG), rho=myRho)
   rates=data.frame(rates)
   
   #simulate resident population
   seeds_res = rep(NA, tot_time)
   r_inv = rep(NA, tot_time)
   seeds_res[1] = 2   # initial population
   for(i in 2:tot_time){
     out = grow_res(seeds_res=seeds_res[i-1],Fec=rates$Fec[i],alpha=myAlpha,seedSurv=mySeedSurv,G_res=rates$G1[i],G_inv=rates$G2[i])
     seeds_res[i] = out$seeds
     r_inv[i] = out$r_inv
   } 
   res_grid[iR,iI] = mean(seeds_res[burn_in:tot_time])
   r_bar = mean(r_inv[burn_in:(tot_time-1)])
   rbar_grid[iR,iI] = r_bar
   
   # report progress
   counter = counter + 1
   print(paste(counter,"out of",length(gRes)*length(gInv),"complete"))
   
  }
}

image(inv.logit(gRes),inv.logit(gInv),rbar_grid,xlab="Resident g",ylab="Invader g")

plot(inv.logit(gInv),rbar_grid[24,])
#matplot(t(rbar_grid),type="l")



