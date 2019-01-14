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
myFecSigma = 30
myRho = 0.8
myAlpha = 1
mySeedSurv = 0.9
tot_time = 5000 # length of each simulation
burn_in = tot_time/5 + 1

# set germination variance
mySigmaG = 1

# loop through germination fraction grid
gRes = logit(seq(0.1,0.9,0.02))
gInv_step = 0.01
res_grid = rbar_grid = matrix(NA,5,length(gRes))
counter=0
for(iR in 1:length(gRes)){
  gInv = logit(seq(inv.logit(gRes[iR]) - 2*gInv_step, inv.logit(gRes[iR]) + 2*gInv_step, gInv_step))
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
   res_grid[iI,iR] = mean(seeds_res[burn_in:tot_time])
   r_bar = mean(r_inv[burn_in:(tot_time-1)])
   rbar_grid[iI,iR] = r_bar
   
   # report progress
   counter = counter + 1
   print(paste(counter,"out of",length(gRes)*length(gInv),"complete"))
   
  }
}

delta_g = seq(-2*gInv_step,2*gInv_step,gInv_step)
contour(inv.logit(gRes),delta_g,t(rbar_grid),xlab="Resident g",ylab="delta g")
image(inv.logit(gRes),delta_g,t(rbar_grid),xlab="Resident g",ylab="delta g")

plot(delta_g,rbar_grid[,3])
matplot(delta_g,rbar_grid,type="l")
abline(h=0)

slopes = numeric(length(gRes))
for(i in 1:length(gRes)){
  slopes[i] = coef(lm(rbar_grid[,i]~delta_g))[2]
}
slope_of_slopes = lm(slopes~gRes)

plot(gRes,slopes)
abline(slope_of_slopes)

ESS = inv.logit(-1*(coef(slope_of_slopes)[1]/coef(slope_of_slopes)[2]))
