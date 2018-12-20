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
source("model_functions.R")

# set constant parameters
myFec=50
myFecSigma = 10
myRho = 0
myAlpha = 1
mySeedSurv = 0.8
tot_time = 10000 # length of each simulation
burn_in = tot_time/5 + 1

# set germination variance
mySigmaG = 0.2

# loop through germination fraction grid
gRes = gInv = seq(-4,3,0.5)
res_grid = rbar_grid = matrix(NA,length(gRes),length(gInv))
for(iR in 1:length(gRes)){
  for(iI in 1:length(gInv)){
    
   # simulate rates
   rates = get_F_G(tot_time, mu=c(myFec,gRes[iR],gInv[iI]),sigma=c(myFecSigma,mySigmaG,mySigmaG), rho=myRho)
   rates=data.frame(rates)
   
   #simulate resident population
   seeds_res = rep(NA, tot_time)
   plants_res = rep(NA, tot_time)
   r_inv = rep(NA, tot_time)
   seeds_res[1] = 2   # initial population
   plants_res[1] = 2   # initial population
   for(i in 2:tot_time){
     out = grow_res(seeds_res=seeds_res[i-1],Fec=rates$Fec[i],alpha=myAlpha,seedSurv=mySeedSurv,G_res=rates$G1[i])
     seeds_res[i] = out$seeds
     plants_res[i] = out$plants
   } 
   res_grid[iR,iI] = mean(seeds_res[burn_in:tot_time])
     
   # get invader growth rate (if the residents persist)
   if(res_grid[iR,iI]>0.01){
      r_inv = grow_inv(plants_res,Fec=rates$Fec[i],alpha=myAlpha,seedSurv=mySeedSurv,rates$G2[i])
      r_bar = mean(r_inv[burn_in:tot_time])
      rbar_grid[iR,iI] = r_bar
   }else{
      rbar_grid[iR,iI] = NA
   }

  }
}







