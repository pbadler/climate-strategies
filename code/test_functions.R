# bits of code to test the functoins

rm(list=ls())

# read in functions
source("model_functions.R")

# simulate resident population (deterministic, constant)
tot_time = 500
seeds_res = rep(NA, tot_time)
plants_res = rep(NA, tot_time)
r_inv = rep(NA, tot_time)
seeds_res[1] = 2   # initial population
plants_res[1] = NA   # initial population
for(i in 2:tot_time){
  out = grow_res(seeds_res=seeds_res[i-1],Fec=50,alpha=1,seedSurv=0.9,G_res=0.5)
  seeds_res[i] = out$seeds
  plants_res[i] = out$plants
}
plot(seeds_res,type="l")
plot(plants_res, type="l")

# test invader function
r_inv = grow_inv(plants_res,Fec=50,alpha=1,seedSurv=0.9,G_inv=0.5)
r_bar = mean(r_inv[(tot_time/4):tot_time])

# test generation of fecundity and germination rates
rates = get_F_G(tot_time, mu=c(50,0,0),sigma=c(10,0.2,0.1), rho=0.9)
plot(rates[,c(1,2)])
plot(rates[,c(1,3)])
plot(rates[,c(2,3)])
rates=data.frame(rates)

# do population growth with stochastic environment
tot_time = 500
burn_in = tot_time/4 + 1
seeds_res = rep(NA, tot_time)
plants_res = rep(NA, tot_time)
r_inv = rep(NA, tot_time)
seeds_res[1] = 2   # initial population
plants_res[1] = NA   # initial population
for(i in 2:tot_time){
  out = grow_res(seeds_res=seeds_res[i-1],Fec=rates$Fec[i],alpha=1,seedSurv=0.9,G_res=rates$G1[i])
  seeds_res[i] = out$seeds
  plants_res[i] = out$plants
}
plot(seeds_res,type="l")
plot(plants_res, type="l")

# test invader function
r_inv = grow_inv(plants_res,Fec=50,alpha=1,seedSurv=0.9,rates$G1[i])
r_bar = mean(r_inv[burn_in:tot_time])

# plot "precip: production" function
plot(rates$Fec[burn_in:tot_time], plants_res[burn_in:tot_time], 
     xlab="Fecundity=Precip",ylab="Production=Density")
summary(lm(plants_res[burn_in:tot_time]~rates$Fec[burn_in:tot_time]))
