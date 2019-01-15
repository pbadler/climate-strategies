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
  out = grow_res(seeds_res=seeds_res[i-1],Fec=50,alpha=1,seedSurv=0.9,G_res=0.5,G_inv=0.5)
  seeds_res[i] = out$seeds
}
plot(seeds_res,type="l")


tot_time = 1000

# test generation of fecundity and germination rates
rates = get_F_G(tot_time, mu=c(50,0,0),sigma=c(10,0.01,0.01), rho=0)
plot(rates[,c(1,2)])
plot(rates[,c(1,3)])
plot(rates[,c(2,3)])
rates=data.frame(rates)

# do population growth with stochastic environment
burn_in = tot_time/4 + 1
seeds_res = rep(NA, tot_time)
plants_res = rep(NA, tot_time)
r_inv = r_res = rep(NA, tot_time)
seeds_res[1] = 2   # initial population
for(i in 2:tot_time){
  out = grow_res(seeds_res=seeds_res[i-1],Fec=rates$Fec[i],alpha=1,seedSurv=0.9,G_res=rates$G1[i],G_inv=rates$G2[i])
  seeds_res[i] = out$seeds
  r_res[i] = out$r_res
  r_inv[i] = out$r_inv
}
plot(seeds_res,type="l")


# plot "precip: production" function
production=seeds_res[1:(tot_time-1)]*rates$G1[2:tot_time]
plot(rates$Fec[burn_in:(tot_time-1)], production[burn_in:(tot_time-1)], 
     xlab="Fecundity=Precip",ylab="Production=Density")
summary(lm(production[burn_in:tot_time]~rates$Fec[burn_in:tot_time]))
