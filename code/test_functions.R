# bits of code to test the functoins

# read in functions
source("model_functions.R")


# test population growth models
tot_time = 500
N = rep(NA, tot_time)
r_inv = rep(NA, tot_time)
N[1] = 2   # initial population
for(i in 2:tot_time){
  out = grow(N_res=N[i-1],Fec=50,alpha=1,seedSurv=0.9,G_res=0.5,G_inv=0.5)
  N[i] = out[1]
  r_inv[i] = out[2]
}
plot(N,type="l")
plot(r_inv,type="l")

# test generation of fecundity and germination rates
rates = get_F_G(tot_time, c(50,0,0),c(10,0.5,0.1), 0.9)
plot(rates[,c(1,2)])
plot(rates[,c(1,3)])
plot(rates[,c(2,3)])
rates=data.frame(rates)

# do population growth with stochastic environment
tot_time = 500
N = rep(NA, tot_time)
lambda_inv = rep(NA, tot_time)
N[1] = 2   # initial population
for(i in 2:tot_time){
  out = grow(N_res=N[i-1],Fec=rates$Fec[i],alpha=1,seedSurv=0.9,G_res=rates$G1[i],G_inv=rates$G2[i])
  N[i] = out[1]
  lambda_inv[i] = out[2]
}
plot(N,type="l")
plot(lambda_inv,type="l")
mean(lambda_inv[(tot_time/4):tot_time])
