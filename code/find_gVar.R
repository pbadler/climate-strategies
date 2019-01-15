
# For one germination fraction and given environment, find ESS germination variance
# Call from master.R

# loop through germination variance grid
sigmaRes = c(0.05,0.1,0.5,1,2,4)
sigma_steps = c(0.5,1,2)  # factors to multiply sigmaRes
res_grid = rbar_grid = matrix(NA,3,length(sigmaRes))
counter=0
for(iR in 1:length(sigmaRes)){
  sigmaInv = sigma_steps*sigmaRes[iR]
  for(iI in 1:length(sigmaInv)){
    
   # simulate rates
   rates = get_F_G(tot_time, mu=c(myFec,myGfrac,myGfrac),sigma=c(myFecSigma,sigmaRes[iR],sigmaInv[iI]), rho=myRho)
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
   print(paste(counter,"out of",length(sigmaRes)*length(sigmaInv),"complete"))
   
  }
}

delta_var = sigma_steps
contour(sigmaRes,delta_var,t(rbar_grid),xlab="Resident var(g)",ylab="delta var(g)")
image(sigmaRes,delta_var,t(rbar_grid),xlab="Resident var(g)",ylab="delta var(g)")

plot(delta_var,rbar_grid[,1])
matplot(delta_var,rbar_grid,type="l")
#abline(h=0)

slopes = numeric(length(sigmaRes))
for(i in 1:length(sigmaRes)){
  slopes[i] = coef(lm(rbar_grid[,i]~delta_var))[2]
}
slope_of_slopes = lm(slopes~sqrt(sigmaRes) + sigmaRes)

# START HERE
polyroot(matrix(coef(slope_of_slopes),ncol=1))

plot(sqrt(sigmaRes),slopes)
lines(sqrt(sigmaRes),predict(slope_of_slopes))

if(sum(slopes>0)==0){
  ESS = 0  # always better to have less var(g)
}else{
  ESS = -1*(coef(slope_of_slopes)[1]/coef(slope_of_slopes)[2])
}



