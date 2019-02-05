
# For one germination fraction and given environment, find ESS germination variance
# Call from master.R

###
### first, determine if ESS gVar > 0 by setting resident gVar to zero and invader > 0
###

# simulate rates
rates = get_F_G(tot_time*2, mu=c(myFec,myGfrac,myGfrac),sigma=c(myFecSigma,0,0.02), rho=myRho)
rates=data.frame(rates)
   
#simulate resident population
seeds_res = rep(NA, tot_time*2)
r_inv = rep(NA, tot_time*2)
seeds_res[1] = 2   # initial population
for(i in 2:(tot_time*2)){
  out = grow_res(seeds_res=seeds_res[i-1],Fec=rates$Fec[i],alpha=myAlpha,seedSurv=mySeedSurv,G_res=rates$G1[i],G_inv=rates$G2[i])
  seeds_res[i] = out$seeds
  r_inv[i] = out$r_inv
}
nonzero = mean(r_inv[burn_in:(2*tot_time-1)]) > 0.0005  

if(nonzero == FALSE ){
  
  ESS = 0
  
}else{ # if ESS gVar > 0, then find optimum
  
  # loop through germination variance grid
  sigmaRes = c(0.05,0.1,0.5,1,2,3,4,5,6)
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
  # contour(sigmaRes,delta_var,t(rbar_grid),xlab="Resident var(g)",ylab="delta var(g)")
  #image(sigmaRes,delta_var,t(rbar_grid),xlab="Resident var(g)",ylab="delta var(g)")
  
  #plot(delta_var,rbar_grid[,1])
  #matplot(delta_var,rbar_grid,type="l")
  
  slopes = numeric(length(sigmaRes))
  for(i in 1:length(sigmaRes)){
    slopes[i] = coef(lm(rbar_grid[,i]~delta_var))[2]
  }
  
  if(sum(slopes<0)==0){
    stop("No negative slopes")
  }else if(sum(slopes>0)==0){
    stop("No positive slopes")
  }else{
    # fit spline through slopes
    ss = smooth.spline(x=sigmaRes,y=slopes,df=4)
    # plot(sigmaRes,slopes)
    # lines(predict(ss))
    ss_fun = function(x) predict(ss, x)$y 
    # find where slope spline crosses zero
    ESS = uniroot(ss_fun,interval=c(min(sigmaRes),max(sigmaRes)))$root  
    
    # alternative approach: find where max rbar is closest to zero
    rbar_max = apply(rbar_grid,2,FUN="max")
    ESS2 = sigmaRes[which(abs(rbar_max)==min(abs(rbar_max)))]
    
    # plot both
   old_par = par(no.readonly=TRUE)
   par(mfrow=c(2,1),mar=c(2.5,4,0,1),oma=c(2,0,1,0))
   plot(sigmaRes,slopes)
   lines(predict(ss))
   abline(v=ESS,col="red")
   abline(h=0,lty="dotted")
   plot(sigmaRes,rbar_max)
   abline(h=0,lty="dotted")
   abline(v=ESS2,col="red")
   mtext("sigmaRes",side=1,outer=T)
   par(old_par)
   
  }
  
}


