
# For one germination variance and given environment, find ESS germination fraction
# Call from master.R

# loop through germination fraction grid
gRes = logit(seq(0.05,0.95,0.05))
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
   if(counter%%10==0) print(paste(counter,"out of",length(gRes)*length(gInv),"complete"))
   
  }
}

delta_g = seq(-2*gInv_step,2*gInv_step,gInv_step)
#contour(inv.logit(gRes),delta_g,t(rbar_grid),xlab="Resident g",ylab="delta g")
#image(inv.logit(gRes),delta_g,t(rbar_grid),xlab="Resident g",ylab="delta g")

#plot(delta_g,rbar_grid[,3])
#matplot(delta_g,rbar_grid,type="l")
#abline(h=0)

slopes = numeric(length(gRes))
for(i in 1:length(gRes)){
  slopes[i] = coef(lm(rbar_grid[,i]~delta_g))[2]
}

if(sum(slopes<0)==0){
  ESS = 1
}else if(sum(slopes>0)==0){
  ESS = 0
}else{
  # fit spline through slopes
  ss = smooth.spline(x=gRes,y=slopes,df=5)
  #plot(gRes,slopes)
  #lines(predict(ss))
  ss_fun = function(x) predict(ss, x)$y 
  # find where slope spline crosses zero
  tmp = uniroot(ss_fun,interval=c(min(gRes),max(gRes)))$root
  ESS = inv.logit(tmp)
}





