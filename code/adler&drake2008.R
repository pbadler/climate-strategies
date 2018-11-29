# Main script used for Adler and Drake 2008 Am Nat


# PARAMETERS -----------------------------------------
# Name of output file:
outfile = "results_new_1.csv"

# Environment (determines germination rates):
envMean = 0   # mean value of the environment (0 gives 0.5 germination probability)
nugget = 0.25  # proportion of env. variation that is random
Nyears = 1000    # length of initial simulated time series (but each run goes to extinction)
reps = 1000       # replication of each parameter combination

# Plant competition model:
meanLambda = 100  # species 1 and 2 are deviations from this value
seedSurv = 0.5  # survival of ungerminated seeds
alpha = 1      # all competition coefficients (intra- and inter-) are equal

# Set up factorial experiment on parameters envVar, scale, rho, lambda
envVar = c(0,0.4,1,3,5)  # determines variance of germination rates
scale =0 # c(0,1,2)     # determines length of temporal autocorrelation (0=none)
rho =c(1,0.8,0.5,-0.5)   # between-species correlation in germination rates
lamDiff = c(1,2,5)    # difference in lambda between species
initMethod = c("mean","hi","low")       # how to initialize each run

#-----------------------------------------------------

experiment = expand.grid("envVar"=envVar,"scale"=scale,"rho"=rho,"lamDiff"=lamDiff,
  "initMethod"=initMethod)
rm(envVar,scale,rho,lamDiff)

# Set up output of simulation results
experiment$meanTimeCoexist=NA
experiment$medianTimeCoexist=NA
experiment$meanN1=NA
experiment$medianN1=NA
experiment$meanN2=NA
experiment$medianN2=NA
experiment$meanTimePersists1=NA
experiment$medianTimePersists1=NA
experiment$meanTimePersists2=NA
experiment$medianTimePersists2=NA
experiment$fractionExtinct1=NA
experiment$fractionExtinct2=NA
saveReps = array(NA,dim=c(reps,2,dim(experiment)[1]))

# -----THE MODEL------------------------------------------------
  updateN = function(N,lambda,alpha,seedSurv,G){
  	# This is an annual plant model, so it tracks
  	# populations as numbers of seeds

  	# N = vector of length 2 giving each population at time t
  	# lambda = vector of the 2 species' fecundities
    # alpha and seedSurv are  scalars
  	# G = vector of length 2 of each species' germination rate
  	# output = vector (2) giving each population at time t+1
  	tiny = 1      # extinction threshold
  	tmp1 = G[1]*N[1]
  	tmp2 = G[2]*N[2]
  	out1 = seedSurv*(1-G[1])*N[1]+lambda[1]*tmp1/(1+alpha*tmp1+alpha*tmp2)
  	out2 = seedSurv*(1-G[2])*N[2]+lambda[2]*tmp2/(1+alpha*tmp1+alpha*tmp2)
  	out = c(out1,out2)
  	out[out<tiny] = 0
  	# Or replace above line by drawing Poisson variates
  	# from distributions with means out 1 and out2 ?
  	return(out)
  }

getG = function(Nyears,envVar,scale,rho){
    # Simulate the environment -------------------------------
    # Generate 2 independent vectors with desired autocorrelation
    if(scale==0 | envVar==0){
      # no environmental variation OR no autocorrelation
      A =  matrix(rnorm(2*Nyears,0,sqrt(envVar)) ,Nyears,2)
    }else{
      model = "stable"
      a = 1   ## see RandomFields help("CovarianceFct") for additional parameters of the covariance functions
      x = seq(1, Nyears)
      A = GaussRF(x=x, model=model, grid=F,n=2,param=c(0, envVar, nugget, scale, a))
    }
    if(rho==1){
      # spp respond identically to environment
      B=cbind(A[,1],A[,1])
    }else{
      B=A
      # Create between species correlation
      B[,2]= B[,1]*rho+sqrt(1-rho^2)*B[,2]
    }
    # adjust for envMean if different from zero
    B = B + envMean
    # Inverse logit transform
    out = exp(B)/(1+exp(B))
    out
}

#-----MAIN LOOP-----------------------------------------------------
library(RandomFields)

for(iExp in 1:dim(experiment)[1]){
  # get parameters
  envVar = experiment$envVar[iExp]
  scale = experiment$scale[iExp]
  rho =  experiment$rho[iExp]
  lambda = c(meanLambda+experiment$lamDiff[iExp], meanLambda-experiment$lamDiff[iExp])
  # set up output
  persistence = matrix(NA,reps,2)
  meanN = matrix(NA,reps,2)
  for(iRep in 1:reps){

    germination=getG(Nyears,envVar,scale,rho)

    # Do one run of the model -------------------------------------
    N=data.frame("spp1"=1,"spp2"=1)
    # assign initial population sizes
    if(experiment$initMethod[iExp]=="mean"){
      N[1,] = mean(lambda)/2
    }else if(experiment$initMethod[iExp]=="hi"){
      N[1,] = mean(lambda)
    }else if(experiment$initMethod[iExp]=="low"){
      N[1,] = mean(lambda)/10
    }else if(experiment$initMethod[iExp]=="hilo"){
      N[1,] = c(mean(lambda),mean(lambda)/10)
    }else if(experiment$initMethod[iExp]=="lohi"){
      N[1,] = c(mean(lambda)/10,mean(lambda))
    }
    counter=1
    stopRun=F
    while(stopRun==F){
      counter=counter+1
    	N[counter,] = updateN(N[counter-1,],lambda,alpha,seedSurv,germination[counter-1,])
    	if(sum(N[counter,]==0)>=1) stopRun=T
    	# make new environmental sequence if necessary
      if((counter-1)==dim(germination)[1]){
         germination=rbind(germination,getG(Nyears,envVar,scale,rho))
      }
    }   # next i
    tmp = rep(dim(N)[1],2)
    tmp[which(N[dim(N)[1],]>0)]=NA
    persistence[iRep,] = tmp
    meanN[iRep,]=colMeans(N)

  } #next iRep
  saveReps[,,iExp]=persistence
  coexist=apply(persistence,MARGIN=1,FUN=min,na.rm=T)
  experiment$meanTimeCoexist[iExp]=mean(coexist,na.rm=T)
  experiment$medianTimeCoexist[iExp]=median(coexist,na.rm=T)
  experiment$meanN1[iExp]=mean(meanN[,1])
  experiment$medianN1[iExp]=median(meanN[,1])
  experiment$meanN2[iExp]=mean(meanN[,2])
  experiment$medianN2[iExp]=median(meanN[,2])
  experiment$meanTimePersists1[iExp]=mean(persistence[,1],na.rm=T)
  experiment$medianTimePersists1[iExp]=median(persistence[,1],na.rm=T)
  experiment$meanTimePersists2[iExp]=mean(persistence[,2],na.rm=T)
  experiment$medianTimePersists2[iExp]=median(persistence[,2],na.rm=T)
  experiment$fractionExtinct1[iExp]=sum(is.na(persistence[,1])==F)/reps
  experiment$fractionExtinct2[iExp]=sum(is.na(persistence[,2])==F)/reps

  tmp = paste("Parameter combination",iExp,"of",dim(experiment)[1],sep=" ")
  print(tmp)
  print(date())
  flush.console()
} # next iExp

# go back and get the standard error of the mean coexistence times
# and 95% confidence limits on the fraction of spp 1 extinctions
experiment$coexistTimeSE=NA
experiment$fracExt1Lower=NA
experiment$fracExt1Upper=NA
for(i in 1:dim(experiment)[1]){
  tmp1=saveReps[,,i]
  tmp2=apply(tmp1,MARGIN=1,min,na.rm=T)
  experiment$coexistTimeSE[i]=sd(tmp2)/sqrt(reps)

  tmp2=1*(is.na(tmp1[,1])==F)
  if(sum(tmp2)>0){
    tmp3=glm(tmp2~1,binomial)
    tmp4=inv.logit(confint(tmp3))
    experiment$fracExt1Lower[i]=tmp4[1]
    experiment$fracExt1Upper[i]=tmp4[2]
  }
}

write.table(experiment,outfile,row.names=F,sep=",")

# if only one run performed, make figures
if(dim(experiment)[1]==1 & reps==1){
  par(mfrow=c(1,3),tcl=-0.2)
  # plot germination rates
  matplot(germination,type="l",xlab="Time",ylab="Germination")
  # plot density time series
  matplot(N,type="l",xlab="Time",ylab="N")
  # plot frequency dependence in growth
  Nfreq = N/rowSums(N)  # calculate frequencies
  Nfreq[Nfreq==1]=NA # remove values after one species goes exinct
  growth = log(N[2:NROW(N),])-log(N[1:(NROW(N)-1),])
  growth[growth==-Inf]=NA
  myLims = c(min(growth,na.rm=T)-0.05,max(growth,na.rm=T)+0.05)
  plot(Nfreq[1:NROW(Nfreq)-1,1],growth[,1],xlab="Frequency",ylab="Growth rate",
    xlim=c(0,1),ylim=myLims,col="black")
  abline(lm(growth[,1]~ Nfreq[1:NROW(Nfreq)-1,1] ),col="black")
  par(new=T)
  plot(Nfreq[1:NROW(Nfreq)-1,2],growth[,2],xlab="",ylab="",xaxt="n",yaxt="n",
    xlim=c(0,1),ylim=myLims,col="red")
  abline(lm(growth[,2]~ Nfreq[1:NROW(Nfreq)-1,2] ),col="red")
  abline(0,0,lty="dotted")
  legend("topright",legend=c("Spp 1", "Spp 2"),lty="solid",col=c("black","red"),bty="n")
}

