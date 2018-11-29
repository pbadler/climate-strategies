
library(mvtnorm)

# simulate Fecundity and Germination rates for two populations;
# both have the same fecundity, but can have differences in germination
# mean and variance. The F-G covariance is fixed for both.
get_F_G = function(Nyears, mu, sigma , rho){
    # mu is a vector of dim 3: mean rates for fecundity, germ pop 1, germ pop 2
    # sigma is a vector of dim 3: variance for fecundity, germ pop 1, germ pop 2
    # rho is the correlation of fecundity and germination
    
    # set up variance-covariance matrix
    vcov = matrix(0,3,3)
    diag(vcov) = sigma^2
    vcov[2,1] = vcov[1,2] = rho*sigma[1]*sigma[2]  # turn correlation into covariance
    vcov[3,1] = vcov[1,3] = rho*sigma[1]*sigma[3]
    vcov[2,3] = vcov[3,2] = sigma[2]*sigma[3]  # make G1 and G2 perfectly correlated
    
    # simulate values
    out = rmvnorm(Nyears, mu, vcov)
    
    # inverse logit transform on germination columns
    out[,2:3] = exp(out[,2:3])/(1+exp(out[,2:3])) 
    
    colnames(out)=c("Fec","G1","G2")
    return(out)
}


# The population model. This tracks the population size
# of a resident, and also calculates growth rate of an
# invader at each time step
grow = function(N_res,Fec,alpha,seedSurv,G_res,G_inv){
  	# This is an annual plant model, so it tracks
  	# populations as numbers of seeds

  	# N = population at time t
  	# Fec = fecundity
    # alpha and seedSurv are  scalars
  	# G_res and G_inv = germination rates
  	# output =  population at time t+1
  
    #update resident
  	tmp1 = G_res*N_res
  	N_res_new = seedSurv*(1-G_res)*N_res[1]+Fec*tmp1/(1+alpha*tmp1)
  	# N_res_new = rpois(1,N_res_new)
  	
  	# get invader's growth rate
  	N_inv_init = 1
  	tmp2 = G_inv*N_inv_init
  	N_inv_new = seedSurv*(1-G_res)*N_inv_init+Fec*tmp2/(1+alpha*tmp1)
  	lambda_inv = log(N_inv_new/N_inv_init)
  	
  	return(c(N_res_new,lambda_inv))
  }





