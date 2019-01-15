
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
    
    # prevent negative values in fecundity
    out[out<0] = 0
    
    colnames(out)=c("Fec","G1","G2")
    return(out)
}


# The population model. This tracks the population size
# of a resident, and outputs the new number of seeds,
# as well as the number of germinated plants.
# It also reports the population growth rate of the resident,
# and the population growth rate of an invader with a different germination rate.
grow_res = function(seeds_res,Fec,alpha,seedSurv,G_res,G_inv){
  	# This is an annual plant model, so it tracks populations as numbers of seeds

  	# N = population at time t
  	# Fec = fecundity
    # alpha and seedSurv are  scalars
  	# G_res = germination rates
  	# output =  list: seeds at time t+1, germinated plants at t+1
  
    #update resident
  	seedbank_carryover = seedSurv*(1-G_res)*seeds_res
  	seed_production = Fec*G_res*seeds_res/(1+alpha*G_res*seeds_res)
  	seeds_updated = seedbank_carryover + seed_production
  	# seeds_updated = rpois(1,seeds_updated)  # demographic stochasticity
  	r_res = log(seeds_updated/seeds_res)
  	
  	#update invader
  	seeds_init=1
  	seeds_inv = seedSurv*(1-G_inv)*seeds_init + Fec*G_inv*seeds_init/(1+alpha*G_res*seeds_res)
  	r_inv = log(seeds_inv/seeds_init)
  	
  	return(list("seeds" = seeds_updated, "yield" = seed_production, "r_res" = r_res, "r_inv" = r_inv))
  }
# 
# grow_inv = function(seeds_res,Fec,alpha,seedSurv,G_res,G_inv){
#   	seeds_init = 1
#   	seeds_new = seedSurv*(1-G_inv)*seeds_init+(Fec*seeds_init*G_inv)/(1+alpha*seeds_res*G_res)
#     r_inv = log(seeds_new/seeds_init)
# }


