### XXX

# !!!!
# Simplest formulation of the models, e.g. without random slopes.
# These models are not used as it in the simulations nor in the sandwich terns example


# -------------------##
##### RSF  NIMBLE #####
# ------------------## 

# Bayesian RSF w/ NIMBLE
library(nimble)

rsf <- nimbleCode({
  
  ### PRIORS ###
  ## habitat cov
  beta_pop ~ dnorm(0,1)
  
  # intercept
  beta0  ~ dnorm(0,1)
  
  # fit a weighted logistic regression 
  
  for(t in 1:ntot){
    # logistic regression RSF
    logit(omega[t]) <- beta0 + beta_pop * veget[t]
    
    kase[t] ~ dbinom(omega[t],w[t])
  }
  
})


# ----------------- #
### Poisson GLM ####
# ----------------- #

#code
library(nimble)
poiglm.ni <- nimbleCode({
  
  # priors
  for(i in 1:2){
    a[i] ~ dnorm(0,1)
  }
  
  # glm 
  log(lambda[1:nsites]) <- a[1] + a[2] * (hab[1:nsites]-1) 
  
  for(i in 1:nsites){
    # likelihood
    N[i] ~ dpois(lambda[i])
  }
  
})

# ----------------------#
#### INTEGRATED MODEL ####
# ----------------------#

# with NIMBLE
intmod <- nimbleCode({
  
  # log-linear regression for the Poisson GLM
  log(lambda[1:nsites]) <- intPoi + beta_pop * (hab[1:nsites]-1)
  
  for(i in 1:nsites){
    # likelihood
    N[i] ~ dpois(lambda[i])
  }
  
  ### PRIORS ###
  ## habitat cov
  beta_pop ~ dnorm(0,1)
  # intercept
  intRSF  ~ dnorm(0,1)
  intPoi  ~ dnorm(0,1)
  
  # fit a wieghted logistic regression
  for(t in 1:ntot){
    # logistic regression RSF
    logit(omega[t]) <- intRSF + beta_pop * veget[t]
    
    kase[t] ~ dbinom(omega[t], w[t])
  }
  
})

# END