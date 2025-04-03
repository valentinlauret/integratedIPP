# load packages 

library(tidyverse)
library(ggplot2)
library(sf)
library(amt)
library(nimble)
library(here)
library(patchwork)
library(lubridate)
library(tidyverse)
library(stars)

# defin coolors palette  https://coolors.co/5f0f40-9a031e-fb8b24-e36414-0f4c5c
blue <-  "#3d405b"
orange <- "#f2cc8f"
red <- "#e07a5f"
purple <-  "#81b29a"

# ------------------#
#### Poisson GLM #####
# ------------------#

load("tern_count.rdata")

##### Poisson GLM model ####

#code
poiglm.ni <- nimbleCode({
  
  # priors
  for(i in 1:4){
    a[i] ~ dnorm(0,1)
  }
  
  for(i in 1:2){
    b[i] ~ dnorm(0,1)
  }
  
  # latent occu
  log(lambda[1:nsites]) <- a[1] + a[2] * prof[1:nsites] +
    a[3] * prof[1:nsites] * prof[1:nsites] +
    a[4] * dcol[1:nsites]
  
  # observation process
  logit(p[1:nsites,1:nocc]) <- b[1] + b[2] * seff[1:nsites, 1:nocc]
  
  for(i in 1:nsites){
    # likelihood
    N[i] ~ dpois(lambda[i])
    
    for(j in 1:nocc){
      nobs[i,j] ~ dbin(p[i,j],N[i])
    }
  } 
  
})

# format data 
constants.o <- list(prof = scale(ypel$depth.x)[,1],
                    dcol = scale(ypel$dist_coast)[,1],
                    seff = scale(cbind(ypel$eff2017,
                                       ypel$eff2018,
                                       ypel$eff2019,
                                       ypel$eff2020,
                                       ypel$eff2021)),
                    nsites = nrow(ypel),
                    nocc = 5)

yy <- ypel %>% 
  dplyr::select(starts_with("y20")) %>% 
  st_drop_geometry()

data.o <- list(nobs = yy)

Ninit <- apply(yy,1,sum)+1
inits.o <- list(a = rnorm(4,0,1), b = rnorm(2,0,1), N = Ninit)

Rmodelo<- nimbleModel(code= nmix.ni, constants = constants.o,
                      data = data.o, inits = inits.o)

# start the nimble process
Rmodelo$initializeInfo()
Rmodelo$calculate() # - 3974

# configure model
confo <- configureMCMC(Rmodelo)

## Build and compile MCMC
Rmcmco <- buildMCMC(confo)
Cmodelo <- compileNimble(Rmodelo)
Cmcmco <- compileNimble(Rmcmco, project = Cmodelo)

# Run
resPGLM <- runMCMC(Cmcmco, niter = 110000, nburnin = 10000, nchains = 2, samplesAsCodaMCMC = TRUE)

# ------------#
#### RSF #####
# ------------#

# load data
load("tern_telemetry.rdata")

# code
rsf.terns <- nimbleCode({
  
  # PRIORS# 
  
  # at the population
  betaprof_pop ~ dnorm(0,1)
  tauprof_pop ~ dunif(0,sd = 1e2)
  betaqprof_pop ~ dnorm(0,1)
  tauqprof_pop ~ dunif(0,sd = 1e2)
  betadcol_pop ~ dnorm(0,1)
  taudcol_pop ~ dunif(0,sd = 1e2)
  
  # population intercept
  beta0_pop ~ dnorm(0,1)
  
  # individual random effect priors
  for( i in 1:nindividual){
    
    ### PRIORS ###
    beta_prof[i] ~ dnorm(betaprof_pop, sd = tauprof_pop)
    beta_qprof[i] ~ dnorm(betaqprof_pop, sd = tauqprof_pop)
    beta_dcol[i] ~ dnorm(betadcol_pop, sd = taudcol_pop)
    # intercept individual specific
    beta_0[i] ~ dnorm(beta0_pop,sd = 1e6)
  }
  
  # likelihood of a weighted logistic regression
  for(t in 1:npts){
    
    logit(omega[t]) <- beta_0[idind[t]] +
      beta_prof[idind[t]] * prof[t] + 
      beta_qprof[idind[t]] *prof [t] * prof[t] +
      beta_dcol[idind[t]] * dcol[t]
    
    kase[t] ~ dbinom(omega[t], w[t])
  }
  
})

# bundle data and run
w <- terns.data$case
w[w==0] <- 1000

nindividual = terns.data %>% pull(idind) %>% max() 
npts = terns.data %>% nrow()

# constants
constants.ni <-  list(npts = npts,
                      idind = terns.data$idind,
                      prof =  terns.data$depth.sc,
                      dcol =  terns.data$dcoast.sc,
                      nindividual = nindividual,
                      w = w)

# data
data.ni <- list(kase = terns.data$case)

# Inits
inits.ni <-  list(beta_0 = rep(0, nindividual),
                  beta_prof =  rep(0, nindividual),
                  beta_qprof=  rep(0, nindividual), 
                  beta_dcol =  rep(0, nindividual),
                  betaprof_pop = 0 ,
                  beta0_pop = 0 ,
                  tauprof_pop = runif(1,0,5),
                  betaqprof_pop = 0 ,
                  tauqprof_pop = runif(1,0,5),
                  betadcol_pop = 0 ,
                  taudcol_pop = runif(1,0,5))


# Nimble pre run
# ART: 20min on Intel i5 processor
Rmodel2 <- nimbleModel(code= rsf.terns, constants = constants.ni, data = data.ni, inits = inits.ni)
Rmodel2$initializeInfo()
Rmodel2$calculate()

# configure model
conf2 <- configureMCMC(Rmodel2)
conf2$printMonitors()
#
## Build and compile MCMC
Rmcmc2 <- buildMCMC(conf2)
Cmodel2 <- compileNimble(Rmodel2)
Cmcmc2 <- compileNimble(Rmcmc2, project = Cmodel2)


# Run
samplesRSF <- runMCMC(Cmcmc2, niter = 110000, nburnin = 10000, nchains = 2, samplesAsCodaMCMC = TRUE)  ## DT: use runMCMC

# ----------------------- #
#### INTEGRATED MODEL #####
# ----------------------- #

#  code 
# with NIMBLE package
tern.modint <- nimbleCode({
  
  # Poisson GLM
  # IPP intensity
  log(lambda[1:nsites]) <- intPoi + betaprof_pop * prof[1:nsites] +
    betaqprof_pop * prof[1:nsites] * prof[1:nsites] +
    betadcol_pop * dcol[1:nsites]
  
  # observation process
  logit(p[1:nsites,1:nocc]) <- b[1] + b[2] * seff[1:nsites, 1:nocc]
  
  
  for(i in 1:nsites){
    # likelihood
    N[i] ~ dpois(lambda[i])
    
    for(j in 1:nocc){
      nobs[i,j] ~ dbin(p[i,j],N[i])
    }
  } 
  
  # PRIORS 
  
  # Poisson GLM specific priors
  intPoi ~ dnorm(0,1)
  
  for(i in 1:2){
    b[i] ~ dnorm(0,1)
  }
  
  # shared coefficients
  betaprof_pop ~ dnorm(0,1)
  tauprof_pop ~ dunif(0,sd = 1e2)
  betaqprof_pop ~ dnorm(0,1)
  tauqprof_pop ~ dunif(0,sd=1e2)
  betadcol_pop ~ dnorm(0,1)
  taudcol_pop ~ dunif(0,sd=1e2)
  
  beta0_pop ~ dnorm(0,sd=1e6)
  
  # individual random effects
  for( i in 1:nindividual){
    
    ### PRIORS ###
    ## habitat cov
    beta_prof[i] ~ dnorm(betaprof_pop,sd = tauprof_pop)
    beta_qprof[i] ~ dnorm(betaqprof_pop,sd = tauqprof_pop)
    beta_dcol[i] ~ dnorm(betadcol_pop,sd = taudcol_pop)
    # intercept cov
    beta_0[i] ~ dnorm(beta0_pop,1e6)
  }
  
  # likelihood
  for(t in 1:npts){
    
    logit(omega[t]) <- beta_0[idind[t]] +
      beta_prof[idind[t]] * profRSF[t] + 
      beta_qprof[idind[t]] *profRSF [t] * profRSF[t] +
      beta_dcol[idind[t]] * dcolRSF[t]
    
    kase[t] ~ dbinom(omega[t], w[t])
  }
  
})

load("tern_count.rdata")
load("tern_telemetry.rdata")

# data for integrated model 
w <- terns.data$case
w[w == 0] <- 1000
nindividual = terns.data %>% pull(idind) %>% max() 
npts = terns.data %>% nrow()

constants.int <-  list(
  # RSF
  npts = npts,
  idind =   terns.data$idind,
  profRSF = terns.data$depth.sc,
  dcolRSF = terns.data$dcoast.sc,
  nindividual = nindividual,
  w = w,
  # Poisson GLM
  prof = scale(ypel$depth.x)[,1],
  dcol = scale(ypel$dist_coast)[,1],
  seff = scale(cbind(ypel$eff2017,
                     ypel$eff2018,
                     ypel$eff2019,
                     ypel$eff2020,
                     ypel$eff2021)),
  nsites = nrow(ypel),
  nocc = 5)

# data
yy <- ypel %>% 
  dplyr::select(starts_with("y20")) %>% 
  st_drop_geometry()

data.int <- list(kase = terns.data$case,
                 nobs = yy)

# Inits
Ninit <- apply(yy,1,sum)+1
inits.int <-  list(beta_0 =     rep(0, nindividual),
                   beta_prof =  rep(0, nindividual),
                   beta_qprof=  rep(0, nindividual), 
                   beta_dcol =  rep(0, nindividual),
                   betaprof_pop = 0 , tauprof_pop = runif(1,0,5),
                   betaqprof_pop = 0 , tauqprof_pop = runif(1,0,5),
                   betadcol_pop = 0 , taudcol_pop = runif(1,0,5),
                   beta0_pop = 0,
                   intPoi = rnorm(1,0,1), b = rnorm(2,0,1), N = Ninit)

# Nimble pre run
Rmodel3 <- nimbleModel(code= tern.modint, constants = constants.int,
                       data = data.int, inits = inits.int)

Rmodel3$initializeInfo()
Rmodel3$calculate()

# configure model
conf3 <- configureMCMC(Rmodel3)

## Build and compile MCMC
Rmcmc3 <- buildMCMC(conf3)
Cmodel3 <- compileNimble(Rmodel3)
Cmcmc3 <- compileNimble(Rmcmc3, project = Cmodel3)

# Run
samplesint <- runMCMC(Cmcmc3, niter = 110000, nburnin = 10000, nchains = 2, samplesAsCodaMCMC = TRUE) 

# END