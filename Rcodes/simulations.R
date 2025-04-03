### XXX

# Long piece of code to simulate landscape, tracks, transect, and fit Poisson GLM, RSF, and integrated model

##### packages ####
library(dagR)
library(tidyverse)
library(ggplot2)
library(spatialEco)
library(raster)
library(sf)
library(amt)
library(localGibbs)
library(tidyverse)
library(stars)
library(patchwork)

# define coolors palette  https://coolors.co/5f0f40-9a031e-fb8b24-e36414-0f4c5c
blue <-  "#3d405b"
orange <- "#f2cc8f"
red <- "#e07a5f"
purple <-  "#81b29a" # not purple in fact

#--------------------------------#
#### LANDSCAPE & TRANSECTS ######
# ----------------------------- #

# Simulate landscape

covlist <- list()
covlist[[1]] <- localGibbs::simRaster(rho=10, lim=c(-500,500,-500,500), res=1)
covlist[[2]] <- localGibbs::simRaster(rho=30, lim=c(-500,500,-500,500), res=1)

rr1 <- st_as_stars(covlist[[1]])
rr2 <- st_as_stars(covlist[[2]])

# discretize one covariate, scale the second.
rr1$layer <- as_factor(ifelse(rr1$layer > 0.50, 1,0))
rr2.sc <- st_as_stars(scale(covlist[[2]]))

# plot landscape
if(FALSE){
  theme_set(theme_minimal())
  
  p_habd <- ggplot() +
    geom_stars(data= rr1) + 
    scale_fill_viridis_d(name = "Discrete", labels = c("A","B"))+
    theme(legend.position = 'top')
  
  p_habc <- ggplot() +
    geom_stars(data= rr2) +
    labs(fill = "Conitnuous")+
    theme(legend.position = 'top')
  
  p_habd + p_habc
}

#------------------------#
####### SIM TRACKS #######
#-----------------------#


#biased CRW towards colony

# Parameter settings
nindividual=50  # number of individuals
nstep=300       # number of steps by individual
p.sel=4         # preference of hab2 over hab1
s0=0.8          # param 1 for biased-CRW
k=0.4           # param 2 for biased-CRW

# optionnal RSF coefficient for continuous covariate
# beta = 1.5 

# number of possible locations to prospect at each step
ntmp = 10

# pre-allocate memory
data_tot=data.frame(matrix(NA,nrow=0,ncol=10))
colnames(data_tot)=c('id','time','x','y','dist2col','step','angle','direction','veget','cov')

# optionnal: add individual heterogeneity
# coef <- rnorm(nindividual, log(p.sel), 0.25) # individual heterogeneity 

# without individual heterogeneity
coef <- rep(log(p.sel), nindividual)
sd(coef)
# note: without heterogeneity sd(coef) = 0.17 over 1000 tracks. 

# Approximate run time: 12-15 hours for 1000 individual and 300 steps on a Intel i5 processor
for(j in 1:nindividual){ # j =1
  
  kcpt= coef[j]
  
  # Create dataset for individual i
  data=data.frame(matrix(NA,nrow=nstep,ncol=10))
  colnames(data)=c('id','time','x','y','dist2col','step','angle','direction','veget','cov')
  data$id=j
  data[1,]$time=0
  data[1,]$x=0
  data[1,]$y=0
  data[1,]$dist2col=0
  
  ## First step from the colony
  
  # pre-allocate memory
  tmp <- data.frame(matrix(NA,nrow=ntmp,ncol=7))
  colnames(tmp)=c('new.coord.x','new.coord.y','veg.tmp','cov.tmp', 'p.tmp','step','angle')
  
  # choose random locations and calculate the likelihood to select each
  for(i in 1:ntmp){
    tmp$step[i]=rgamma(n=1,shape=5,scale=2)
    tmp$angle[i]=runif(n=1,0,2*pi)
    tmp[i,1:2]=anglePoint(c(0,0),angl=tmp$angle[i]+pi,len=tmp$step[i])
    tmp$veg.tmp[i] =st_extract(rr1,cbind(tmp[i,1],tmp[i,2]))
    tmp$cov.tmp[i] =st_extract(rr2.sc,cbind(tmp[i,1],tmp[i,2]))
    # calculate a likelihood to select this location
    # discrete effect only
    tmp$p.tmp[i] <-  exp(kcpt * (as.numeric(tmp$veg.tmp[i][[1]])-1) ) # for continuous covariate beta * tmp$cov.tmp[i][[1]] 
  }
  
  # choose the next location form a categorical draw in candidate locations
  loc.f <- which(rmultinom(1,1, tmp$p.tmp) ==1)[1]
  
  # store step data
  data[2,]$time = 1
  data[2,]$x=tmp$new.coord.x[loc.f]
  data[2,]$y=tmp$new.coord.y[loc.f]
  data[2,]$dist2col=sqrt(tmp$new.coord.x[loc.f]^2+tmp$new.coord.y[loc.f]^2)
  data[2,]$step=tmp$step[loc.f]
  data[2,]$angle=tmp$angle[loc.f]
  data[2,]$direction=tmp$angle[loc.f]
  data[2,]$veget=tmp$veg.tmp[loc.f][[1]]
  data[2,]$cov=tmp$cov.tmp[loc.f][[1]]
  
  ## all other steps 
  for(i in 3:nstep){
    tmp <-data.frame(matrix(NA,nrow=ntmp,ncol=9))
    colnames(tmp)=c('new.coord.x','new.coord.y','veg.tmp','cov.tmp', 'p.tmp','step','angle','direction','dist2col')
    
    for( n in 1:ntmp){
      tmp$step[n]=rgamma(n=1,shape=5,scale=2)
      alpha=atan2(data[i-2,]$y-0,data[i-2,]$x-0) - atan2(data[i-2,]$y-data[i-1,]$y,data[i-2,]$x-data[i-1,]$x)
      alpha=2*((abs(alpha)-0)/(pi-0))-1
      angle=rnorm(n=1,mean=0,sd=s0*(1+k*alpha))
      angle=ifelse(angle>pi,angle-pi,ifelse(angle<(-pi),angle+pi,angle))
      tmp$angle[n]=angle
      tmp$direction[n]=data[i-1,]$direction+angle
      tmp$direction[n]=ifelse(tmp$direction[n]>2*pi,tmp$direction[n]-2*pi,ifelse(tmp$direction[n]<0,tmp$direction[n]+2*pi,tmp$direction[n]))
      tmp[n,1:2]=anglePoint(c(data[i-1,]$x,data[i-1,]$y),angl=tmp$direction[n]+pi,len=tmp$step[n])
      tmp$dist2col[n]=sqrt(tmp[n,1]^2+tmp[n,2]^2)
      tmp$veg.tmp[n] =st_extract(rr1,cbind(tmp[n,1],tmp[n,2]))
      tmp$cov.tmp[n] =st_extract(rr2.sc,cbind(tmp[n,1],tmp[n,2]))
      # calculate a likelihood to select this location
      # discrete effect only 
      tmp$p.tmp[n] <- exp(kcpt * (as.numeric(tmp$veg.tmp[n][[1]])-1)) # beta * tmp$cov.tmp[i][[1]] +
      # continuous cov tmp$p.tmp[n] <- exp(beta * tmp$cov.tmp[n][[1]])
      
    }
    # if NAs because candidate location are outside the study area 
    # then choose any locations inside the study area
    # it does not happen very often
    if(TRUE %in% is.na(tmp$p.tmp)){
      loc.f <-  which(is.na(tmp$p.tmp) ==F)[1] 
    }else{
      loc.f <-  which(rmultinom(1,1, na.omit(tmp$p.tmp)) ==1)[1] # which(tmp$p.tmp == max(tmp$p.tmp))[1]
    }
    
    data[i,]$step=tmp$step[loc.f]
    data[i,]$angle=tmp$angle[loc.f]
    data[i,]$direction=tmp$direction[loc.f]
    data[i,]$x=tmp[loc.f,1]
    data[i,]$y=tmp[loc.f,2]
    data[i,]$dist2col=tmp$dist2col[loc.f]
    data[i,]$time=i-1
    data[i,]$veget =tmp$veg.tmp[loc.f][[1]]
    data[i,]$cov =tmp$cov.tmp[loc.f][[1]]
    
  }
  
  # bind data from individual i to previous
  data_tot=rbind(data_tot,data)
  
}

# true value to estimate 
p.est = table(data_tot$veget)[1] / table(data_tot$veget)[2]
log(p.est)

# ---------------- #
###### FIT RSF #####
# ---------------- #

head(data_tot)

# set parameters for the RSF

# how many individuals
nind <- 5
# how many points to keep among the entire track
nptInd <- 100
# how many available points to generate, should be > nptInd * 10
nptRand <- 1100

# make a loop for each individual
dfRSF <- tibble()
for(i in 1:nind){
  
  # subset data for RSF
  # this take the nind firsts individual, it does randomly select individuals.
  dat_g <- data_tot %>%
    mutate(id = as.numeric(id)) %>% 
    filter(id == i) %>% #sample(unique(data_tot$id),1)
    mutate(id = as_factor(i),
           dist = dist2col) 
  
  # select the nptInd used points
  dat_samp <- dat_g[sample(nrow(dat_g),size = nptInd),] %>% 
    dplyr::select(x,y,id,time,veget,dist) %>% 
    mutate(case = 1,
           veget = as_factor(veget -1))
  
  # generat available data 
  meanDistanceColony <- 250 # distance parameter that controls the shape of the negative exponential
  nullCoord <- data.frame(
    dist = rexp(nptRand, 1/meanDistanceColony),
    angle = runif(nptRand, 0, 2*pi)
  ) %>% 
    mutate(
      x = cos(angle)*dist,
      y = sin(angle)*dist
    )
  
  nullCoord <- nullCoord %>% 
    as_tibble() %>% 
    arrange(dist) 
  
  # remove outlier points
  nullCoord <- nullCoord[1:(0.95*nrow(nullCoord)),]
  
  rpts <- st_as_sf(nullCoord, coords= c("x","y")) %>% 
    st_crop(rr1)
  
  # extract covariate for available points
  rpts2 <- rpts %>%
    mutate(veget = as_factor(st_extract(rr1, at = rpts)$layer),
           case = 0,
           id = as_factor(i),
           time = 0,
           x = st_coordinates(rpts)[,"X"],
           y = st_coordinates(rpts)[,"Y"]) %>% 
    dplyr::select(veget,case,time,dist,x,y, id)
  
  dfRSF <- bind_rows(dfRSF, rpts2 %>% st_drop_geometry(), dat_samp)
  
}

## ----  RSF with NIMBLE ----
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

# function to run the RSF with Nimble
run.rsf <- function(dfRSF = dfRSF){
  
  w = dfRSF$case
  w[w == 0] <- 1000
  
  # constants
  constants.ni <-  list(ntot = nrow(dfRSF), 
                        veget = dfRSF$veget[],
                        w = w)
  
  # data
  data.ni <- list(kase = dfRSF$case)
  
  # Inits
  inits.ni <-  list(beta0 = 0, beta_pop = 0)
  
  Rmodel2 <- nimbleModel(code= rsf, constants = constants.ni, data = data.ni, inits = inits.ni)
  print(Rmodel2$calculate())
  
  # configure model
  conf2 <- configureMCMC(Rmodel2)
  
  ## Build and compile MCMC
  Rmcmc2 <- buildMCMC(conf2)
  Cmodel2 <- compileNimble(Rmodel2)
  Cmcmc2 <- compileNimble(Rmcmc2, project = Cmodel2)
  
  # Run
  samples <- runMCMC(Cmcmc2, niter = 110000, nburnin = 10000, nchains = 2, samplesAsCodaMCMC = TRUE)
  return(samples)
}

# to run with 50 individuals
# ART: 15min with Intel i5 processor
RSF50 <- run.rsf(dfRSF = dfRSF %>%
                   mutate(idd= as.numeric(id)) %>%
                   filter(idd <=50))

# run 10 simulations for an RSF with 20, 5, and 2 individuals
nsim = 10
RSF20 <- RSF5 <- RSF2 <- list()
for(i in 1:nsim){
  # first have to create the available points, i.e. the dfRSF, for 300 individuals
  samp_ind <- sample(1:300, 20, replace = F)
  
  RSF20[[i]] <- run.rsf(dfRSF = dfRSF %>%
                          mutate(idd= as.numeric(id)) %>%
                          filter(idd %in% samp_ind[1:20]))
  RSF5[[i]] <- run.rsf(dfRSF = dfRSF %>%
                         mutate(idd= as.numeric(id)) %>%
                         filter(idd %in% samp_ind[1:5]))
  
  RSF2[[i]] <- run.rsf(dfRSF = dfRSF %>%
                         mutate(idd= as.numeric(id)) %>%
                         filter(idd %in% samp_ind[1:2]))
}

# ----------------------#
# ##### SIM TRANSECTS ###
# --------------------- #

# TRANSECTS
# install dssd package
#devtools::install_github("DistanceDevelopment/dssd")
#remotes::install_github("DistanceDevelopment/dsims")

library(dssd)
library(dsims)

# set the central location
col <- tibble( x = 0, y = 0) %>% 
  st_as_sf(coords = c("x","y"))

# set a frame for the transects
study.area <- st_buffer(col, dist= 450, endCapStyle = "SQUARE")

region <- make.region(region.name = "study area",
                      shape = study.area)

# function to simulate transect and prepare data
simTransect <- function(trlength = trlength){
  
  # set parameters of the transects
  design <- make.design(region = region,
                        transect.type = "line",
                        design = c("eszigzagcom"),
                        line.length = trlength,
                        design.angle = c(1),
                        edge.protocol = "minus",
                        truncation = 1)
  
  # Create a single set of transects to check
  survey <- generate.transects(design)
  
  x <- survey@samplers %>%
    as.data.frame() %>%
    dplyr::select("transect","geometry") %>% 
    st_as_sf() %>%
    st_cast("LINESTRING")
  
  # make a grid around the transects
  grid <- st_make_grid(st_buffer(st_union(x), 1), cellsize = 1, square = T) %>% 
    st_as_sf() 
  
  ee <- unlist(st_intersects(st_union(x), grid))
  
  grid_c <- grid[ee,]
  
  # and extract covariates
  hab1 <- rr1 %>%
    st_as_stars() %>% 
    st_extract(at = st_centroid(grid_c), mean()) 
  
  hab2 <- rr2.sc %>%
    st_as_stars() %>% 
    st_extract(at = st_centroid(grid_c), mean() ) 
  
  # crop grid for available habitat
  grid2 <- grid_c %>%  
    mutate(hab1 = as_factor(hab1$layer),
           hab2 = hab2$layer)  %>% 
    filter(is.na(hab1) ==F)
  
  # create dist to col covariate
  grid3 <- grid2 %>% 
    mutate(distcol = unlist(purrr::map(grid2$x,st_distance,st_point(x=c(0,0))))) %>% 
    mutate(sc.dist = scale(distcol))
  
  # took only used points that intersects wih the grid
  points <- data_tot %>% 
    as_tibble() %>% 
    filter(!is.na(x)) %>% 
    st_as_sf(coords = c("x", "y"), na.fail = F) 
  
  # calculate density of used locations in each gridcell
  pt_count <-  lengths(st_intersects(grid3, points)) 
  
  # create the results table
  data.p <- grid3 %>% 
    mutate(pt_count = pt_count,
           detect = ifelse(pt_count >0,1,0))
  
  return(data.p)
}

# simulate multiple transects
data.p2k <- list()
data.p5k <- list()
data.p10k <- list()
data.p15k <- list()

# for 20 simulations
# ART: 30min
nsim = 20
for(i in 1:nsim){
  data.p10k[[i]]<- simTransect(trlength = 10000)
  data.p15k[[i]]<- simTransect(trlength = 15000)
  data.p2k[[i]]<- simTransect(trlength = 2000)
  data.p5k[[i]]<- simTransect(trlength = 5000)
  print(i)
}

# or do it once
#data.p10k<- simTransect(trlength = 10000)
#data.p15k<- simTransect(trlength = 15000)
#data.p2k <- simTransect(trlength = 2000)
#data.p5k <- simTransect(trlength = 5000)

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

# function to run it 
run.ipp <- function(data.p = data.p2k[[1]]){
  
  constants.c <- list(hab = as.numeric(data.p$hab1),
                      nsites = nrow(data.p))
  
  data.c <- list(N = data.p$pt_count)
  
  inits.c <- list(a = rnorm(2,0,1), N = data.p$pt_count)
  
  Rmodelc<- nimbleModel(code= poiglm.ni, constants = constants.c,
                        data = data.c, inits = inits.c)
  
  print(Rmodelc$calculate()) 
  # configure model
  confc <- configureMCMC(Rmodelc)
  
  ## Build and compile MCMC
  Rmcmcc <- buildMCMC(confc)
  Cmodelc <- compileNimble(Rmodelc)
  Cmcmcc <- compileNimble(Rmcmcc, project = Cmodelc)
  
  # Run
  samples2 <- runMCMC(Cmcmcc, niter = 110000, nburnin = 10000, nchains = 2, samplesAsCodaMCMC = TRUE)  ## DT: use runMCMC
  
  nimble:::clearCompiled(Cmcmcc) 
  return(samples2)
}

# fit multiple transects
# resIPP15k<- run.ipp(data.p = data.p15k)
# resIPP10k<- run.ipp(data.p = data.p10k)
# resIPP5k<- run.ipp(data.p  = data.p5k)
# resIPP2k<- run.ipp(data.p  = data.p2k)

# or for multiple simulations
resIPP5k <- resIPP2k <- resIPP15k <- resIPP10k <- list()
for(i in 1:nsim){
  resIPP15k[[i]] <- run.ipp(data.p = data.p15k[[i]])
  resIPP10k[[i]] <- run.ipp(data.p = data.p10k[[i]])
  resIPP2k[[i]] <- run.ipp(data.p = data.p2k[[i]])
  resIPP5k[[i]] <- run.ipp(data.p = data.p5k[[i]])
  print(paste0("sim n°",i))
}

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

## function to run the integrated model
runInt <- function(data.p = data.p, dfRSF = dfRSF){
  
  # weights for the RSF
  w = dfRSF$case
  w[w == 0] <- 1000
  
  constants.int <-  list(
    # RSF
    ntot = nrow(dfRSF), 
    veget = dfRSF$veget[],
    w = w,
    # Poisson GLM
    hab = as.numeric(data.p$hab1),
    nsites = nrow(data.p)) 
  
  # data
  data.int <- list(kase = dfRSF$case,
                   N = data.p$pt_count)
  
  # Inits
  inits.int <-  list(N = data.p$pt_count,
                     intPoi = 0, intRSF = 0,
                     beta_pop = 0)#, betadist = 0)
  
  
  Rmodel3 <- nimbleModel(code= intmod, constants = constants.int,
                         data = data.int, inits = inits.int)
  
  print(Rmodel3$calculate())
  # configure model
  conf3 <- configureMCMC(Rmodel3)
  
  ## Build and compile MCMC
  Rmcmc3 <- buildMCMC(conf3)
  Cmodel3 <- compileNimble(Rmodel3)
  Cmcmc3 <- compileNimble(Rmcmc3, project = Cmodel3)
  
  # Run
  samples30 <- runMCMC(Cmcmc3, niter = 110000, nburnin = 10000, nchains = 2, samplesAsCodaMCMC = TRUE) ## DT: use runMCMC
  nimble:::clearCompiled(Cmcmc3) 
  return(samples30)
}

# run the configuration you want
rsf2ipp5 <-  runInt(data.p = data.p5k, dfRSF = dfRSF %>%
                      mutate(idd = as.numeric(id)) %>% 
                      filter(idd <=2))

rsf2ipp2 <-  runInt(data.p = data.p2k, dfRSF = dfRSF %>%
                      mutate(idd = as.numeric(id)) %>% 
                      filter(idd <=2))

rsf20ipp2 <-  runInt(data.p = data.p2k, dfRSF = dfRSF %>%
                       mutate(idd = as.numeric(id)) %>% 
                       filter(idd <=20))

rsf20ipp5 <-  runInt(data.p = data.p5k, dfRSF = dfRSF %>%
                       mutate(idd = as.numeric(id)) %>% 
                       filter(idd <=20))

rsf5ipp2 <- runInt(data.p = data.p2k, dfRSF = dfRSF %>%
                     mutate(idd = as.numeric(id)) %>% 
                     filter(idd <=5))
rsf5ipp5 <- runInt(data.p = data.p5k, dfRSF = dfRSF %>%
                     mutate(idd = as.numeric(id)) %>% 
                     filter(idd <=5))
rsf5ipp10 <- runInt(data.p = data.p10k, dfRSF = dfRSF %>%
                      mutate(idd = as.numeric(id)) %>% 
                      filter(idd <=5))

# END