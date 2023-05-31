## PACKAGES:
library(rjags)
library(R2WinBUGS)
library(parallel)



setwd("....") 

################################################
## DATA:

## Structure of the data to be passed on:
## 1) We have three data frames: dat.long -- longitudinal data, dat.terminal -- terminal event data, and dat.recurrent -- recurrent events data
## 2) For each data, define the following columns:
### a) For longitudinal (dat.long) 
#### --> subj.id -- Subject ID, t -- repeated measurement times for each subject, 
####     Y -- longitudinal outcome measured at each measurement time, rest is covariates -- First baseline covariates (X), then covariates for event-varying effects (Z)
##
### b) For terminal event data (dat.terminal) -- this data has 1 row per subject
#### --> subj.id -- Subject ID, T.terminal --  terminal event times for each subject, 
####     event.t -- terminal event status for each subject, rest is covariates -- First baseline covariates (X), then covariates for event-varying effects (Z)
##
### c) For recurrent event data (dat.recurrent) 
#### --> subj.id -- Subject ID, T.gap -- gap times for each subject for each recurrent event, event.r --  recurrent event status for each subject for each recurrent event 
####     r.obs.id -- index for the number of observations/rows for each subject (0,1,2,3; this is only for better visualization of the data, not used in the analysis), 
####     rest is covariates -- First baseline covariates (X), then covariates for event-varying effects (Z)

## DOWNLOAD the 3 data frames to the workspace

load("dat.long.rdata") 
load("dat.terminal.rdata")
load("dat.recurrent.rdata")

## SPECIFY: the number of baseline (X) and event-varying (Z) covariates, and the number of maximum recurrent events per subject

ncX <- 1 ## number of baseline covariates (X)
ncZ <- 1 ## number of covariates for event-varying effects (Z)
max.nr <- 3 ## maximum number of recurrent events across all subjects

BTJM_fit = function(dat.long, dat.terminal, dat.recurrent, ncX, ncZ, n.iter = 5000, n.burn = 1000, n.thin = 5, n.chains = 3, parallel = TRUE)
{
  ################################################
  ## The following script prepares data for jags
  
  source("prepare_data.R")
  
  ################################################
  ## The following part creates the txt file for the jags model
  
  load.module("glm") ## this module should be downloaded for better convergence
  
  source("jags.R")
  filename <- file.path("Trivariate_jags.txt")
  write.model(model, filename)
  
  ## Parameters -- notation follows the paper, only sigmasq.long stands for sigmasq.epsilon
  params <- c("beta0.l", "beta.l", "gamma", "phi.l", "beta.r", "beta.t", "eta.r0", "eta.r1", "zeta",
              "eta.t0", "eta.t1","alpha", "phi.t","phi.r",
              "sigmasq.long", "sigmasq.nu", "sigmasq.b0", "sigmasq.b1", "rho.01",
              "Bs.gammas.r", "Bs.gammas.t")
  
  ###############################################
  ## Run the jags model -- the following part would run the chains serially, see the code below to run the chains in parallel 
  
  if(parallel == FALSE)
  {
    model.fit <- jags.model(file = "Trivariate_jags.txt", data = Data, n.chains = n.chains, n.adapt = 1000)
    
    update(model.fit, n.burn)
    
    post.samples <- coda.samples(model.fit, params,  n.iter = n.iter, thin = n.thin)
  }
  

  
  
  ###############################################
  ## Run jags model -- the following part would run the chains in parallel
  
  if(parallel == TRUE)
  {
    coda.samples.wrapper <- function(j)
    {
      model.fit <- jags.model(file = "Trivariate_jags.txt",
                              inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = j), ## to make sure each chain starts from a different set of initial values
                              data = Data, n.chains = 1, n.adapt = 1000)
      update(model.fit, n.burn)
      
      coda.samples(model.fit, params,  n.iter = n.iter, thin = n.thin)
    }
    
    post.samples <- mclapply(1:n.chains, coda.samples.wrapper,  mc.cores = n.chains)
    
    for(ii in 1:length(post.samples))
    { post.samples[[ii]] <- post.samples[[ii]][[1]]}
    class(post.samples) <- "mcmc.list"
    
  }
  
  post.samples

}



##################
## Run BTJM_fit
##################

post.samples <- BTJM_fit(dat.long, dat.terminal, dat.recurrent, ncX, ncZ, n.iter = 5000, 
                         n.burn = 2000, n.thin = 5, n.chains = 3, parallel = TRUE)

  

##################
## Results
##################

summary(post.samples)








