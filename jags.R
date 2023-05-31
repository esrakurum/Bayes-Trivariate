model<- function ()
{
  ## Cumulative sum for alpha -- used in recurrent event code
  
  for(i in 1:max.ni) ## max.ni is the maximum number of rows per subject 
  {  alpha.sum[i] <- sum(alpha1[1:i])}
  
  for (i in 1:n) 
  {
    ######
    # Longitudinal part
    ######
    b[i, 1:2] ~ dmnorm(mu0[], inv.D[, ]) ## subject-level random effects linking all three outcomes
    
    
    for(j in (sum(n.long[1:i])+1): (sum(n.long[1:(i+1)])))
    {
      mu[j] <- beta0.l + inprod(beta.l[1:ncX], x1[i, ])  + gamma * times[j]  + b[i, 1] + b[i, 2]*times[j] + inprod(phi.l[1:ncZ], z1[i, ])
      y[j] ~ dnorm(mu[j], tau.long)
    }  
    
    ######
    ## Recurrent event part
    ######
    
    nu[i] ~ dnorm(0, tau.nu) ## random effects linking survival outcomes
    
    for(j in (sum(n.r[1:i])+1): (sum(n.r[1:(i+1)]))) ## n.r[i] is the number of gap times per subject = number of events
    {
 
      log.h0.r[j] <-  inprod(Lam.r[j, ], Bs.gammas.r[]) ## baseline hazard function
      
      ## hazard function
      
      haz.r[j] <- exp(log.h0.r[j] + inprod(beta.r[1:ncX], x1[i, ]) +  inprod(phi.r1[j - sum(n.r[1:i]), ], z1[i, ]) + eta.r0*b[i, 1] +  eta.r1*b[i, 2] +  nu[i] + alpha.sum[j - sum(n.r[1:i])])
      
      ## survival function -- computed using Gauss quadrature approximation
      
      for(k in 1:K)
      {
        log.h0.rs[j, k] <- inprod(Lam.sr[K * (j - 1) + k, ], Bs.gammas.r[])
        surv.r[j, k] <- exp(log.h0.rs[j, k]  + inprod(beta.r[1:ncX], x1[i, ]) +  inprod(phi.r1[j - sum(n.r[1:i]), ], z1[i, ]) + eta.r0*b[i, 1] +  eta.r1*b[i, 2] +  nu[i] + alpha.sum[j - sum(n.r[1:i])])
        
      }
      
      log.survival.r[j] <-  -T1.r[j]/2*inprod(wk[], surv.r[j, ])
      phi.re[j] <- 100000 - (event.r[j] * log(haz.r[j])) - log.survival.r[j] 
      zeros.r[j] ~ dpois(phi.re[j])
      
    }
    ######
    ## Terminal event part 
    ######
    
    log.h0.t[i] <-  inprod(Lam.t[i, ], Bs.gammas.t[]) ## baseline hazard function
    
    ## hazard function
    haz.t[i] <- exp(log.h0.t[i] + inprod(beta.t[1:ncX], x1[i, ]) + inprod(phi.t[1:ncZ], z1[i, ]) + eta.t0*b[i, 1] +  eta.t1*b[i, 2] +  zeta*nu[i])
    
    ## survival function -- computed using Gauss quadrature approximation
    
    for(k in 1:K)
    {
      log.h0.ts[i, k] <- inprod(Lam.st[K * (i - 1) + k, ], Bs.gammas.t[])
      surv.t[i, k] <- exp(log.h0.ts[i, k]  + inprod(beta.t[1:ncX], x1[i, ]) + inprod(phi.t[1:ncZ], z1[i, ]) +  eta.t0*b[i, 1] +  eta.t1*b[i, 2] +  zeta*nu[i])
    }
    log.survival.t[i] <-  -T1.t[i]/2*inprod(wk, surv.t[i ,])
    
    phi.te[i] <- 100000 - (event.t[i] * log(haz.t[i])) - log.survival.t[i] 
    zeros.t[i] ~ dpois(phi.te[i]) 
  }
  
  ######  
  ## Priors
  ######  
  
  ## Longitudinal -- Coefficients
  
  beta0.l ~ dnorm(0, priorInvTau.beta)
  beta.l[1:ncX] ~ dmnorm(priorMean.beta[], priorInvTau.beta*identity.beta)
  
  priorInvTau.beta ~ dgamma(1, 0.005)
  
  
  phi.l[1:ncZ] ~ dmnorm(priorMean.phi[], priorInvTau.phi*identity.phi)
  
  priorInvTau.phi ~ dgamma(1, 0.005)
  
  gamma ~ dnorm(0, priorInvTau.gamma)
  priorInvTau.gamma ~ dgamma(1, 0.005)
  
  ## Recurrent -- Coefficients
  
  beta.r[1:ncX] ~ dmnorm(priorMean.beta[], priorInvTau.beta*identity.beta)
  
  alpha1[1] = 0
  for(i in 2:(max.ni))
  {  alpha1[i] ~ dnorm(0, priorInvTau.alpha)}
  
  priorInvTau.alpha ~ dgamma(1, 0.005)
  
  ## phi.r a special case, dimensions: max.ni x ncZ (maximum number of recurrent events x number of Zs )
  ## Each row is the \phi_i for all Zs, i = 1, ..., max.ni
  for(i in 1:(max.ni)) 
  {  phi.r1[i, 1:ncZ] ~ dmnorm(priorMean.phi[], priorInvTau.phi*identity.phi)}
  
  alpha[1:(max.ni-1)] <- alpha1[1:(max.ni-1)]
  phi.r[1:(max.ni-1), 1:ncZ] <- phi.r1[1:(max.ni-1), 1:ncZ] 
  
  eta.r0 ~ dnorm(0, priorInvTau.eta)
  eta.r1 ~ dnorm(0, priorInvTau.eta)
  
  priorInvTau.eta ~ dgamma(1, 0.005)
  
  
  ## Terminal -- Coefficients
  
  beta.t[1:ncX] ~ dmnorm(priorMean.beta[], priorInvTau.beta*identity.beta)
  
  phi.t[1:ncZ] ~ dmnorm(priorMean.phi[], priorInvTau.phi*identity.phi)
  
  eta.t0 ~ dnorm(0, priorInvTau.eta)
  eta.t1 ~ dnorm(0, priorInvTau.eta)
  
  zeta ~ dnorm(0, priorInvTau.zeta)
  priorInvTau.zeta ~ dgamma(1, 0.005)
  
  
  
  ## Variance and correlations
  
  sigmasq.long <- 1.0/tau.long
  tau.long ~ dgamma(1.0E-3, 1.0E-3)
  
  sigmasq.nu <- 1.0/tau.nu
  tau.nu ~ dgamma(1.0E-3, 1.0E-3)
  
  sigmasq.b0 <- 1.0/tau.b0
  tau.b0 ~ dgamma(1.0E-3, 1.0E-3)
  
  sigmasq.b1 <- 1.0/tau.b1
  tau.b1 ~ dgamma(1.0E-3, 1.0E-3)
  
  rho.01 ~ dunif(0.3, 1)
  
  Sigma.b[1,1] <- sigmasq.b0
  Sigma.b[2,2] <- sigmasq.b1
  Sigma.b[1,2] <- sqrt(sigmasq.b0)*sqrt(sigmasq.b1)*rho.01
  Sigma.b[2,1] <- Sigma.b[1,2]
  
  inv.D <- inverse(Sigma.b)
  
  
  ## Baseline hazard functions
  
  Bs.gammas.r[1:nK] ~ dmnorm(priorMean.Bs.gammas.r[], Tau.gammas.r * priorTau.gammas.r[, ])
  
  Tau.gammas.r ~ dgamma(1, 0.005)
  
  Bs.gammas.t[1:nK] ~ dmnorm(priorMean.Bs.gammas.t[], Tau.gammas.t * priorTau.gammas.t[, ])
  
  Tau.gammas.t ~ dgamma(1, 0.005)
}
