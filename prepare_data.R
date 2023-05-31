## We have three data frames: dat.long -- longitudinal data, dat.terminal -- terminal event data, and dat.recurrent -- recurrent events data

############################################
## First, gather information on the recurrent events

Time.r <- dat.recurrent$T.gap ## these are the gap times for all subjects at all recurrent event times (stacked)

event.r <- dat.recurrent$event.r ## these are the event status for all subjects at all recurrent event times (stacked) 

subjs <- dat.recurrent$subj.id ## subjects with recurrent data

n <- length(unique(subjs)) ## n: number of subjects 

N.r <- dim(dat.recurrent)[1] ## N.r: dimension of the recurrent event data

## data.id has one row per subject
data.id <- dat.recurrent[!duplicated(subjs), ]

x1 <- data.id[, 5: (length(data.id)-(ncZ))] ## secondary covariates

z1 <- data.id[,(length(data.id)-(ncZ-1)): length(data.id)] ## primary covariates 

## Number of observations (rows) for each subject in the recurrent events data, denoted as n.r
nts.r <- as.vector(table(x = subjs))

n.r <- c(0, nts.r) 

#############################################
## Second, gather information on the terminal events
 
Time.t <- dat.terminal$T.terminal ## these are the terminal event times for all subjects (stacked)
event.t <- dat.terminal$event.t ## these are the event status for all subjects (stacked) 

#############################################
## Weights and nodes for Gauss quadrature to be used in approximation of survival function -- both in terminal and recurrent events
 
wk <- c(0.0630920926299786, 0.140653259715526, 0.190350578064785,
        0.209482141084728, 0.190350578064785, 0.140653259715526,
        0.0630920926299786, 0.0229353220105292, 0.10479001032225,
        0.169004726639268, 0.204432940075299, 0.204432940075299,
        0.169004726639268, 0.10479001032225, 0.0229353220105292)
sk <- c(-0.949107912342758, -0.741531185599394, -0.405845151377397,
        0, 0.405845151377397, 0.741531185599394, 0.949107912342758,
        -0.991455371120813, -0.864864423359769, -0.586087235467691,
        -0.207784955007898, 0.207784955007898, 0.586087235467691,
        0.864864423359769, 0.991455371120813)

ordsk <- order(sk)
sk <- sk[ordsk]
wk <- wk[ordsk]

K <- length(sk) ## number of quadrature points

## For recurrent
P.r <- Time.r/2
st.r <- outer(P.r, sk + 1) ## dimension: nxlength(sk)

## For terminal
P.t <- Time.t/2
st.t <- outer(P.t, sk + 1) ## dimension:nxlength(sk)
  

############################################
## Functions for P-splines -- used in fitting the baseline hazard functions

## Function for truncated p-th power function

tpower <- function(x, t, p)
  (x - t) ^ p * (x > t)

# Function for B-spline basis
bbase1 <- function(x, xl, xr, ndx, deg){
  dx <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(D)
  B
}

nknots = 10 ## 10 knots are used in each baseline function fit
tlo = 0
nseg = nknots-2
bdeg = 2

############################################ 
## P-spline terms to estimate the baseline hazard function in the recurrent events model 

thi = floor(max(Time.r)) + 1

B0 = bbase1(Time.r, tlo, thi, nseg, bdeg)
Lam.r <- B0 ## dimension: nxnknots, one row for each subject -- this is for the hazard function


B0 = bbase1(c(t(st.r)), tlo, thi, nseg, bdeg)
Lam.sr <- B0 

DDal <- diag(ncol(Lam.r)) 
priorTau.gammas.r <- crossprod(diff(DDal, diff = 2)) + 1e-04 * DDal

############################################
## P-spline terms to estimate the baseline hazard function in the terminal event model 

thi = floor(max(Time.t)) + 1

B0 = bbase1(Time.t, tlo, thi, nseg, bdeg)
Lam.t <- B0 


B0 = bbase1(c(t(st.t)), tlo, thi, nseg, bdeg)
Lam.st <- B0 

DDal <- diag(ncol(Lam.t)) 
priorTau.gammas.t <- crossprod(diff(DDal, diff = 2)) + 1e-04 * DDal

############################################
## Lastly, gather information on the longitudinal outcome

nts.long <- table(x = dat.long$subj_id) ## Number of repeated measurements per subject

n.long <- c(0, as.vector(nts.long)) 

times <- dat.long$t ## repeated measurement times
############################################
## Mean for the priors

Bs.gammas.r <- rep(0, nknots)
Bs.gammas.t <- rep(0, nknots)

beta.mean <- rep(0, ncX)
phi.mean <- rep(0, ncZ)

############################################
# Data that will be passed on to jags

Data <- list(n = n, 
             K = K, 
             wk = wk, ## weights used in approximating survival function
             y = dat.long$Y, ## longitudinal outcome
             times = times, 
             n.long = n.long,
             mu0 = c(0, 0), ## prior mean for subject-level random effects
             x1 = matrix(x1, ncol = ncX),
             z1 = matrix(z1, ncol = ncX),
             ncX = ncX, 
             ncZ = ncZ,
             n.r = n.r, 
             max.ni = max(n.r),
             event.r = event.r,
             event.t = event.t,
             T1.r = Time.r,
             T1.t = Time.t,
             zeros.r = rep(0, N.r), 
             zeros.t = rep(0, n),
             priorMean.beta = beta.mean,
             priorMean.phi = phi.mean,
             Lam.sr = Lam.sr,
             Lam.r = Lam.r,
             priorMean.Bs.gammas.r = Bs.gammas.r,
             priorTau.gammas.r = priorTau.gammas.r,
             Lam.st = Lam.st,
             Lam.t = Lam.t,
             priorMean.Bs.gammas.t = Bs.gammas.t,
             priorTau.gammas.t = priorTau.gammas.t,
             nK = nknots,
             identity.beta = diag(1, ncX, ncX),
             identity.phi = diag(1, ncZ, ncZ)
             )


