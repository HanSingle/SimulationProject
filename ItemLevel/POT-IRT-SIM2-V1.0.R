library(psych)
library(MASS)
library(lavaan)
library(semPlot)


## Global Settings

# Please change the working directory to your preferred folder
setwd("C:/Users/henry/OneDrive - Claremont Graduate University/Simulation Project/IRT Simulation")

# Setting up the maximum number of warnings being presented
options(nwarnings = 1000)


## Setting up the parameters

set.seed(5)
# Sample size, num of subjects
np <- 400
# Number of tasks
nt <- 9
# Items per task
ni <- 100

# Executive Functioning processes
ne <- 50
# Fluid only processes
nf <- 50
# Domain specific processes, "v" for "Verbal" and "s" for "Spatial(Visual)"
nv <- 50
ns <- 50
# Total number of all processes
nc <- ne + nf + nv + ns

# Probabilities
pef <- .14 # The probability of an executive functioning process being activated in a domain general task
pff <- .06 # The probability of a fluid-only process being activated in a domain general task

pde <- .06 # The probability of an executive functioning process being activated in a domain specific task
pds <- .14 # The probability of a domain-specific process being activated in a domain specific task (corresponding domain)

pgs <- .05


## Creating Probablity Array and observation Array

# Both are 3-dimensional arrays that store info per item, task, and person.
# P stores the aggregate probability that someone gives the correct response on an item in a task;
# OBS stores the acctual responses.
P.pot <- array(rep(0, np*nt*ni), dim = c(np,nt,ni))
OBS.pot <- array(rep(0, np*nt*ni), dim = c(np,nt,ni))

P.fcm <- array(rep(0, np*nt*ni), dim = c(np,nt,ni))
OBS.fcm <- array(rep(0, np*nt*ni), dim = c(np,nt,ni))


## G & S are data storage arrays for thetas in POT, E is a data storage array for theta in GSM

G <- array(rep(0, np*nt*ni), dim = c(np,nt,ni))
S <- array(rep(0, np*nt*ni), dim = c(np,nt,ni))

E <- array(rep(0, np*nt*ni), dim = c(np,nt,ni))

## Creating Subjects

# Mean
mu <- rep(0, nc)
# SD
sig <- rep(1, nc)
psi <- diag(sig)
# Process Matrix (np*nc)
dat <- mvrnorm(np,mu,psi)


## Item difficulty parameters

# Creating a varing b matrix
bmu <- rep(0, ni)
bsig <- rep(1, ni)
bpsi <- diag(bsig)

vg <- as.vector(mvrnorm(nt, bmu, bpsi))# difficulty parameters for domain general processes (nt*ni elements)
vs <- as.vector(mvrnorm(nt, bmu, bpsi))# difficulty parameters for domain specific processes(nt*ni elements)

bg <- as.vector(t(replicate(np, vg)))# Replicate and reformate the nt*ni elements np times
bs <- as.vector(t(replicate(np, vs)))#

Bg <- array(bg, dim = c(np,nt,ni)) # difficulty parameters array for domain general processes (np*nt*ni)
Bs <- array(bs, dim = c(np,nt,ni)) # difficulty parameters array for domain specific processes(np*nt*ni)


## POT & FCM
## Data generation

for(id in 1:np)
{
  for(task in 1:nt)
  {
    # Choose domains and probabilities based on task
    if (task <= nt/3)
    {
      pg = pef
      ps = pff
      DSP = (ne+1):(ne+nf)
    } else if (task <= 2*nt/3)
    {
      pg = pde
      ps = pds
      DSP = (ne+nf+1):(ne+nf+nv)
    } else
    {
      pg = pde
      ps = pds
      DSP = (ne+nf+nv+1):nc
    }
    
    # Looping through every item of every task for every subject
    for(item in 1:ni)
    {
      thetaG <- dat[id, 1:ne]%*%rbinom(ne, 1, pg)
      thetaS <- dat[id, DSP]%*%rbinom(nf, 1, ps)
      G[id,task,item] <- thetaG
      S[id,task,item] <- thetaS
    }
  }
}

# Scaling the theta values of each item across subjects
thG <- apply(G, c(2,3), scale)
thS <- apply(S,c(2,3),scale)


# IRT functions, different algorisms included, uncomment corresponding line of codes to utilize

# original product of Ps, non-compensatory with varying b parameters across items
#P <- ((1+exp(Bg - thG))*(1+exp(Bs - thS)))^(-1)
# original product of Ps, non-compensatory with b = 0
#P <- ((1+exp(-thG))*(1+exp(-thS)))^(-1)
# square-rooted product of Ps, non-compensatory with varying b parameters across items
P.pot <- sqrt(((1+exp(Bg-thG))*(1+exp(Bs-thS)))^(-1))

# Compensatory, p as a function of the average of two thetas 
P.fcm <- (1+ exp((Bg + Bs - thG - thS)/2))^(-1)
# Compensatory, p as an average of the functions of two thetas 
#P <- ((1+exp(Bg-thG))^(-1) + (1+exp(Bs-thS))^(-1))/2

# Generating the observations (0 & 1s) according to the p values of each item
OBS.pot <- apply(P.pot, c(1,2,3), function(x) rbinom(1,1,p = x))
OBS.fcm <- apply(P.fcm, c(1,2,3), function(x) rbinom(1,1,p = x))

resp.pot <- matrix(c(73, 31, 100, 73, 154, 218), nrow = 2, ncol = 3,dimnames = list(c("Pre","Post"),c("SD-D","N","A-SA")))
# Organizing the data frame
obsdata <- as.data.frame(Scores)
names(obsdata) <- c("F1","F2","F3","V1","V2","V3","S1","S2","S3")


