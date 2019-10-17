
## Global Settings

library(psych)
library(MASS)
library(lavaan)
library(semPlot)

# Please change the working directory to your preferred folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Setting up the maximum number of warnings being presented
#options(nwarnings = 1000)


## Setting up the parameters

set.seed(020)

# Sample size, num of subjects
np <- 1000
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
pef <- .28 # The probability of an executive functioning process being activated in a domain general task
pff <- .12 # The probability of a fluid-only process being activated in a domain general task

pde <- .12 # The probability of an executive functioning process being activated in a domain specific task
pds <- .28 # The probability of a domain-specific process being activated in a domain specific task (corresponding domain)

pgs <- .10 # The probability of a process being activated in general sampling model (GSM)


## Creating Probablity Array and observation Array for both POT and GSM

# Both are 3-dimensional arrays that store info per item, task, and person.
# P stores the aggregate probability that someone gives the correct response on an item in a task;
# OBS stores the acctual responses.
P.pot <- array(rep(0, np*nt*ni), dim = c(np,nt,ni))
OBS.pot <- array(rep(0, np*nt*ni), dim = c(np,nt,ni))

P.gsm <- array(rep(0, np*nt*ni), dim = c(np,nt,ni))
OBS.gsm <- array(rep(0, np*nt*ni), dim = c(np,nt,ni))


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


## Model Specification

o.model <- '

  o =~ V1 + V2 + V3 + F1 + F2 + F3 + S1 + S2 + S3

'

g.model <- ' 

  Verbal =~ V1 + V2 + V3 
  Fluid  =~ F1 + F2 + F3
  Spatial =~ S1 + S2 + S3 
  
  g =~ Fluid + Verbal + Spatial

'

## POT

# Data generation

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

# Generating the observations (0 & 1s) according to the p values of each item
OBS.pot <- apply(P.pot, c(1,2,3), function(x) rbinom(1,1,p = x))

# Calculating the total scores of each task
Scores.pot <- apply(OBS.pot, c(1,2), sum)

# Organizing the data frame
obsdata.pot <- as.data.frame(Scores.pot)
names(obsdata.pot) <- c("F1","F2","F3","V1","V2","V3","S1","S2","S3")
cortable.pot <- cor(obsdata.pot)
cortable.pot

# Model fitting
pot.o <- sem(o.model, data = obsdata.pot, orthogonal=TRUE, fixed.x=FALSE)
pot.g <- sem(g.model, data = obsdata.pot, orthogonal=TRUE, fixed.x=FALSE)

# Model summary
summary(pot.o, fit.measures = TRUE, standardized = TRUE)
summary(pot.g, fit.measures = TRUE, standardized = TRUE)

# Plotting
semPaths(pot.o, what = "std", layout = "tree2",edge.label.cex = 1.5)
semPaths(pot.g, what = "std", layout = "tree2",edge.label.cex = 1.5)


## FCM

## Data generation

# Full Compensatory, p as a function of the average of two thetas 
P.fcm <- (1+ exp((Bg + Bs - thG - thS)/2))^(-1)
# Full Compensatory, p as an average of the functions of two thetas 
#P <- ((1+exp(Bg-thG))^(-1) + (1+exp(Bs-thS))^(-1))/2

OBS.fcm <- apply(P.fcm, c(1,2,3), function(x) rbinom(1,1,p = x))

Scores.fcm <- apply(OBS.fcm, c(1,2), sum)
obsdata.fcm <- as.data.frame(Scores.fcm)
names(obsdata.fcm) <- c("F1","F2","F3","V1","V2","V3","S1","S2","S3")
cortable.fcm <- cor(obsdata.fcm)
cortable.fcm

# Model fitting
fcm.o <- sem(o.model, data = obsdata.fcm, orthogonal=TRUE, fixed.x=FALSE)
fcm.g <- sem(g.model, data = obsdata.fcm, orthogonal=TRUE, fixed.x=FALSE)

# Model summary
summary(fcm.o, fit.measures = TRUE, standardized = TRUE)
summary(fcm.g, fit.measures = TRUE, standardized = TRUE)

# Plotting
semPaths(fcm.o, what = "std", layout = "tree2",edge.label.cex = 1.5)
semPaths(fcm.g, what = "std", layout = "tree2",edge.label.cex = 1.5)


# Creating and saving pot data at the item level
ItemResponses.pot <- cbind(OBS.pot[,1,],
                           OBS.pot[,2,],
                           OBS.pot[,3,],
                           OBS.pot[,4,],
                           OBS.pot[,5,],
                           OBS.pot[,6,],
                           OBS.pot[,7,],
                           OBS.pot[,8,],
                           OBS.pot[,9,])

pot.responses <- as.data.frame(ItemResponses.pot)

pot.var.names <- vector()

for(i in names(obsdata.pot)){
  for(j in 1:100){
    pot.var.names <- c(pot.var.names,paste0(i,".",j))
  }

}

names(pot.responses) <- pot.var.names

write.csv(pot.responses, file = "POT Item Level Responses.csv")


# Creating and saving fcm data at the item level
ItemResponses.fcm <- cbind(OBS.fcm[,1,],
                           OBS.fcm[,2,],
                           OBS.fcm[,3,],
                           OBS.fcm[,4,],
                           OBS.fcm[,5,],
                           OBS.fcm[,6,],
                           OBS.fcm[,7,],
                           OBS.fcm[,8,],
                           OBS.fcm[,9,])

fcm.responses <- as.data.frame(ItemResponses.fcm)

fcm.var.names <- vector()

for(i in names(obsdata.fcm)){
  for(j in 1:100){
    fcm.var.names <- c(fcm.var.names,paste0(i,".",j))
  }
  
}

names(fcm.responses) <- fcm.var.names

write.csv(fcm.responses, file = "FCM Item Level Responses.csv")



