library(psych)
library(MASS)
library(lavaan)


#Change the working directory to your preferred folder
setwd("C:/Users/henry/OneDrive - Claremont Graduate University/Simulation Project")

options(nwarnings = 1000)

# Sample size
ni <- 400
# Executive Functioning processes
ne <- 50
# Fluid only processes
nf <- 50
# Domain specific processes
nv <- 50
ns <- 50

# Total
nc <- ne + nf + nv + ns


pef <- .14 # Executive functioning prob in domain general
pff <- .06 # Fluid only prob in domain general

pde <- .06 # Executive functioning in domain specific
pds <- .14 # Domain specific abilities in domain specific


# Iteration num
x <- 1000


# Creating Subjects
fm <- rep(0,8)
MG_O <- rep(list(fm),x)
MG_G <- rep(list(fm),x)

mu <- rep(0,nc)
sig <- rep(1,nc)
psi <- diag(sig)


dat <- mvrnorm(ni,mu,psi)

i = 1

while(i < x+1)
{
  tryCatch({
    Gf_scores<-dat[,1:ne]%*%replicate(3,rbinom(ne,1,p=pef)) + dat[,(ne+1):(ne+nf)]%*%replicate(3,rbinom(nf,1,p=pff));
    Gv_scores<-dat[,1:ne]%*%replicate(3,rbinom(ne,1,p=pde)) + dat[,(ne+nf+1):(ne+nf+nv)]%*%replicate(3,rbinom(nv,1,p=pds));
    Gs_scores<-dat[,1:ne]%*%replicate(3,rbinom(ne,1,p=pde)) + dat[,(ne+nf+nv+1):nc]%*%replicate(3,rbinom(ns,1,p=pds));
    scores<-cbind(Gf_scores,Gv_scores,Gs_scores);
  
    sumGf<-rowSums(dat[,1:(ne+nf)]);
    sumGv<-rowSums(dat[,(ne+nf+1):(ne+nf+nv)]);
    sumGs<-rowSums(dat[,(ne+nf+nv+1):nc]);
  
    obsdata<-as.data.frame(scores);
    alldata<-as.data.frame(cbind(scores,sumGf,sumGv,sumGs));
  
    o.model <- '
      o =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
    ';

    g.model <- ' 
      Verbal =~ V4 + V5 + V6 
      Fluid  =~ V1 + V2 + V3
      Visual =~ V7 + V8 + V9 
  
      g =~ Fluid + Verbal + Visual

    ';
  
  
    fit.o <- sem(o.model, data = alldata, orthogonal=TRUE, fixed.x=FALSE);
    fit.g <- sem(g.model, data = obsdata, orthogonal=TRUE,fixed.x=FALSE);
  
    if (fit.o@optim$converged && fit.g@optim$converged == TRUE){
      FL_O <- standardizedSolution(fit.o);
      MG_O[i] <- list(c(as.vector(fitmeasures(fit.o, c("chisq", "pvalue","cfi","rmsea","srmr"))),FL_O[c(1:9), "est.std"]));
  
      FL_G <- standardizedSolution(fit.g);
      MG_G[i] <- list(c(as.vector(fitmeasures(fit.g, c("chisq", "pvalue","cfi","rmsea","srmr"))),FL_G[c(1:12), "est.std"]));

      i = i+1;
      } else {i = i};
  }, error=function(e){})
}


# Two data frames with fit indices and factor loadings.
O <- t(as.data.frame(MG_O, row.names = c("chisq", "pvalue","cfi","rmsea","srmr","V1","V2","V3","V4","V5","V6","V7","V8","V9"), col.names = list(1:x)))
G <- t(as.data.frame(MG_G, row.names = c("chisq", "pvalue","cfi","rmsea","srmr","V1","V2","V3","V4","V5","V6","V7","V8","V9","Fluid","Verbal","Visual"),col.names = c(1:x)))

write.csv(O, file = paste("POT_O_model_", format(Sys.time(), "%H%M%S-%m-%d-%Y"), ".csv", sep = ""))
write.csv(G, file = paste("POT_G_model_", format(Sys.time(), "%H%M%S-%m-%d-%Y"), ".csv", sep = ""))


CHI.G <- describe(G[,'chisq'])

if(file.exists("Chi2 Summary of G_POT.csv")){
  write.table(CHI.G, "Chi2 Summary of G_POT.csv", append = TRUE, col.names = FALSE, sep=',')
}else{
  write.csv(CHI.G, file = "Chi2 Summary of G_POT.csv")
}

