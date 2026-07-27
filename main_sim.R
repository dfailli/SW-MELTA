## Simulation Study ####

Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)
Sys.setenv(VECLIB_MAXIMUM_THREADS = 1)

options(mc.cores = 1)

if ("data.table" %in% installed.packages()) {
  data.table::setDTthreads(1)
}

rm(list=ls())
setwd("C:/Users/faill/OneDrive/Desktop/Mixture of experts LTA/SIMULATION")
library(MASS)
library(mclust)
library(sampling)


# Codes

load("sw.melta.RData")
load("mlta.conc.RData")
load("mlta.RData")


# Cores 

parallel::detectCores()
library(foreach)
library(parallel)
library(doParallel)
ncore = 50
registerDoParallel(cores=ncore)
getDoParWorkers()


## Setting ####

S <- 500
G <- 2
D <- 1

N <- 10000
n <- 1000
L <- 10 

P <- 2

beta0=c(0,2)

set.seed(1)
b0 = round(rbind(rnorm(L, -1.5, 1), rnorm(L,  1.5, 1)),2)

set.seed(1)
w0 = round(rnorm(L, 0.3, 0.6),2)

set.seed(1)
gamma0 = matrix(c(rep(0,L), round(rnorm(L, 0.04, 0.1),2)),P, L, byrow = T)


## Store estimates ####

b.sw.melta.nSRS <- b.sw.melta.strat <- 
  b.melta.nSRS <- b.melta.strat <- 
  b.mlta.conc.nSRS <- b.mlta.conc.strat <- 
  b.mlta.nSRS <- b.mlta.strat <- 
  b.pop <- 
  b.meltaSRS <- b.mlta.concSRS <- b.mltaSRS <- array(,c(G,L,S))

w.sw.melta.nSRS <- w.sw.melta.strat <- 
  w.melta.nSRS <- w.melta.strat <- 
  w.mlta.conc.nSRS <- w.mlta.conc.strat <-
  w.mlta.nSRS <- w.mlta.strat <- 
  w.pop <- 
  w.meltaSRS <- w.mlta.concSRS <- w.mltaSRS <- array(,c(D,L,S))

beta.sw.melta.nSRS <- beta.sw.melta.strat <- 
  beta.melta.nSRS <- beta.melta.strat <-
  beta.mlta.conc.nSRS <- beta.mlta.conc.strat <- 
  beta.pop <- 
  beta.meltaSRS <- beta.mlta.concSRS <- beta.popSRS <- array(,c(P,G-1,S))

gamma.sw.melta.nSRS <- gamma.sw.melta.strat <- 
  gamma.melta.nSRS <- gamma.melta.strat <- 
  gamma.pop <- 
  gamma.meltaSRS <- array(,c(P,L,S))

rand.sw.melta.nSRS = rand.sw.melta.strat = 
  rand.melta.nSRS = rand.melta.strat = 
  rand.mlta.conc.nSRS = rand.mlta.conc.strat = 
  rand.mlta.nSRS = rand.mlta.strat = 
  rand.pop = 
  rand.meltaSRS = rand.mlta.concSRS = rand.mltaSRS = c()

time.sw.melta.nSRS = time.sw.melta.strat = 
  time.melta.nSRS = time.melta.strat = 
  time.mlta.conc.nSRS = time.mlta.conc.strat =
  time.mlta.nSRS = time.mlta.strat = 
  time.pop = 
  time.meltaSRS = time.mlta.concSRS = time.mltaSRS = c()


## Simulation ####

block_size =  ncore
n_blocks =  ceiling(S / block_size)

all_results =  list()

for (b in 1:n_blocks) {
  
  cat("Running block", b, "of", n_blocks, "\n")
  start =  (b - 1) * block_size + 1
  end =  min(b * block_size, S)
  
  res_block =  foreach(s = start:end, .packages = c("MASS","mclust","sampling"), .errorhandling = "pass") %dopar% {
    
    set.seed(s)
    
    Xpop = cbind(rep(1,N), rnorm(N, 1, 1))
    
    exb =  exp(Xpop %*% beta0)
    eta =  cbind(1,exb)/(rowSums(exb)+1)  
    
    z =  matrix(NA, nrow=N, ncol=G)
    for(i in 1:N){
      z[i,] =  t(rmultinom(1, size = 1, prob = eta[i,]))
    }
    
    ord =  order(apply(b0, 1, mean))
    cl.pop = apply(z[,ord], 1, which.max)
    
    u =  mvrnorm(N, mu=rep(0,D), Sigma=diag(D))
    p =  matrix(NA, nrow=N, ncol=L)
    for (i in 1:N){
      for (m in 1:L){
        p[i,m] =  1/(1+exp(-(b0[which.max(z[i,]),m] + w0[m]*u[i,] + Xpop[i,] %*% gamma0[,m])))       }
    }
    
    Ypop =  matrix(NA, nrow=N, ncol=L)
    for (i in 1:N){
      for (m in 1:L){
        Ypop[i,m] =  rbinom(1, size = 1, prob = p[i,m])
      }
    }
    
    
    ### POP ####
    
    start = proc.time()
    set.seed(s)
    mod.pop <- f_MELTA_vfix(Y=Ypop, X=Xpop, G=2, D=1, maxiter=1000, npoints=7, sw=rep(1,N))
    end = proc.time()
    time.pop = (end-start)[3]
    
    b.pop =  mod.pop$b
    w.pop =  mod.pop$v
    beta.pop =  mod.pop$beta
    gamma.pop =  mod.pop$gamma
    eta.pop =  mod.pop$eta
    z.pop =  mod.pop$z
    ord.pop = order(apply(b.pop, 1, mean))
    cl.pop2 =  apply(mod.pop$z[,ord.pop],1,which.max)
    rand.pop =  adjustedRandIndex(cl.pop, cl.pop2)
    
    
    ### SRS ####
      
    set.seed(s)
    id =  sample(1:N, n)
    y =  Ypop[id, ]
    X =  Xpop[id, ]
    cl =  cl.pop[id]
    sw =  rep(1, n)
    
    
    #### SW-MELTA ####
    
    start = proc.time()
    set.seed(s)
    mod.ccs =  f_MELTA_vfix(Y=y, X=X, G=2, D=1, maxiter=1000, npoints=7, sw=sw)
    end = proc.time()
    time.ccs = (end-start)[3]
    
    b.ccs =  mod.ccs$b
    w.ccs =  mod.ccs$v
    beta.ccs =  mod.ccs$beta
    gamma.ccs =  mod.ccs$gamma
    eta.ccs =  mod.ccs$eta
    z.ccs =  mod.ccs$z
    ord.ccs = order(apply(b.ccs, 1, mean))
    cl.ccs =  apply(mod.ccs$z[,ord.ccs],1,which.max)
    rand.ccs =  adjustedRandIndex(cl, cl.ccs)
    
    
    #### MLTA with concomitant ####
    
    start = proc.time()
    set.seed(s)
    mod.mlta.conc.srs =  f_mlta_wfix.conc(X=y, DM=X, G=2, D=1, maxiter=1000, beta0=NULL, pdGH = 7)
    end = proc.time()
    time.mlta.conc.srs = (end-start)[3]
    
    b.mlta.conc.srs =  mod.mlta.conc.srs$b
    w.mlta.conc.srs =  mod.mlta.conc.srs$w
    beta.mlta.conc.srs =  mod.mlta.conc.srs$beta
    eta.mlta.conc.srs =  mod.mlta.conc.srs$eta
    z.mlta.conc.srs =  mod.mlta.conc.srs$z
    ord.mlta.conc.srs = order(apply(b.mlta.conc.srs, 1, mean))
    cl.mlta.conc.srs =  apply(mod.mlta.conc.srs$z[,ord.mlta.conc.srs],1,which.max)
    rand.mlta.conc.srs =  adjustedRandIndex(cl, cl.mlta.conc.srs)
    
    
    #### MLTA ####
    
    start = proc.time()
    set.seed(s)
    mod.mlta.srs = f_mlta_wfix(S=y, counts=rep(1,n) , G=2, D=1, maxiter=1000, pdGH = 7)
    end = proc.time()
    time.mlta.srs = (end-start)[3]
    
    b.mlta.srs =  mod.mlta.srs$b
    w.mlta.srs =  mod.mlta.srs$w
    eta.mlta.srs =  mod.mlta.srs$eta
    z.mlta.srs =  mod.mlta.srs$z
    ord.mlta.srs = order(apply(b.mlta.srs, 1, mean))
    cl.mlta.srs =  apply(mod.mlta.srs$z[,ord.mlta.srs],1,which.max)
    rand.mlta.srs =  adjustedRandIndex(cl, cl.mlta.srs)
    
    
    ### nSRS ####
      
    Ypop_str = apply(Ypop, 1, paste, collapse="")
    freq_i = as.numeric(table(Ypop_str)[Ypop_str])
    pi_i = 1 / freq_i
    pi_i = pi_i / sum(pi_i)  
    set.seed(s)
    pik = inclusionprobabilities(1/freq_i, n)
    sum(pik)
    id = which(UPsystematic(pik) == 1)
    y = Ypop[id, ]
    X = Xpop[id, ]
    cl = cl.pop[id]
    sw = 1/pik[id]
    sw = sw * N / sum(sw)
      
    
    #### SW-MELTA #### 
    
    start = proc.time()
    set.seed(s)
    mod.sw.melta.nSRS =  f_MELTA_vfix(Y=y, X=X, G=2, D=1, maxiter=1000, npoints=7, sw=sw)
    end = proc.time()
    time.sw.melta.nSRS = (end-start)[3]
    
    b.sw.melta.nSRS =  mod.sw.melta.nSRS$b
    w.sw.melta.nSRS =  mod.sw.melta.nSRS$v
    beta.sw.melta.nSRS =  mod.sw.melta.nSRS$beta
    gamma.sw.melta.nSRS =  mod.sw.melta.nSRS$gamma
    eta.sw.melta.nSRS =  mod.sw.melta.nSRS$eta
    z.sw.melta.nSRS =  mod.sw.melta.nSRS$z
    ord.sw.melta.nSRS = order(apply(b.sw.melta.nSRS, 1, mean))
    cl.sw.melta.nSRS =  apply(mod.sw.melta.nSRS$z[,ord.sw.melta.nSRS],1,which.max)
    rand.sw.melta.nSRS =  adjustedRandIndex(cl, cl.sw.melta.nSRS)
    
    
    #### MELTA ####
    
    start = proc.time()
    set.seed(s)
    mod.melta.nSRS = f_MELTA_vfix(Y=y, X=X , G=2, D=1, maxiter=1000, npoints=7, sw = rep(1,n))
    end = proc.time()
    time.melta.nSRS = (end-start)[3]
    
    b.melta.nSRS =  mod.melta.nSRS$b
    w.melta.nSRS =  mod.melta.nSRS$v
    beta.melta.nSRS =  mod.melta.nSRS$beta
    gamma.melta.nSRS =  mod.melta.nSRS$gamma
    eta.melta.nSRS =  mod.melta.nSRS$eta
    z.melta.nSRS =  mod.melta.nSRS$z
    ord.melta.nSRS = order(apply(mod.melta.nSRS$b, 1, mean))
    cl.melta.nSRS =  apply(mod.melta.nSRS$z[,ord.melta.nSRS],1,which.max)
    rand.melta.nSRS =  adjustedRandIndex(cl,cl.melta.nSRS)
    
    
    #### MLTA with concomitant ####
    
    start = proc.time()
    set.seed(s)
    mod.mlta.conc.nSRS = f_mlta_wfix.conc(X=y, DM=X, G=2, D=1, maxiter=1000, beta0=NULL, pdGH = 7)
    end = proc.time()
    time.mlta.conc.nSRS = (end-start)[3]
    
    b.mlta.conc.nSRS =  mod.mlta.conc.nSRS$b
    w.mlta.conc.nSRS =  mod.mlta.conc.nSRS$w
    beta.mlta.conc.nSRS =  mod.mlta.conc.nSRS$beta
    eta.mlta.conc.nSRS =  mod.mlta.conc.nSRS$eta
    z.mlta.conc.nSRS =  mod.mlta.conc.nSRS$z
    ord.mlta.conc.nSRS = order(apply(mod.mlta.conc.nSRS$b, 1, mean))
    cl.mlta.conc.nSRS =  apply(mod.mlta.conc.nSRS$z[,ord.mlta.conc.nSRS],1,which.max)
    rand.mlta.conc.nSRS =  adjustedRandIndex(cl,cl.mlta.conc.nSRS)
    
    
    #### MLTA ####
    
    start = proc.time()
    set.seed(s)
    mod.mlta.nSRS = f_mlta_wfix(S=y, counts=rep(1,n) , G=2, D=1, maxiter=1000, pdGH = 7)
    end = proc.time()
    time.mlta.nSRS = (end-start)[3]
    
    b.mlta.nSRS =  mod.mlta.nSRS$b
    w.mlta.nSRS =  mod.mlta.nSRS$w
    eta.mlta.nSRS =  mod.mlta.nSRS$eta
    z.mlta.nSRS =  mod.mlta.nSRS$z
    ord.mlta.nSRS = order(apply(mod.mlta.nSRS$b, 1, mean))
    cl.mlta.nSRS =  apply(mod.mlta.nSRS$z[,ord.mlta.nSRS],1,which.max)
    rand.mlta.nSRS =  adjustedRandIndex(cl,cl.mlta.nSRS)
    
    
    ### strat nSRS ####
    
    q90 <- quantile(Xpop[,2], 0.90)
    strata <- ifelse(Xpop[,2] <= q90, 1, 2)
    N1 <- sum(strata == 1)  
    N2 <- sum(strata == 2)  
    n1 <- round(n * N1/N)
    n2 <- n - n1
    Ypop_str <- apply(Ypop, 1, paste, collapse="")
    set.seed(s)
    idx1  <- which(strata == 1)
    freq1 <- as.numeric(table(Ypop_str[idx1])[Ypop_str[idx1]])
    pik1  <- inclusionprobabilities(1/freq1, n1)
    selected1 <- UPsystematic(pik1)
    id1       <- idx1[which(selected1 == 1)]
    pik_samp1 <- pik1[which(selected1 == 1)] 
    idx2  <- which(strata == 2)
    freq2 <- as.numeric(table(Ypop_str[idx2])[Ypop_str[idx2]])
    pik2  <- inclusionprobabilities(1/freq2, n2)
    selected2 <- UPsystematic(pik2)
    id2       <- idx2[which(selected2 == 1)]
    pik_samp2 <- pik2[which(selected2 == 1)]
    id  <- c(id1, id2)
    y   <- Ypop[id, ]
    X   <- Xpop[id, ]
    cl  <- cl.pop[id]
    pik_samp <- c(pik_samp1, pik_samp2)
    sw       <- 1 / pik_samp
    sw       <- sw * N / sum(sw)
    
    
    #### SW-MELTA ####
    
    start = proc.time()
    set.seed(s)
    mod.sw.melta.strat =  f_MELTA_vfix(Y=y, X=X, G=2, D=1, maxiter=1000, npoints=7, sw=sw)
    end = proc.time()
    time.sw.melta.strat = (end-start)[3]
    
    b.sw.melta.strat =  mod.sw.melta.strat$b
    w.sw.melta.strat =  mod.sw.melta.strat$v
    beta.sw.melta.strat =  mod.sw.melta.strat$beta
    gamma.sw.melta.strat =  mod.sw.melta.strat$gamma
    eta.sw.melta.strat =  mod.sw.melta.strat$eta
    z.sw.melta.strat =  mod.sw.melta.strat$z
    ord.sw.melta.strat = order(apply(b.sw.melta.strat, 1, mean))
    cl.sw.melta.strat =  apply(mod.sw.melta.strat$z[,ord.sw.melta.strat],1,which.max)
    rand.sw.melta.strat =  adjustedRandIndex(cl, cl.sw.melta.strat)
    
    
    # MELTA 
    
    start = proc.time()
    set.seed(s)
    mod.melta.strat = f_MELTA_vfix(Y=y, X=X , G=2, D=1, maxiter=1000, npoints=7, sw=rep(1,n))
    end = proc.time()
    time.melta.strat = (end-start)[3]
    
    b.melta.strat =  mod.melta.strat$b
    w.melta.strat =  mod.melta.strat$v
    beta.melta.strat =  mod.melta.strat$beta
    gamma.melta.strat =  mod.melta.strat$gamma
    eta.melta.strat =  mod.melta.strat$eta
    z.melta.strat =  mod.melta.strat$z
    ord.melta.strat = order(apply(mod.melta.strat$b, 1, mean))
    cl.melta.strat =  apply(mod.melta.strat$z[,ord.melta.strat],1,which.max)
    rand.melta.strat =  adjustedRandIndex(cl,cl.melta.strat)
    
    
    # MLTA with concomitant
    
    start = proc.time()
    set.seed(s)
    mod.mlta.conc.strat = f_mlta_wfix.conc(X=y, DM=X, G=2, D=1, maxiter=1000, beta0=NULL, pdGH = 7)
    end = proc.time()
    time.mlta.conc.strat = (end-start)[3]
    
    b.mlta.conc.strat =  mod.mlta.conc.strat$b
    w.mlta.conc.strat =  mod.mlta.conc.strat$w
    beta.mlta.conc.strat =  mod.mlta.conc.strat$beta
    eta.mlta.conc.strat =  mod.mlta.conc.strat$eta
    z.mlta.conc.strat =  mod.mlta.conc.strat$z
    ord.mlta.conc.strat = order(apply(mod.mlta.conc.strat$b, 1, mean))
    cl.mlta.conc.strat =  apply(mod.mlta.conc.strat$z[,ord.mlta.conc.strat],1,which.max)
    rand.mlta.conc.strat =  adjustedRandIndex(cl,cl.mlta.conc.strat)
    
    
    # MLTA 
    
    start = proc.time()
    set.seed(s)
    mod.mlta.strat = f_mlta_wfix(S=y, counts=rep(1,n) , G=2, D=1, maxiter=1000, pdGH = 7)
    end = proc.time()
    time.mlta.strat = (end-start)[3]
    
    b.mlta.strat =  mod.mlta.strat$b
    w.mlta.strat =  mod.mlta.strat$w
    eta.mlta.strat =  mod.mlta.strat$eta
    z.mlta.strat =  mod.mlta.strat$z
    ord.mlta.strat = order(apply(mod.mlta.strat$b, 1, mean))
    cl.mlta.strat =  apply(mod.mlta.strat$z[,ord.mlta.strat],1,which.max)
    rand.mlta.strat =  adjustedRandIndex(cl,cl.mlta.strat)
    
    list(b.pop, w.pop, beta.pop, gamma.pop, rand.pop, time.pop,
         b.ccs, w.ccs, beta.ccs, gamma.ccs, rand.ccs, time.ccs,
         b.mlta.conc.srs, w.mlta.conc.srs, beta.mlta.conc.srs, rand.mlta.conc.srs, time.mlta.conc.srs,
         b.mlta.srs, w.mlta.srs, rand.mlta.srs, time.mlta.srs,
         b.sw.melta.nSRS, w.sw.melta.nSRS, beta.sw.melta.nSRS, gamma.sw.melta.nSRS, rand.sw.melta.nSRS, time.sw.melta.nSRS,
         b.melta.nSRS, w.melta.nSRS, beta.melta.nSRS, gamma.melta.nSRS, rand.melta.nSRS, time.melta.nSRS,
         b.mlta.conc.nSRS, w.mlta.conc.nSRS, beta.mlta.conc.nSRS, rand.mlta.conc.nSRS, time.mlta.conc.nSRS,
         b.mlta.nSRS, w.mlta.nSRS, rand.mlta.nSRS, time.mlta.nSRS,
         b.sw.melta.strat, w.sw.melta.strat, beta.sw.melta.strat, gamma.sw.melta.strat, rand.sw.melta.strat, time.sw.melta.strat,
         b.melta.strat, w.melta.strat, beta.melta.strat, gamma.melta.strat, rand.melta.strat, time.melta.strat,
         b.mlta.conc.strat, w.mlta.conc.strat, beta.mlta.conc.strat, rand.mlta.conc.strat, time.mlta.conc.strat,
         b.mlta.strat, w.mlta.strat, rand.mlta.strat, time.mlta.strat)
  }
  
  all_results =  c(all_results, res_block)
  
  save(all_results, file = paste0("sim_results_up_to_", end, ".RData"))
}


## Results ####

load("C:/Users/faill/OneDrive/Desktop/Mixture of experts LTA/SIMULATION/Results/new/n_1000.RData")
load("C:/Users/faill/OneDrive/Desktop/Mixture of experts LTA/SIMULATION/Results/new/n_2000.RData")

res = all_results

for(s in 1:S){
  
  b.pop[,,s] = res[[s]][[1]]
  w.pop[,,s] = res[[s]][[2]]
  beta.pop[,,s] = res[[s]][[3]]
  gamma.pop[,,s] = res[[s]][[4]]
  rand.pop[s] = res[[s]][[5]]
  time.pop[s] = res[[s]][[6]]
  
  b.meltaSRS[,,s] = res[[s]][[7]]
  w.meltaSRS[,,s] = res[[s]][[8]]
  beta.meltaSRS[,,s] = res[[s]][[9]]
  gamma.meltaSRS[,,s] = res[[s]][[10]]
  rand.meltaSRS[s] = res[[s]][[11]]
  time.meltaSRS[s] = res[[s]][[12]]
   
  b.mlta.concSRS[,,s] = res[[s]][[13]]
  w.mlta.concSRS[,,s] = res[[s]][[14]]
  beta.mlta.concSRS[,,s] = res[[s]][[15]]
  rand.mlta.concSRS[s] = res[[s]][[16]]
  time.mlta.concSRS[s] = res[[s]][[17]]
  
  b.mltaSRS[,,s] = res[[s]][[18]]
  w.mltaSRS[,,s] = res[[s]][[19]]
  rand.mltaSRS[s] = res[[s]][[20]]
  time.mltaSRS[s] = res[[s]][[21]]
  
  b.sw.melta.nSRS[,,s] = res[[s]][[22]]
  w.sw.melta.nSRS[,,s] = res[[s]][[23]]
  beta.sw.melta.nSRS[,,s] = res[[s]][[24]]
  gamma.sw.melta.nSRS[,,s] = res[[s]][[25]]
  rand.sw.melta.nSRS[s] = res[[s]][[26]]
  time.sw.melta.nSRS[s] = res[[s]][[27]]
  
  b.melta.nSRS[,,s] = res[[s]][[28]]
  w.melta.nSRS[,,s] = res[[s]][[29]]
  beta.melta.nSRS[,,s] = res[[s]][[30]]
  gamma.melta.nSRS[,,s] = res[[s]][[31]]
  rand.melta.nSRS[s] = res[[s]][[32]]
  time.melta.nSRS[s] = res[[s]][[33]]
  
  b.mlta.conc.nSRS[,,s] = res[[s]][[34]]
  w.mlta.conc.nSRS[,,s] = res[[s]][[35]]
  beta.mlta.conc.nSRS[,,s] = res[[s]][[36]]
  rand.mlta.conc.nSRS[s] = res[[s]][[37]]
  time.mlta.conc.nSRS[s] = res[[s]][[38]]
  
  b.mlta.nSRS[,,s] = res[[s]][[39]] 
  w.mlta.nSRS[,,s] = res[[s]][[40]] 
  rand.mlta.nSRS[s] = res[[s]][[41]]
  time.mlta.nSRS[s] = res[[s]][[42]]
  
  b.sw.melta.strat[,,s] = res[[s]][[43]]
  w.sw.melta.strat[,,s] = res[[s]][[44]]
  beta.sw.melta.strat[,,s] = res[[s]][[45]]
  gamma.sw.melta.strat[,,s] = res[[s]][[46]]
  rand.sw.melta.strat[s] = res[[s]][[47]]
  time.sw.melta.strat[s] = res[[s]][[48]]
  
  b.melta.strat[,,s] = res[[s]][[49]]
  w.melta.strat[,,s] = res[[s]][[50]]
  beta.melta.strat[,,s] = res[[s]][[51]]
  gamma.melta.strat[,,s] = res[[s]][[52]]
  rand.melta.strat[s] = res[[s]][[53]]
  time.melta.strat[s] = res[[s]][[54]]
  
  b.mlta.conc.strat[,,s] = res[[s]][[55]]
  w.mlta.conc.strat[,,s] = res[[s]][[56]]
  beta.mlta.conc.strat[,,s] = res[[s]][[57]]
  rand.mlta.conc.strat[s] = res[[s]][[58]]
  time.mlta.conc.strat[s] = res[[s]][[59]]
  
  b.mlta.strat[,,s] = res[[s]][[60]]
  w.mlta.strat[,,s] = res[[s]][[61]]
  rand.mlta.strat[s] = res[[s]][[62]]
  time.mlta.strat[s] = res[[s]][[63]]
}


## Time ####

summary(time.pop)/60

summary(time.meltaSRS)/60
summary(time.mlta.concSRS)/60
summary(time.mltaSRS)/60

summary(time.sw.melta.nSRS)/60
summary(time.melta.nSRS)/60
summary(time.mlta.conc.nSRS)/60
summary(time.mlta.nSRS)/60

summary(time.sw.melta.strat)/60
summary(time.melta.strat)/60
summary(time.mlta.conc.strat)/60
summary(time.mlta.strat)/60


## Clustering ####

summary(rand.pop)

summary(rand.meltaSRS)
summary(rand.mlta.concSRS)
summary(rand.mltaSRS)

summary(sort(rand.sw.melta.nSRS)[-c(1:25)])
summary(sort(rand.melta.nSRS)[-c(1:25)])
summary(rand.mlta.conc.nSRS)
summary(rand.mlta.nSRS)

summary(sort(rand.sw.melta.strat)[-c(1:25)])
summary(sort(rand.melta.strat)[-c(1:25)])
summary(rand.mlta.conc.strat)
summary(rand.mlta.strat)


## Parameter estimates ####

reorder_b <- function(b_arr, s) {
  b <- b_arr[,,s]
  ord <- order(apply(b, 1, mean))
  b[ord, ]
}

reorder_beta <- function(b_arr, beta_arr, s) {
  b <- b_arr[,,s]
  ord <- order(apply(b, 1, mean))
  be <- cbind(0, beta_arr[,,s])
  be <- be[, ord]
  be[, 2:G] - be[, 1]
}

models_b <- list(
  ccs       = b.meltaSRS,
  mlta.concSRS = b.mlta.concSRS,
  mltaSRS = b.mltaSRS,
  sw.melta.nSRS  = b.sw.melta.nSRS,
  melta.nSRS     = b.melta.nSRS,
  mlta.conc.nSRS = b.mlta.conc.nSRS,
  mlta.nSRS      = b.mlta.nSRS,
  sw.melta.strat  = b.sw.melta.strat,
  melta.strat     = b.melta.strat,
  mlta.conc.strat = b.mlta.conc.strat,
  mlta.strat      = b.mlta.strat
)

models_w <- list(
  ccs       = w.meltaSRS,
  mlta.concSRS = w.mlta.concSRS,
  mltaSRS = w.mltaSRS,
  sw.melta.nSRS  = w.sw.melta.nSRS,
  melta.nSRS     = w.melta.nSRS,
  mlta.conc.nSRS = w.mlta.conc.nSRS,
  mlta.nSRS      = w.mlta.nSRS,
  sw.melta.strat  = w.sw.melta.strat,
  melta.strat     = w.melta.strat,
  mlta.conc.strat = w.mlta.conc.strat,
  mlta.strat      = w.mlta.strat
)

models_beta <- list(
  ccs          = list(b=b.meltaSRS, beta=beta.meltaSRS),
  mlta.concSRS = list(b=b.mlta.concSRS, beta=beta.mlta.concSRS),
  sw.melta.nSRS     = list(b=b.sw.melta.nSRS, beta=beta.sw.melta.nSRS),
  melta.nSRS        = list(b=b.melta.nSRS, beta=beta.melta.nSRS),
  mlta.conc.nSRS   = list(b=b.mlta.conc.nSRS, beta=beta.mlta.conc.nSRS),
  sw.melta.strat     = list(b=b.sw.melta.strat, beta=beta.sw.melta.strat),
  melta.strat        = list(b=b.melta.strat, beta=beta.melta.strat),
  mlta.conc.strat   = list(b=b.mlta.conc.strat, beta=beta.mlta.conc.strat)
)

models_gamma <- list(
  ccs      = gamma.meltaSRS,
  sw.melta.nSRS = gamma.sw.melta.nSRS,
  melta.nSRS    = gamma.melta.nSRS,
  sw.melta.strat = gamma.sw.melta.strat,
  melta.strat    = gamma.melta.strat
)

beta_pop_reord <- sapply(1:S, function(s) reorder_beta(b.pop, beta.pop, s))
beta_pop_reord <- array(beta_pop_reord, c(P, G-1, S))


# BIAS ####

# b
cat("b:\n")
for(nm in names(models_b)){
  all_estimates <- array(NA, dim = c(G, L, S))
  all_pop       <- array(NA, dim = c(G, L, S))
  for(s in 1:S) {
    all_estimates[,,s] <- reorder_b(models_b[[nm]], s)
    all_pop[,,s]       <- reorder_b(b.pop, s)
  }
  bias_per_cella <- apply(all_estimates - all_pop, c(1, 2), mean)
  bias_global <- mean(abs(bias_per_cella))
  cat(nm, ":", round(bias_global, 4), "\n")
}


# w
cat("\nw:\n")
for(nm in names(models_w)){
  bias_per_cella <- apply(abs(models_w[[nm]]) - abs(w.pop), c(1, 2), mean)
  bias_global    <- mean(abs(bias_per_cella))
  cat(nm, ":", round(bias_global, 4), "\n")
}


# beta
cat("\nbeta:\n")
for(nm in names(models_beta)){
  all_estimates <- array(NA, dim = c(P, G-1, S))
  for(s in 1:S) {
    all_estimates[,,s] <- reorder_beta(models_beta[[nm]]$b, models_beta[[nm]]$beta, s)
  }
  bias_per_cella <- apply(all_estimates - beta_pop_reord, c(1, 2), mean)
  bias_global    <- mean(abs(bias_per_cella))
  cat(nm, ":", round(bias_global, 4), "\n")
}


# gamma
cat("\ngamma:\n")
for(nm in names(models_gamma)){
  bias_per_cella <- apply(models_gamma[[nm]] - gamma.pop, c(1, 2), mean)
  bias_global    <- mean(abs(bias_per_cella))
  cat(nm, ":", round(bias_global, 4), "\n")
}


# relative BIAS ####

# b
cat("b:\n")
for(nm in names(models_b)){
  all_estimates <- array(NA, dim = c(G, L, S))
  all_pop       <- array(NA, dim = c(G, L, S))
  for(s in 1:S) {
    all_estimates[,,s] <- reorder_b(models_b[[nm]], s)
    all_pop[,,s]       <- reorder_b(b.pop, s)
  }
  bias_per_cella <- apply(all_estimates - all_pop, c(1, 2), mean)
  bias_global <- mean(abs(bias_per_cella/apply(all_pop, c(1,2), mean)))
  cat(nm, ":", round(bias_global, 4), "\n")
}


# w
cat("\nw:\n")
for(nm in names(models_w)){
  bias_per_cella <- apply(abs(models_w[[nm]]) - abs(w.pop), c(1, 2), mean)
  bias_global    <- mean(abs(bias_per_cella/apply(abs(w.pop), c(1,2), mean)))
  cat(nm, ":", round(bias_global, 4), "\n")
}


# beta
cat("\nbeta:\n")
for(nm in names(models_beta)){
  all_estimates <- array(NA, dim = c(P, G-1, S))
  for(s in 1:S) {
    all_estimates[,,s] <- reorder_beta(models_beta[[nm]]$b, models_beta[[nm]]$beta, s)
  }
  bias_per_cella <- apply(all_estimates - beta_pop_reord, c(1, 2), mean)
  bias_global    <- mean((abs(bias_per_cella)/apply(beta_pop_reord, c(1,2), mean)))
  cat(nm, ":", round(bias_global, 4), "\n")
}


# gamma
cat("\ngamma:\n")
for(nm in names(models_gamma)){
  bias_per_cella <- apply(models_gamma[[nm]] - gamma.pop, c(1, 2), mean)
  bias_global    <- mean((abs(bias_per_cella/apply(gamma.pop, c(1,2), mean))))
  cat(nm, ":", round(bias_global, 4), "\n")
}


# MSE ####

# b
cat("b:\n")
for(nm in names(models_b)){
  all_estimates <- array(NA, dim = c(G, L, S))
  all_pop       <- array(NA, dim = c(G, L, S))
  for(s in 1:S) {
    all_estimates[,,s] <- reorder_b(models_b[[nm]], s)
    all_pop[,,s]       <- reorder_b(b.pop, s)
  }
  mse_per_cella <- apply((all_estimates - all_pop)^2, c(1, 2), mean)
  mse_global    <- mean(mse_per_cella)
  cat(nm, ":", round(mse_global, 4), "\n")
}


# w

cat("\nw:\n")
for(nm in names(models_w)){
  mse_per_cella <- apply((abs(models_w[[nm]]) - abs(w.pop))^2, c(1, 2), mean)
  mse_global    <- mean(mse_per_cella)
  cat(nm, ":", round(mse_global, 4), "\n")
}


# beta

cat("\nbeta:\n")
for(nm in names(models_beta)){
  all_estimates <- array(NA, dim = c(P, G-1, S))
  for(s in 1:S) {
    all_estimates[,,s] <- reorder_beta(models_beta[[nm]]$b, models_beta[[nm]]$beta, s)
  }
  mse_per_cella <- apply((all_estimates - beta_pop_reord)^2, c(1, 2), mean)
  mse_global    <- mean(mse_per_cella)
  cat(nm, ":", round(mse_global, 4), "\n")
}

# gamma
cat("\ngamma:\n")
for(nm in names(models_gamma)){
  mse_per_cella <- apply((models_gamma[[nm]] - gamma.pop)^2, c(1, 2), mean)
  mse_global    <- mean(mse_per_cella)
  cat(nm, ":", round(mse_global, 4), "\n")
}


# RRMSE ####

# b
cat("b:\n")
for(nm in names(models_b)){
  all_estimates <- array(NA, dim = c(G, L, S))
  all_pop       <- array(NA, dim = c(G, L, S))
  for(s in 1:S) {
    all_estimates[,,s] <- reorder_b(models_b[[nm]], s)
    all_pop[,,s]       <- reorder_b(b.pop, s)
  }
  rmse_per_cella <- sqrt(apply((all_estimates - all_pop)^2, c(1, 2), mean))
  denom_per_cella <- apply(abs(all_estimates), c(1, 2), mean)
  rrmse_global <- mean(rmse_per_cella / denom_per_cella)
  cat(nm, ":", round(rrmse_global, 4), "\n")
}


# w
cat("\nw:\n")
for(nm in names(models_w)){
  rmse_per_cella  <- sqrt(apply((abs(models_w[[nm]]) - abs(w.pop))^2, c(1, 2), mean))
  denom_per_cella <- apply(abs(models_w[[nm]]), c(1, 2), mean)
  rrmse_global <- mean(rmse_per_cella / denom_per_cella)
  cat(nm, ":", round(rrmse_global, 4), "\n")
}


# beta
cat("\nbeta:\n")
for(nm in names(models_beta)){
  all_estimates <- array(NA, dim = c(P, G-1, S))
  for(s in 1:S) {
    all_estimates[,,s] <- reorder_beta(models_beta[[nm]]$b, models_beta[[nm]]$beta, s)
  }
  rmse_per_cella  <- sqrt(apply((all_estimates - beta_pop_reord)^2, c(1, 2), mean))
  denom_per_cella <- apply(abs(all_estimates), c(1, 2), mean)
  rrmse_global <- mean((rmse_per_cella / denom_per_cella))
  cat(nm, ":", round(rrmse_global, 4), "\n")
}


# gamma
cat("\ngamma:\n")
for(nm in names(models_gamma)){
  rmse_per_cella  <- sqrt(apply((models_gamma[[nm]] - gamma.pop)^2, c(1, 2), mean))
  denom_per_cella <- apply(abs(models_gamma[[nm]]), c(1, 2), mean)
  rrmse_global <- mean((rmse_per_cella / denom_per_cella))
  cat(nm, ":", round(rrmse_global, 4), "\n")
}
    
cat("\ngamma:\n")
for(nm in names(models_gamma)){
  rmse_per_cella  <- sqrt(apply((models_gamma[[nm]] - gamma.pop)^2, c(1, 2), mean))
  denom_per_cella <- apply(abs(models_gamma[[nm]]), c(1, 2), mean)
  rrmse_global <- mean((rmse_per_cella / denom_per_cella))
  cat(nm, ":", round(rrmse_global, 4), "\n")
}
