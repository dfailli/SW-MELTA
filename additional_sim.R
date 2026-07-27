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
library(MASS)
library(mclust)


# Codes

load("sw.melta.RData")


# Cores 

parallel::detectCores()
library(foreach)
library(parallel)
library(doParallel)
ncore=50
registerDoParallel(cores=ncore)
getDoParWorkers()


## Scenarios ####

S <- 500

# benchmark
  
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
  
cluster_sep=TRUE; cluster_sep_med=FALSE; cont_plus_cat=FALSE; G3=FALSE; D2=FALSE; L20=FALSE; probit=FALSE
cluster_sep_med=TRUE; cluster_sep=FALSE; cont_plus_cat=FALSE; G3=FALSE; D2=FALSE; L20=FALSE; probit=FALSE
cont_plus_cat=TRUE; cluster_sep_med=FALSE; cluster_sep=FALSE; G3=FALSE; D2=FALSE; L20=FALSE; probit=FALSE
G3=TRUE; cluster_sep=FALSE; cluster_sep_med=FALSE; cont_plus_cat=FALSE; D2=FALSE; L20=FALSE; probit=FALSE
D2=TRUE; G3=FALSE; cluster_sep=FALSE; cluster_sep_med=FALSE; cont_plus_cat=FALSE; L20=FALSE; probit=FALSE
L20=TRUE; D2=FALSE; G3=FALSE; cluster_sep=FALSE; cluster_sep_med=FALSE; cont_plus_cat=FALSE; probit=FALSE
probit=TRUE; L20=FALSE; D2=FALSE; G3=FALSE; cluster_sep=FALSE; cluster_sep_med=FALSE; cont_plus_cat=FALSE

if (cluster_sep==TRUE){
  
  set.seed(1)
  b0 = round(rbind(rnorm(L, -0.3, 1), rnorm(L,  0.3, 1)),2)
  set.seed(1)
  w0 = round(rnorm(L, 0.1, 0.6),2)
  
} else if (cluster_sep_med==TRUE){
  
  set.seed(1)
  b0 = round(rbind(rnorm(L, -1, 1), rnorm(L,  1, 1)),2)
  
} else if (cont_plus_cat==TRUE){
  
  P=3
  beta0=c(0,2,1)
  gamma0 = matrix(c(rep(0,L), round(rnorm(L, 0.04, 0.1),2),
                    round(rnorm(L, 0.08, 0.1),2)),P, L, byrow = T)
  
} else if (G3==TRUE){
  
  G <- 3
  beta0 = rbind(c(1, -0.4), c(1.5, -0.9))
  set.seed(1)
  b0 = round(rbind(rnorm(L, -2, 0.5), 
                   rnorm(L,  0.0, 0.5), 
                   rnorm(L,  2, 0.5)), 2)
  set.seed(1)
  w0 = round(rnorm(L, 0.1, 0.6), 2)

} else if (D2==TRUE){
  
  D <- 2
  set.seed(1)
  w0 <- t(round(MASS::mvrnorm(L, mu=c(0,0), 
                              Sigma=diag(2)*(0.3)^2), 2))
  
} else if (L20==TRUE){
  
  L <- 20
  set.seed(1)
  b0_bench <- round(rbind(rnorm(10, -1.5, 1), rnorm(10, 1.5, 1)), 2)
  set.seed(2)
  b0_extra <- round(rbind(rnorm(10, -1.5, 1), rnorm(10, 1.5, 1)), 2)
  b0 <- cbind(b0_bench, b0_extra)
  
  set.seed(1)
  w0_bench <- round(rnorm(10, 0.3, 0.6), 2)
  set.seed(2)
  w0_extra <- round(rnorm(10, 0.3, 0.6), 2)
  w0 <- c(w0_bench, w0_extra)
  
  set.seed(1)
  g0_bench <- round(rnorm(10, 0.04, 0.1), 2)
  set.seed(2)
  g0_extra <- round(rnorm(10, 0.04, 0.1), 2)
  gamma0 <- matrix(c(rep(0, L), c(g0_bench, g0_extra)), P, L, byrow=TRUE)
  
} else if (probit==TRUE){
  
  b0 = b0/1.6
  w0 = w0/1.6
  gamma0 = gamma0/1.6
  
}


## Store estimates ####

b.sw.melta <- b.pop <- array(,c(G,L,S))
w.sw.melta <- w.pop <- array(,c(D,L,S))
beta.sw.melta <- beta.pop <- array(,c(P,G-1,S))
gamma.sw.melta <- gamma.pop <- array(,c(P,L,S))
eta.sw.melta <- list()
eta.pop <- array(,c(N,G,S))
rand.sw.melta = rand.pop = c()
time.sw.melta = time.pop = c()
z.sw.melta = list()
z.pop = array(,c(N,G,S))


## Simulation ####

block_size =  ncore
n_blocks =  ceiling(S / block_size)

all_results =  list()

for (b in 1:n_blocks) {
  
  cat("Running block", b, "of", n_blocks, "\n")
  start =  (b - 1) * block_size + 1
  end =  min(b * block_size, S)
  
  res_block =  foreach(s = start:end, .packages = c("MASS","mclust"), .errorhandling = "pass") %dopar% {
    
    set.seed(s)
    
    Xpop = cbind(rep(1,N), rnorm(N, 1, 1))
    
    if (cont_plus_cat==TRUE){
      Xpop = cbind(rep(1,N), rnorm(N, 1, 1), rbinom(N, 1, 0.6))
    }
    
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
    if(D2==TRUE){
      for (i in 1:N){
        for (m in 1:L){
          p[i,m] =  1/(1+exp(-(b0[which.max(z[i,]),m] + t(w0[,m])%*%u[i,] + Xpop[i,] %*% gamma0[,m])))       
        }
      }
    } else if (probit==TRUE){
      for (i in 1:N){
        for (m in 1:L){
          p[i,m] = pnorm(b0[which.max(z[i,]),m] + w0[m]*u[i,] + Xpop[i,] %*% gamma0[,m])
        }
      }
    } else {
      for (i in 1:N){
        for (m in 1:L){
          p[i,m] =  1/(1+exp(-(b0[which.max(z[i,]),m] + w0[m]*u[i,] + Xpop[i,] %*% gamma0[,m])))       
        }
      }
    }
    
    
    Ypop =  matrix(NA, nrow=N, ncol=L)
    for (i in 1:N){
      for (m in 1:L){
        Ypop[i,m] =  rbinom(1, size = 1, prob = p[i,m])
      }
    }
    
    
    # POP
    
    start = proc.time()
    set.seed(s)
    mod.pop <- f_MELTA_vfix(Y=Ypop, X=Xpop, G=G, D=D, maxiter=1000, npoints=7, sw=rep(1,N))
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
    
    
    # nSRS
    
    Ypop_str = apply(Ypop, 1, paste, collapse="")
    freq_i = as.numeric(table(Ypop_str)[Ypop_str])
    pi_i = 1 / freq_i
    pi_i = pi_i / sum(pi_i)  
    library(sampling)
    set.seed(s)
    pik = inclusionprobabilities(1/freq_i, n)
    sum(pik)
    id = which(UPsystematic(pik) == 1)
    y = Ypop[id, ]
    X = Xpop[id, ]
    cl = cl.pop[id]
    sw = 1/pik[id]
    sw = sw * N / sum(sw)
    
    
    # SW-MELTA
    
    start = proc.time()
    set.seed(s)
    mod.sw.melta =  f_MELTA_vfix(Y=y, X=X, G=G, D=D, maxiter=2000, npoints=7, sw=sw)
    end = proc.time()
    time.sw.melta = (end-start)[3]
    
    b.sw.melta =  mod.sw.melta$b
    w.sw.melta =  mod.sw.melta$v
    beta.sw.melta =  mod.sw.melta$beta
    gamma.sw.melta =  mod.sw.melta$gamma
    eta.sw.melta =  mod.sw.melta$eta
    z.sw.melta =  mod.sw.melta$z
    
    ord.sw.melta = order(apply(b.sw.melta, 1, mean))
    cl.sw.melta =  apply(mod.sw.melta$z[,ord.sw.melta],1,which.max)
    rand.sw.melta =  adjustedRandIndex(cl, cl.sw.melta)
    
    list(b.sw.melta, w.sw.melta, beta.sw.melta, gamma.sw.melta, eta.sw.melta,
         rand.sw.melta, time.sw.melta, z.sw.melta,
         b.pop, w.pop, beta.pop, gamma.pop, eta.pop, rand.pop, time.pop, z.pop)
  }
  
  all_results =  c(all_results, res_block)
  
  save(all_results, file = paste0("sim_results_up_to_", end, ".RData"))
}
