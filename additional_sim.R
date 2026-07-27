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


## Results ####

setwd("C:/Users/faill/OneDrive/Desktop/Mixture of experts LTA/SIMULATION/Results/new")
load("cluster_sep_med.RData")
load("cont_plus_cat.RData")
load("G3.RData")
load("D2.RData")
load("R20.RData")
load("probit.RData")

res = all_results

for(s in 1:S){
  
  b.sw.melta[,,s] = res[[s]][[1]]
  w.sw.melta[,,s] = res[[s]][[2]]
  beta.sw.melta[,,s] = res[[s]][[3]]
  gamma.sw.melta[,,s] = res[[s]][[4]]
  eta.sw.melta[[s]] = res[[s]][[5]]
  rand.sw.melta[s] = res[[s]][[6]]
  time.sw.melta[s] = res[[s]][[7]]
  z.sw.melta[[s]] = res[[s]][[8]]
  
  b.pop[,,s] = res[[s]][[9]]
  w.pop[,,s] = res[[s]][[10]]
  beta.pop[,,s] = res[[s]][[11]]
  gamma.pop[,,s] = res[[s]][[12]]
  eta.pop[,,s] = res[[s]][[13]]
  rand.pop[s] = res[[s]][[14]]
  time.pop[s] = res[[s]][[15]]
  z.pop[,,s] = res[[s]][[16]]
}


## Time ####

summary(time.pop)/60
summary(time.sw.melta)/60


## Clustering ####

# sample

summary(rand.pop)
summary(rand.sw.melta)


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
  sw.melta  = b.sw.melta
)

models_w <- list(
  sw.melta  = w.sw.melta
)

models_beta <- list(
  sw.melta     = list(b=b.sw.melta, beta=beta.sw.melta)
)

models_gamma <- list(
  sw.melta = gamma.sw.melta
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
  bias_global <- mean(abs(bias_per_cella)[-1,])
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
