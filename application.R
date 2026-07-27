# Data ####
rm(list=ls())
load("data.RData")
sw=mysample$PESI
sw <- sw * (length(sw) / sum(sw))
hist(sw)
sum(sw)


# Code ####

source("code_application.R")


# Model selection ####

# G=1
set.seed(1)
mod11 <- f_MELTA_g1_nstarts_vfix(Y=y, X=X, D=1, nstarts=1, tol=1e-4, maxiter=500, npoints=7, sw=sw)
bic11 = -2*mod11$LL + mod11$npar + log(nrow(y))
set.seed(1)
mod12 <- f_MELTA_g1_nstarts_vfix(Y=y, X=X, D=2, nstarts=1, tol=1e-4, maxiter=500, npoints=7, sw=sw)
bic12 = -2*mod12$LL + mod12$npar + log(nrow(y))
set.seed(1)
mod13 <- f_MELTA_g1_nstarts_vfix(Y=y, X=X, D=3, nstarts=1, tol=1e-4, maxiter=500, npoints=7, sw=sw)
bic13 = -2*mod13$LL + mod13$npar + log(nrow(y))
set.seed(1)
mod14 <- f_MELTA_g1_nstarts_vfix(Y=y, X=X, D=4, nstarts=1, tol=1e-4, maxiter=500, npoints=7, sw=sw)
bic14 = -2*mod14$LL + mod14$npar + log(nrow(y))

# Parallelizzazione per G=2,..,4

Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)
Sys.setenv(VECLIB_MAXIMUM_THREADS = 1)
options(mc.cores = 1)

rm(list=ls())
load("data.RData")
sw <- mysample$PESI
sw <- sw * (length(sw) / sum(sw))
source("code_application.R")

library(foreach)
library(parallel)
library(doParallel)
ncore <- 12
registerDoParallel(cores=ncore)

grid <- expand.grid(G=2:4, D=1:4)  

results <- foreach(i = 1:nrow(grid), .packages = c("MASS", "glmmML", "mvtnorm"),
                   .errorhandling = "pass") %dopar% {
                     
                     G <- grid$G[i]
                     D <- grid$D[i]
                     
                     set.seed(1)
                     mod <- f_MELTA_nstarts_vfix(Y=y, X=X, G=G, D=D,
                                                 nstarts=1, tol=1e-4, maxiter=1000,
                                                 npoints=7, sw=sw)
                     bic <- -2*mod$LL + mod$npar*log(nrow(y))
                     
                     list(G=G, D=D, mod=mod, bic=bic)
                   }

save(results, file="results.RData")

# Results model selection ####

load("results.RData")

grid <- expand.grid(G=2:4, D=1:4)  

for(i in 1:nrow(grid)){
  G <- grid$G[i]; D <- grid$D[i]
  assign(paste0("mod", G, D), results[[i]]$mod)
  assign(paste0("bic", G, D), results[[i]]$bic)
}

save(mod11, mod12, mod13, mod14,
     mod21, mod22, mod23, mod24,
     mod31, mod32, mod33, mod34,
     mod41, mod42, mod43, mod44, 
     file="mod.RData")

rm(list=ls())
load("data.RData")

sw=mysample$PESI
# sw <- sw * (length(sw) / sum(sw))
# hist(sw)
# sum(sw)

source("code_application.R")
load("mod.RData")

MOD=mod23

MOD$b

MOD$v

MOD$gamma

MOD$beta

P <- ncol(X)
L <- ncol(y)
N <- nrow(y)

ord=order(apply(MOD$b,1,mean))
ord
b=MOD$b[ord,]
colnames(b)=colnames(y)
b

ord=order(apply(MOD$v,1,mean))
ord
v=MOD$v[ord,]
colnames(v)=colnames(b)
v

gamma=MOD$gamma
gamma

beta=MOD$beta
beta

table(apply(MOD$z,1,which.max))/nrow(y)*100


# SE ####

Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)
Sys.setenv(VECLIB_MAXIMUM_THREADS = 1)
options(mc.cores = 1)

rm(list=ls())
load("data.RData")

sw=mysample$PESI
# sw <- sw * (length(sw) / sum(sw))
# hist(sw)
# sum(sw)

source("code_application.R")

load("mod.RData")

MOD=mod23

B = 500
n = nrow(y)
P = ncol(X)
L = ncol(y)

library(foreach)
library(parallel)
library(doParallel)
ncore=50
registerDoParallel(cores=ncore)

sw_int = pmax(1, round(sw))
id_rep = rep(1:n, sw_int)
Y_pseudo = y[id_rep, ]
X_pseudo = X[id_rep, ]
N_pseudo = nrow(Y_pseudo)

block_size =  ncore
n_blocks =  ceiling(B / block_size)

all_results =  list()

for (bl in 1:n_blocks) {
  
  cat("Running block", bl, "of", n_blocks, "\n")
  start =  (bl - 1) * block_size + 1
  end =  min(bl * block_size, B)
  
  res_block =  foreach(b = start:end, .packages = c("MASS","mclust"), .errorhandling = "pass") %dopar% {
    
    set.seed(b)
    
    id_boot = sample(1:N_pseudo, n, replace=FALSE)
    
    y_boot  = Y_pseudo[id_boot, ]
    x_boot  = X_pseudo[id_boot, ]
    
    sw_boot = rep(N_pseudo/n, n)
    sw_boot <- sw_boot * (length(sw_boot) / sum(sw_boot))
    
    mod = f_MELTA_nstarts_vfix(Y=y_boot, X=x_boot, G=2, D=3, 
                             nstarts=1, tol=1e-4, maxiter=1000,
                             npoints=7, sw=sw_boot)
    
    list(b = mod$b, 
         v = mod$v, 
         gamma = mod$gamma, 
         beta  = mod$beta)
  }
  all_results =  c(all_results, res_block)
  
  save(all_results, file = paste0("sim_results_up_to_", end, ".RData"))
}

# RESULTS SE ####

rm(list=ls())
load("Cdata.RData")
load("mod.RData")
load("se_finite_boot.RData")

N <- nrow(y)

MOD=mod23

ord=order(apply(MOD$b,1,mean))
ord
MOD$b=MOD$b[ord,]
colnames(MOD$b)=colnames(y)
item_names <- colnames(MOD$b)

ord=order(apply(MOD$v,1,mean))
ord
MOD$v=MOD$v[ord,]

B=500

res <- all_results
res <- res[sapply(res, function(x) !inherits(x, "error"))]
B_eff <- length(res)

G <- 2  
D <- 3  
P <- ncol(X)  
L <- ncol(y)  

b_boot <- array(
  unlist(lapply(res, function(x) {
    ord <- order(apply(x$b, 1, mean))
    x$b[ord, ]
  })), c(G, L, B_eff))

v_boot <- array(
  unlist(lapply(res, function(x) {
    ord <- order(apply(x$v, 1, mean))
    x$v[ord, ]
  })), c(D, L, B_eff))

beta_boot <- array(
  unlist(lapply(res, function(x) {
    ord <- order(apply(x$b, 1, mean))
    be <- cbind(rep(0, P), x$beta)
    be[, ord][ , -1]  
  })), c(P, G-1, B_eff))

gamma_boot <- array(
  unlist(lapply(res, function(x) x$gamma)),
  c(5, L, B_eff))  

SE_b     <- apply(b_boot,        c(1,2), sd)
SE_v     <- apply(abs(v_boot),   c(1,2), sd)
SE_beta  <- apply(beta_boot,     c(1,2), sd)
SE_gamma <- apply(gamma_boot,    c(1,2), sd)

z <- qnorm(0.975)

CI <- function(est, SE) list(
  lower = est - z * SE,
  upper = est + z * SE
)

b_est <- MOD$b  
ci_b  <- CI(b_est, SE_b)

v_est <- MOD$v  
ci_v  <- CI((v_est), SE_v)

beta_est <- MOD$beta  
ci_beta  <- CI(beta_est, SE_beta)

gamma_est <- MOD$gamma  
ci_gamma  <- CI(gamma_est, SE_gamma)


make_df <- function(est_row, lower_row, upper_row, set_name, 
                    item_names, skill_labels){
  data.frame(
    item     = factor(item_names, levels = item_names),
    skill    = skill_labels,
    estimate = as.vector(est_row),
    lower    = as.vector(lower_row),
    upper    = as.vector(upper_row),
    set      = set_name
  )
}


skill_labels <- c(rep("Data literacy", 4), rep("Communication", 6),
                  rep("Problem solving", 8), rep("Digital content", 7),
                  rep("Safety skills", 6))

library(dplyr)
library(ggplot2)

colori <- c(
  "Data literacy"    = "#0d0d1a",  # quasi nero
  "Communication"    = "#4d4d70",  # grigio scuro bluastro
  "Problem solving"  = "#7d7da0",  # grigio medio bluastro
  "Digital content"  = "#a8a8c5",  # grigio chiaro bluastro
  "Safety skills"    = "#d0d0e8"   # grigio chiarissimo bluastro
)

# Plot b
df_b <- bind_rows(lapply(1:nrow(b_est), function(g)
  make_df(b_est[g,], ci_b$lower[g,], ci_b$upper[g,],
          paste("Cluster", g), item_names, skill_labels)))

df_b$item <- factor(df_b$item, levels = rev(item_names))

b.plot <- ggplot(df_b, aes(x = estimate, y = item, color = skill)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  facet_wrap(~set, nrow = 1) +
  scale_color_manual(values = colori,
                     breaks = c("Data literacy", "Communication", "Problem solving", "Digital content", "Safety skills")) + 
  # scale_x_continuous(limits = c(-9, 9)) +
  labs(x = "", y = "", color = "Domain") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 8))

ggsave("b.eps", b.plot, device = cairo_ps, width = 10, height = 4)


# Plot v
df_v <- bind_rows(lapply(1:nrow(v_est), function(d)
  make_df((v_est[d,]), ci_v$lower[d,], ci_v$upper[d,],
          paste("Dimension", d), item_names, skill_labels)))

df_v$item <- factor(df_v$item, levels = rev(item_names))

v.plot <- ggplot(df_v, aes(x = estimate, y = item, color = skill)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  facet_wrap(~set, nrow = 1) +
  scale_color_manual(values = colori,
                     breaks = c("Data literacy", "Communication", "Problem solving", "Digital content", "Safety skills")) + 
  # scale_x_continuous(limits = c(-2, 2)) +
  labs(x = "", y = "", color = "Domain") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 8))

ggsave("v.eps", v.plot, device = cairo_ps, width = 10, height = 4)


# Plot gamma
set_names <- c("Intercept", "South and Island", "Female", "Age 60+", 
               "Medium/High education")
set_levels <- set_names[-1]  

df_gamma <- bind_rows(lapply(2:nrow(gamma_est), function(g)
  make_df(gamma_est[g,], ci_gamma$lower[g,], ci_gamma$upper[g,],
          set_names[g], item_names, skill_labels)))
df_gamma$set <- factor(df_gamma$set, levels = set_levels)

df_gamma$item <- factor(df_gamma$item, levels = rev(item_names))

gamma.plot <- ggplot(df_gamma, aes(x = estimate, y = item, color = skill)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  facet_wrap(~set, nrow = 1) +
  scale_color_manual(values = colori,
                     breaks = c("Data literacy", "Communication", "Problem solving", "Digital content", "Safety skills")) + 
  scale_x_continuous(limits = c(-0.45, 0.45)) +
  labs(x = "", y = "", color = "Domain") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 8))

ggsave("gamma.eps", gamma.plot, device = cairo_ps, width = 10, height = 4)


# Plot beta
item_names_beta <- c("Intercept", "South and Islands", "Female", "Age 60+", "Meidum/high education")

df_beta <- bind_rows(lapply(1:ncol(beta_est), function(g)
  make_df(beta_est[,g], ci_beta$lower[,g], ci_beta$upper[,g],
          paste("Cluster", g+1), 
          item_names_beta, rep("", length(item_names_beta)))))

df_beta

beta.plot <- ggplot(df_beta, aes(x = estimate, y = item)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  facet_wrap(~set, nrow = 1) +
  scale_x_continuous(limits = c(-15, 15)) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 8))

ggsave("beta.eps", beta.plot, device = cairo_ps, width = 10, height = 4)
