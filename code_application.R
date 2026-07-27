library(MASS)
library(foreach)
library(mclust)

## Model estimation ####

f_MELTA_vfix <- function(Y, X, G, D, tol, maxiter, npoints, sw)
{
  
  # Function parameters:
  # Y = matrix of response variables
  # X = covariates
  # G = number of mixture components
  # D = latent trait dimension
  # tol = tolerance level for convergence
  # maxiter = maximum number of iterations
  # npoints = quadrature points
  # sw = sampling weights
  
  
  N <- nrow(Y)         # units
  L <- ncol(Y)         # items
  P <- ncol(X)         # covariates (including intercept)
  
  # Remove intercept from X for gamma (identifiability: gamma_0k = 0)
  X_gamma <- as.matrix(X[, -1, drop=FALSE])
  P_gamma <- ncol(X_gamma)  # P - 1
  dim_wh  <- D + G + P_gamma
  
  
  # Parameter initialization
  
  beta <- rep(0,P*(G-1))
  beta <- matrix(beta,ncol=G-1)
  # beta=beta0
  
  gamma <- matrix(rep(0, P_gamma*L), nrow=P_gamma, ncol=L)
  
  exb <- exp(X %*% beta)
  eta <- cbind(1,exb)/(rowSums(exb)+1)  # priors
  
  z <- matrix(NA, nrow=N, ncol=G)
  for(i in 1:N)   z[i,] <- t(rmultinom(1, size = 1, prob = eta[i,]))
  
  v <- matrix(rnorm(L * D), D, L)
  
  b=matrix(rnorm(L*G,0, 0.000001), G, L) # matrix(rnorm(G*L), G, L)
  ord.b=order(apply(b,1,mean))
  b <- b[ord.b,]
  
  
  # Variational approximation initialization
  
  xi <- array(20, c(N, L, G))
  sigma_xi <- 1 / (1 + exp(-xi))
  lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
  
  C <- array(0, c(D*D, N, G))
  mu <- array(0, c(N, D, G))
  
  
  # Iterative process
  
  ll <- -Inf
  diff <- 1
  iter <- 0
  cond <- TRUE
  
  print(c(beta))
  
  lxi <- matrix(0, N, G) # variational approx
  
  YY <- array(0, c(D * D, N, G))
  gam <- matrix(0, dim_wh, L)
  K <- array(diag(dim_wh), c(dim_wh, dim_wh, L))
  wh <- matrix(0, dim_wh, L) # matrix to store model parameters
  
  while (diff > tol & iter < maxiter)
  {
    iter <- iter + 1
    ll_old <- ll
    
    ## D=1 ####
    if (D == 1) {
      for (g in 1:G)
      {
        # Posterior Statistics
        
        C[, , g] <-
          1 / (1 - 2 * rowSums(sweep(lambda_xi[, , g], MARGIN = 2, v ^ 2, `*`)))
        mu[, , g] <-
          C[, , g] * rowSums(sweep(
            Y - 0.5 + 2 * lambda_xi[, , g] * (X_gamma%*%gamma) + 
              2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`),
            MARGIN = 2, v, `*`))
        
        YY[, , g] <- matrix(C[, , g] + mu[, , g] ^ 2, ncol = 1)
        
        # Variational Parameters (xi)
        
        xi[, , g] <-  
          (YY[, , g] %*% v^2) + (mu[, , g] %*% (2 * b[g, ] * v)) +
          (apply((t(apply(2 * (X_gamma%*%gamma),1,function(x) x* v))), 2, function(y) y*mu[,,g])) + 
          (2 * t(apply(X_gamma%*%gamma,1,function(xx)  xx * b[g, ]))) + 
          (matrix(b[g, ]^2, nrow = N, ncol = L, byrow = TRUE)) + 
          (matrix((X_gamma%*%gamma)^2, nrow = N, ncol = L, byrow = FALSE))
      }
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      # Model Parameters
      
      gam[1, ] <- colSums(crossprod(sw * z * mu[, 1, ], Y - 0.5)) # v
      gam[2:(1 + G), ] <- crossprod(sw * z, Y - 0.5) # b
      gam[(D+G+1):dim_wh, ] <- t(apply(apply(array(apply(sw * z,2,function(x) x*X_gamma),c(N,P_gamma,G)), # gamma
                                             c(2,3), function(y) crossprod(y, Y-0.5)),
                                       c(1,2), function(xy) sum(xy)))
      
      
      K[1, 1, ] <- apply(aperm(array(sw * z * YY[1, , ], c(N, G, L)), c(1, 3, 2)) * lambda_xi, 2, sum) # vv
      
      K[2:(G + 1), 1, ] <- t(apply(aperm(array( # b v
        sw * z * mu[, 1, ], c(N, G, L)
      ), c(1, 3, 2)) * lambda_xi, c(2, 3), sum))
      
      K[1, 2:(G + 1), ] <- K[2:(G + 1), 1, ] # v b
      
      K[(D+G+1):dim_wh, 1, ] <- apply(apply(array(apply(lambda_xi,2, function(x) x*(sw * z * mu[, 1, ])), # gamma v
                                                  c(N,G,L)), c(2,3), function(xx) crossprod(xx, X_gamma)),
                                      c(1,3),sum)
      
      K[1, (D+G+1):dim_wh, ] <- K[(D+G+1):dim_wh, 1, ] # v gamma
      
      for(l in 1:L){
        
        diag(K[2:(G + 1),2:(G + 1),l]) <- apply(array(apply(lambda_xi,2,function(x) sw*z*x),c(N,G,L)), # b b
                                                c(2,3), sum)[,l]
        
        K[2:(G + 1), (D+G+1):dim_wh, l] <- t(apply(array(apply(lambda_xi,2,function(x) sw*z*x),c(N,G,L)), # b gamma
                                                   c(2,3), function(xx) crossprod(xx, X_gamma))[,,l])
        
        K[(D+G+1):dim_wh, 2:(G + 1), l] = t(K[2:(G + 1), (D+G+1):dim_wh, l]) # gamma b
        
        valori_gg <- apply(array(apply(array(apply(array(apply(lambda_xi,2,function(x) sw*z*x),c(N,G,L)),
                                                   c(2,3), function(xx) xx*X_gamma),c(N,P_gamma,G,L)),
                                       c(3,4), function(xxx) xxx*X_gamma),c(N,P_gamma,G,L)), # gamma gamma
                           c(2,4), sum)[,l]
        if(P_gamma == 1){
          K[(D+G+1):dim_wh, (D+G+1):dim_wh, l] <- valori_gg
        } else {
          for(jj in seq_len(P_gamma)){
            K[(D+G+jj), (D+G+jj), l] <- valori_gg[jj]
          }
        }
        
        wh[, l] <- -ginv(2 * K[, , l]) %*% gam[, l] # parameter updates
      }
      
      
      v <- wh[1:D, ]
      v <- as.matrix(t(v))
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
      }
      ord.b=order(apply(b,1,mean))
      b <- b[ord.b,]
      
      gamma=wh[(D+G+1):dim_wh,]
      
      
      # log(P(Y|Z,X,xi))
      
      for(g in 1:G){
        lxi[,g] <-
          (0.5 * log(C[1, , g])) + (mu[, 1, g] ^ 2 / (2 * C[1, , g])) + rowSums(
            (log(sigma_xi[, , g])) - (0.5 * xi[, , g]) - (lambda_xi[, , g] * xi[, , g] ^ 2) + 
              (sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`)) +
              (lambda_xi[, , g] * (X_gamma%*%gamma)^2) +
              (2 * sweep(lambda_xi[, , g] * (X_gamma%*%gamma), MARGIN=2, b[g, ], `*`)) + 
              (sweep(Y - 0.5, MARGIN = 2, b[g, ], `*`)) + 
              ((Y-0.5) * (X_gamma%*%gamma)))
      } 
    } # end D=1
    
    ## D=2 ####
    if (D == 2) {
      for (g in 1:G)
      {
        # Latent Posterior Statistics
        
        C[, , g] <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(v), t(v))))
        
        mu[, , g] <-
          t(apply(rbind(C[, , g], 
                        tcrossprod(v, Y - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`)+
                                     2 * (lambda_xi[, , g] * (X_gamma%*%gamma)))), 
                  2, function(x) matrix(x[1:4], nrow = D) %*% x[-(1:4)]))
        
        # Variational Parameters (xi)
        
        YY[, , g] <- C[, , g] + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY[, , g], 2, function(x)
              rowSums(crossprod(v, matrix(x, ncol = D)) * t(v))) + 
              tcrossprod(2 * b[g, ] * t(v), mu[, , g]) + 
              2 * t(X_gamma%*%gamma) * tcrossprod(t(v),mu[,,g]) + 
              2 * t(X_gamma%*%gamma) * b[g,] + 
              + matrix(
                b[g, ] ^ 2,
                nrow = L,
                ncol = N,
                byrow = FALSE
              ) + 
              t((X_gamma%*%gamma)^2)
          )
      }
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      
      # Model Parameters 
      
      gam[3:(2 + G), ] <- crossprod(sw*z, Y - 0.5) # b
      
      gam[(D+G+1):dim_wh, ] <- t(apply(apply(array(apply(sw*z,2,function(x) x*X_gamma),c(N,P_gamma,G)), # gamma
                                             c(2,3), function(y) crossprod(y, Y-0.5)),
                                       c(1,2), function(xy) sum(xy)))
      
      aa <- aperm(array(sw*z, c(N, G, D)), c(1, 3, 2)) * mu
      bb <- aperm(array(sw*z, c(N, G, D * D)), c(3, 1, 2)) * YY
      
      K[3:(2 + G), 3:(2 + G), ] <- # b b
        apply(apply(aperm(array(sw*z, c(
          N, G, L
        )), c(1, 3, 2)) * lambda_xi, c(2, 3), sum), 1, diag)
      
      kk <- 0
      vgamma <- 0
      for (g in 1:G) {
        
        K[2 + g, 1:2, ] <- crossprod(aa[, , g], lambda_xi[, , g]) # b v
        K[1:2, 2 + g, ] <- K[2 + g, (1:2), ] # v b
        
        kk <- kk + bb[, , g] %*% lambda_xi[, , g]
        
        vgamma <- vgamma + apply(array(apply(array(apply(aa[,,g],2,function(x) x*lambda_xi[,,g]),c(N,L,D)),
                                             c(2,3), function(xx) xx*X_gamma),c(N,P_gamma,L,D)),
                                 c(2,3,4), sum)
      }
      
      K[1:2, 1:2, ] <- kk # v v
      
      for (l in 1:L)
      {
        gam[1:D, l] <- apply(aa * (Y[, l] - 0.5), 2, sum) # v
        
        K[(D+G+1):dim_wh, 1:D, l] = vgamma[,l,] # gamma v
        K[1:D, (D+G+1):dim_wh, l] = t(K[(D+G+1):dim_wh, 1:D, l])  # v gamma 
        
        K[(D+1):(G + D), (D+G+1):dim_wh, l] <- t(apply(array(apply(lambda_xi,2,function(x) sw*z*x),c(N,G,L)), # b gamma
                                                       c(2,3), function(xx) crossprod(xx, X_gamma))[,,l])
        K[(D+G+1):dim_wh, (D+1):(G + D), l] = t(K[(D+1):(G + D), (D+G+1):dim_wh, l]) # gamma b
        
        valori_gg <- apply(array(apply(array(apply(array(apply(lambda_xi,2,function(x) sw*z*x),c(N,G,L)),
                                                   c(2,3), function(xx) xx*X_gamma),c(N,P_gamma,G,L)),
                                       c(3,4), function(xxx) xxx*X_gamma),c(N,P_gamma,G,L)), # gamma gamma
                           c(2,4), sum)[,l]
        if(P_gamma == 1){
          K[(D+G+1):dim_wh, (D+G+1):dim_wh, l] <- valori_gg
        } else {
          for(jj in seq_len(P_gamma)){
            K[(D+G+jj), (D+G+jj), l] <- valori_gg[jj]
          }
        }
        
        wh[, l] <- -ginv(2 * K[, , l]) %*% gam[, l]
      }
      
      v <- wh[1:D, ]
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
      }
      ord.b=order(apply(b,1,mean))
      b <- b[ord.b,]
      
      gamma=wh[(D+G+1):dim_wh,]
      
      
      # log(P(Y|Z,X,xi))
      
      for(g in 1:G){
        
        detC <- C[1, , g] * C[4, , g] - C[3, , g] * C[2, , g]
        
        lxi[,g] <-
          (0.5 * log(detC)) + (0.5 * apply(rbind(C[4, , g] / detC, -C[2, , g] / detC, -C[3, , g] /
                                                   detC, C[1, , g] / detC, t(mu[, , g])), 2, function(x)
                                                     t((x[-(1:4)])) %*% matrix(x[1:4], nrow = D) %*% (x[-(1:4)]))) + 
          rowSums(
            (log(sigma_xi[, , g])) - (0.5 * xi[, , g]) - (lambda_xi[, , g] * xi[, , g]^2) + 
              (sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`)) + 
              (lambda_xi[, , g] * (X_gamma%*%gamma)^2) +
              (2 * sweep(lambda_xi[, , g] * (X_gamma%*%gamma), MARGIN=2, b[g, ], `*`)) + 
              (sweep(Y - 0.5, MARGIN = 2, b[g, ], `*`)) +
              ((Y-0.5) * (X_gamma%*%gamma))
          )
      } 
    }# end D=2
    
    ## D>2 ####
    if (D > 2) {
      for (g in 1:G)
      {
        # Latent Posterior Statistics
        
        C[, , g] <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(v), t(v))))
        
        mu[, , g] <-
          t(apply(rbind(C[, , g], v %*% t(
            Y - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`) +
              2 * lambda_xi[, , g] * (X_gamma%*%gamma)
          )), 2, function(x)
            matrix(x[1:(D * D)], nrow = D) %*% x[-(1:(D * D))]))
        
        # Variational Parameters (xi)
        
        YY[, , g] <- C[, , g] + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY[, , g], 2, function(x)
              rowSums((t(v) %*% matrix(x, ncol = D)) * t(v))) + 
              (2 * b[g, ] * t(v)) %*% t(mu[, , g]) + 
              2 * t(X_gamma%*%gamma) * tcrossprod(t(v),mu[,,g]) + 
              2 * t(X_gamma%*%gamma) * b[g,] + 
              + matrix(
                b[g, ] ^ 2,
                nrow = L,
                ncol = N,
                byrow = FALSE
              ) + 
              t((X_gamma%*%gamma)^2)
          )
        
      }
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      
      # Model Parameters 
      
      gam[(D + 1):(D + G), ] <- t(sw*z) %*% (Y - 0.5) # b
      
      gam[(D+G+1):dim_wh, ] <- t(apply(apply(array(apply(sw*z,2,function(x) x*X_gamma),c(N,P_gamma,G)), # gamma
                                             c(2,3), function(y) crossprod(y, Y-0.5)),
                                       c(1,2), function(xy) sum(xy)))
      
      aa <- aperm(array(sw*z, c(N, G, D)), c(1, 3, 2)) * mu
      bb <- aperm(array(sw*z, c(N, G, D * D)), c(3, 1, 2)) * YY
      
      K[(D + 1):(D + G), (D + 1):(D + G), ] <- # b b
        apply(apply(aperm(array(sw*z, c(
          N, G, L
        )), c(1, 3, 2)) * lambda_xi, c(2, 3), sum), 1, diag)
      
      kk <- 0
      vgamma <- 0
      for (g in 1:G) {
        K[D + g, 1:D, ] <- crossprod(aa[, , g], lambda_xi[, , g]) # b v
        K[1:D, D + g, ] <- K[D + g, 1:D, ] # v b 
        kk <- kk + bb[, , g] %*% lambda_xi[, , g]
        vgamma <- vgamma + apply(array(apply(array(apply(aa[,,g],2,function(x) x*lambda_xi[,,g]),c(N,L,D)),
                                             c(2,3), function(xx) xx*X_gamma),c(N,P_gamma,L,D)),
                                 c(2,3,4), sum)
      }
      
      K[1:D, 1:D, ] <- kk # v v
      
      for (l in 1:L)	{
        
        gam[1:D, l] <- apply(aa * (Y[, l] - 0.5), 2, sum) # v
        
        K[(D+G+1):dim_wh, 1:D, l] = vgamma[,l,] # gamma v
        
        K[1:D, (D+G+1):dim_wh, l] = t(K[(D+G+1):dim_wh, 1:D, l])  # v gamma 
        
        K[(D+1):(G + D), (D+G+1):dim_wh, l] <- t(apply(array(apply(lambda_xi,2,function(x) sw*z*x),c(N,G,L)), # b gamma
                                                       c(2,3), function(xx) crossprod(xx, X_gamma))[,,l])
        
        K[(D+G+1):dim_wh, (D+1):(G + D), l] = t(K[(D+1):(G + D), (D+G+1):dim_wh, l]) # gamma b
        
        valori_gg <- apply(array(apply(array(apply(array(apply(lambda_xi,2,function(x) sw*z*x),c(N,G,L)),
                                                   c(2,3), function(xx) xx*X_gamma),c(N,P_gamma,G,L)),
                                       c(3,4), function(xxx) xxx*X_gamma),c(N,P_gamma,G,L)), # gamma gamma
                           c(2,4), sum)[,l]
        if(P_gamma == 1){
          K[(D+G+1):dim_wh, (D+G+1):dim_wh, l] <- valori_gg
        } else {
          for(jj in seq_len(P_gamma)){
            K[(D+G+jj), (D+G+jj), l] <- valori_gg[jj]
          }
        }
        
        wh[, l] <- -ginv(2 * K[, , l]) %*% gam[, l]
      }
      
      
      v <- wh[1:D, ]
      v <- as.matrix(v)
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
      }
      ord.b=order(apply(b,1,mean))
      b <- b[ord.b,]
      
      gamma=wh[(D+G+1):dim_wh,]
      
      
      # log(P(Y|U,Z,X,xi))
      
      for(g in 1:G){
        
        detC <- apply(C[, , g], 2, function(x)
          det(matrix(x, D, D)))
        
        lxi[,g] <-
          (0.5 * log(detC)) + (0.5 * apply(rbind(C[, , g], t(mu[, , g])), 2, function(x)
            t((x[-(1:(D * D))])) %*% solve(matrix(x[1:(D * D)], nrow = D)) %*% (x[-(
              1:(D * D))]))) + rowSums(
                (log(sigma_xi[, , g])) - (0.5 * xi[, , g]) - (lambda_xi[, , g] *
                                                                xi[, , g]^2) + sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + 
                  (lambda_xi[, , g] * (X_gamma%*%gamma)^2) +
                  (2 * sweep(lambda_xi[, , g] * (X_gamma%*%gamma), MARGIN=2, b[g, ], `*`)) + 
                  (sweep(Y - 0.5, MARGIN = 2, b[g, ], `*`)) +
                  ((Y-0.5) * (X_gamma%*%gamma))
              )
      } 
      
    } # end D > 2
    
    
    # M-step
    
    lk = sum(sw*z*log(eta))
    it = 0; lko = lk
    XXdis = array(0,c(G,(G-1)*ncol(X),N))
    for(i in 1:N){
      XXdis[,,i] = diag(G)[,-1]%*%(diag(G-1)%x%t(X[i,]))
    }
    while((lk-lko>10^-6 & it<100) | it==0){
      it = it+1; lko = lk 
      sc = 0; Fi = 0
      for(i in 1:N){
        pdis = eta[i,]
        sc = sc + sw[i] * t(XXdis[,,i])%*%(z[i,]-pdis)
        Fi = Fi + sw[i] * t(XXdis[,,i])%*%(diag(pdis)-pdis%o%pdis)%*%XXdis[,,i]
      }
      
      dbe = as.vector(ginv(Fi)%*%sc)
      mdbe = max(abs(dbe))
      if(mdbe>0.5) dbe = dbe/mdbe*0.5
      be0 = c(beta)
      flag = TRUE
      while(flag){
        beta = be0+dbe
        Eta = matrix(0,N,G)
        for(i in 1:N){
          
          if(ncol(X)==1) Eta[i,] = XXdis[,,i]*beta
          else Eta[i,] = XXdis[,,i]%*%beta
        }	
        if(max(abs(Eta))>100){
          dbe = dbe/2
          flag = TRUE	
        }else{
          flag = FALSE
        }	        	
      }
      if(iter/10 == floor(iter/10))       print(beta)
      
      beta = matrix(beta, P, G-1)    
      exb <- exp(X %*% beta) # updfe priors
      eta <- cbind(1,exb)/(rowSums(exb)+1)
      
      lk = sum(sw*z*log(eta))
    }
    
    
    # E-step
    
    num <- eta * exp(lxi)
    if(any(is.nan(num))) browser()
    #v[is.nan(v)] <- 0
    den <- apply(num, 1, sum)
    z <- num / den
    
    
    # Log-likelihood
    
    ll <- sum(sw * log(den))
    
    
    # Stopping Criteria
    
    diff <- sum(abs(ll - ll_old))
    
    if(sum((ll - ll_old))<0) print(paste(iter, ll, sum((ll - ll_old)), "aaaaaa"))
    else print(paste(iter, ll, sum((ll - ll_old))))
    
  } # end while
  
  
  # Correction to the log-likelihood: Gauss-Hermite Quadrature
  
  Qd <- npoints ^ D       
  GaHer <- glmmML::ghq(npoints, FALSE)                     
  ugh <- as.matrix(expand.grid(rep(list(GaHer$zeros), D))) 
  ugh.star <- sqrt(2)*ugh                                 
  p.gh = as.matrix(expand.grid(rep(list(GaHer$weights), times = D)))	
  Ww =  (2)^(D/2) * exp(apply(ugh, 1, crossprod)) * apply(p.gh, 1, prod)
  W.prod =  apply(p.gh, 1, prod)
  Phi <- apply(ugh.star, 1, mvtnorm::dmvnorm)               
  
  fxy <- array(0, c(N, Qd, G))
  for (g in 1:G)
  {
    for(i in 1:N){
      if(D>1){
        Agh <- t(tcrossprod(t(v), ugh.star)) + b[g, ] + 
          (X_gamma%*%gamma)[i,]
      } else {
        Agh <- t(tcrossprod(t(v), ugh.star)) + b[g, ] + 
          (X_gamma%*%gamma)[i,]
      }
      
      pgh <- 1 / (1 + exp(-Agh))
      fxy[i, , g] <-
        exp(tcrossprod(Y[i,], log(pgh)) + tcrossprod(1 - Y[i,], log(1 - pgh)))
    }
  }
  
  LLva <- ll
  npar = (G-1)*P + (P_gamma * L) + (G * L) + ((L * D) - (D * (D - 1)/2))
  AICva <- -2 * LLva + 2 * npar
  
  llGH1 <- apply(aperm(apply(fxy,c(1,3),function(x) x*Ww*Phi),c(2,1,3)),c(1,3),sum)
  LL = sum(sw * log(rowSums(t(apply(llGH1,1,function(x) x*eta)))))
  AIC <- -2 * LL + 2 * npar
  
  rownames(b) <- NULL
  colnames(b) <- colnames(b, do.NULL = FALSE, prefix = "Item ")
  rownames(b) <- rownames(b, do.NULL = FALSE, prefix = "Cluster ")
  
  rownames(v) <- rownames(v, do.NULL = FALSE, prefix = "Dim ")
  colnames(v) <- colnames(v, do.NULL = FALSE, prefix = "Item ")
  
  gamma_out <- rbind(rep(0, L), gamma)
  colnames(gamma_out) <- colnames(Y)
  rownames(gamma_out) <- colnames(X)
  
  colnames(beta) <- 2:G
  rownames(beta) <- c(colnames(X))
  
  names(LLva) <- c("Log-Likelihood (variational approximation):")
  names(AICva) <- c("AIC (variational approximation):")
  names(LL) <- c("Log-Likelihood (G-H Quadrature correction):")
  names(AIC) <- c("AIC (G-H Quadrature correction):")
  
  out1 <-
    list(
      b = b,
      v = v,
      gamma = gamma_out,
      beta = beta,
      eta = eta,
      mu = mu,
      C = C,
      z = z,
      ugh.star = ugh.star,
      LL = LL,
      AIC = AIC,
      LLva = LLva,
      AICva = AICva,
      npar = npar
    )
  
  return(out1)
  
}



f_MELTA_vfix_final <- function(Y, X, G, D, tol, maxiter, npoints, sw, 
                               b.init, v.init, gamma.init, beta.init)
{
  
  # Function parameters:
  # Y = matrix of response variables
  # X = covariates
  # G = number of mixture components
  # D = latent trait dimension
  # tol = tolerance level for convergence
  # maxiter = maximum number of iterations
  # npoints = quadrature points
  # sw = sampling weigths
  
  
  N <- nrow(Y)         # units
  L <- ncol(Y)         # items
  P <- ncol(X)         # covariates
  
  
  # Parameter initialization
  
  beta <- beta.init
  
  gamma <- gamma.init
  
  exb <- exp(X %*% beta)
  eta <- cbind(1,exb)/(rowSums(exb)+1)  # priors
  
  z <- matrix(NA, nrow=N, ncol=G)
  for(i in 1:N)   z[i,] <- t(rmultinom(1, size = 1, prob = eta[i,]))
  
  v <- v.init
  
  b=b.init
  
  
  # Variational approximation initialization
  
  xi <- array(20, c(N, L, G))
  sigma_xi <- 1 / (1 + exp(-xi))
  lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
  
  C <- array(0, c(D*D, N, G))
  mu <- array(0, c(N, D, G))
  
  
  # Iterative process
  
  ll <- -Inf
  diff <- 1
  iter <- 0
  cond <- TRUE
  
  print(c(beta))
  
  lxi <- matrix(0, N, G) # variational approx
  
  YY <- array(0, c(D * D, N, G))
  gam <- matrix(0, dim_wh, L)
  K <- array(diag(dim_wh), c(dim_wh, dim_wh, L))
  wh <- matrix(0, dim_wh, L) # matrix to store model parameters
  
  while (diff > tol & iter < maxiter)
  {
    iter <- iter + 1
    ll_old <- ll
    
    ## D=1 ####
    if (D == 1) {
      for (g in 1:G)
      {
        # Posterior Statistics
        
        C[, , g] <-
          1 / (1 - 2 * rowSums(sweep(lambda_xi[, , g], MARGIN = 2, v ^ 2, `*`)))
        mu[, , g] <-
          C[, , g] * rowSums(sweep(
            Y - 0.5 + 2 * lambda_xi[, , g] * (X_gamma%*%gamma) + 
              2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`),
            MARGIN = 2, v, `*`))
        
        YY[, , g] <- matrix(C[, , g] + mu[, , g] ^ 2, ncol = 1)
        
        # Variational Parameters (xi)
        
        xi[, , g] <-  
          (YY[, , g] %*% v^2) + (mu[, , g] %*% (2 * b[g, ] * v)) +
          (apply((t(apply(2 * (X_gamma%*%gamma),1,function(x) x* v))), 2, function(y) y*mu[,,g])) + 
          (2 * t(apply(X_gamma%*%gamma,1,function(xx)  xx * b[g, ]))) + 
          (matrix(b[g, ]^2, nrow = N, ncol = L, byrow = TRUE)) + 
          (matrix((X_gamma%*%gamma)^2, nrow = N, ncol = L, byrow = FALSE))
      }
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      # Model Parameters
      
      gam[1, ] <- colSums(crossprod(sw * z * mu[, 1, ], Y - 0.5)) # v
      gam[2:(1 + G), ] <- crossprod(sw * z, Y - 0.5) # b
      gam[(D+G+1):dim_wh, ] <- t(apply(apply(array(apply(sw * z,2,function(x) x*X_gamma),c(N,P_gamma,G)), # gamma
                                             c(2,3), function(y) crossprod(y, Y-0.5)),
                                       c(1,2), function(xy) sum(xy)))
      
      
      K[1, 1, ] <- apply(aperm(array(sw * z * YY[1, , ], c(N, G, L)), c(1, 3, 2)) * lambda_xi, 2, sum) # vv
      
      K[2:(G + 1), 1, ] <- t(apply(aperm(array( # b v
        sw * z * mu[, 1, ], c(N, G, L)
      ), c(1, 3, 2)) * lambda_xi, c(2, 3), sum))
      
      K[1, 2:(G + 1), ] <- K[2:(G + 1), 1, ] # v b
      
      K[(D+G+1):dim_wh, 1, ] <- apply(apply(array(apply(lambda_xi,2, function(x) x*(sw * z * mu[, 1, ])), # gamma v
                                                  c(N,G,L)), c(2,3), function(xx) crossprod(xx, X_gamma)),
                                      c(1,3),sum)
      
      K[1, (D+G+1):dim_wh, ] <- K[(D+G+1):dim_wh, 1, ] # v gamma
      
      for(l in 1:L){
        
        diag(K[2:(G + 1),2:(G + 1),l]) <- apply(array(apply(lambda_xi,2,function(x) sw*z*x),c(N,G,L)), # b b
                                                c(2,3), sum)[,l]
        
        K[2:(G + 1), (D+G+1):dim_wh, l] <- t(apply(array(apply(lambda_xi,2,function(x) sw*z*x),c(N,G,L)), # b gamma
                                                   c(2,3), function(xx) crossprod(xx, X_gamma))[,,l])
        
        K[(D+G+1):dim_wh, 2:(G + 1), l] = t(K[2:(G + 1), (D+G+1):dim_wh, l]) # gamma b
        
        valori_gg <- apply(array(apply(array(apply(array(apply(lambda_xi,2,function(x) sw*z*x),c(N,G,L)),
                                                   c(2,3), function(xx) xx*X_gamma),c(N,P_gamma,G,L)),
                                       c(3,4), function(xxx) xxx*X_gamma),c(N,P_gamma,G,L)), # gamma gamma
                           c(2,4), sum)[,l]
        if(P_gamma == 1){
          K[(D+G+1):dim_wh, (D+G+1):dim_wh, l] <- valori_gg
        } else {
          for(jj in seq_len(P_gamma)){
            K[(D+G+jj), (D+G+jj), l] <- valori_gg[jj]
          }
        }
        
        wh[, l] <- -ginv(2 * K[, , l]) %*% gam[, l] # parameter updates
      }
      
      
      v <- wh[1:D, ]
      v <- as.matrix(t(v))
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
      }
      ord.b=order(apply(b,1,mean))
      b <- b[ord.b,]
      
      gamma=wh[(D+G+1):dim_wh,]
      
      
      # log(P(Y|Z,X,xi))
      
      for(g in 1:G){
        lxi[,g] <-
          (0.5 * log(C[1, , g])) + (mu[, 1, g] ^ 2 / (2 * C[1, , g])) + rowSums(
            (log(sigma_xi[, , g])) - (0.5 * xi[, , g]) - (lambda_xi[, , g] * xi[, , g] ^ 2) + 
              (sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`)) +
              (lambda_xi[, , g] * (X_gamma%*%gamma)^2) +
              (2 * sweep(lambda_xi[, , g] * (X_gamma%*%gamma), MARGIN=2, b[g, ], `*`)) + 
              (sweep(Y - 0.5, MARGIN = 2, b[g, ], `*`)) + 
              ((Y-0.5) * (X_gamma%*%gamma)))
      } 
    } # end D=1
    
    ## D=2 ####
    if (D == 2) {
      for (g in 1:G)
      {
        # Latent Posterior Statistics
        
        C[, , g] <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(v), t(v))))
        
        mu[, , g] <-
          t(apply(rbind(C[, , g], 
                        tcrossprod(v, Y - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`)+
                                     2 * (lambda_xi[, , g] * (X_gamma%*%gamma)))), 
                  2, function(x) matrix(x[1:4], nrow = D) %*% x[-(1:4)]))
        
        # Variational Parameters (xi)
        
        YY[, , g] <- C[, , g] + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY[, , g], 2, function(x)
              rowSums(crossprod(v, matrix(x, ncol = D)) * t(v))) + 
              tcrossprod(2 * b[g, ] * t(v), mu[, , g]) + 
              2 * t(X_gamma%*%gamma) * tcrossprod(t(v),mu[,,g]) + 
              2 * t(X_gamma%*%gamma) * b[g,] + 
              + matrix(
                b[g, ] ^ 2,
                nrow = L,
                ncol = N,
                byrow = FALSE
              ) + 
              t((X_gamma%*%gamma)^2)
          )
      }
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      
      # Model Parameters 
      
      gam[3:(2 + G), ] <- crossprod(sw*z, Y - 0.5) # b
      
      gam[(D+G+1):dim_wh, ] <- t(apply(apply(array(apply(sw*z,2,function(x) x*X_gamma),c(N,P_gamma,G)), # gamma
                                             c(2,3), function(y) crossprod(y, Y-0.5)),
                                       c(1,2), function(xy) sum(xy)))
      
      aa <- aperm(array(sw*z, c(N, G, D)), c(1, 3, 2)) * mu
      bb <- aperm(array(sw*z, c(N, G, D * D)), c(3, 1, 2)) * YY
      
      K[3:(2 + G), 3:(2 + G), ] <- # b b
        apply(apply(aperm(array(sw*z, c(
          N, G, L
        )), c(1, 3, 2)) * lambda_xi, c(2, 3), sum), 1, diag)
      
      kk <- 0
      vgamma <- 0
      for (g in 1:G) {
        
        K[2 + g, 1:2, ] <- crossprod(aa[, , g], lambda_xi[, , g]) # b v
        K[1:2, 2 + g, ] <- K[2 + g, (1:2), ] # v b
        
        kk <- kk + bb[, , g] %*% lambda_xi[, , g]
        
        vgamma <- vgamma + apply(array(apply(array(apply(aa[,,g],2,function(x) x*lambda_xi[,,g]),c(N,L,D)),
                                             c(2,3), function(xx) xx*X_gamma),c(N,P_gamma,L,D)),
                                 c(2,3,4), sum)
      }
      
      K[1:2, 1:2, ] <- kk # v v
      
      for (l in 1:L)
      {
        gam[1:D, l] <- apply(aa * (Y[, l] - 0.5), 2, sum) # v
        
        K[(D+G+1):dim_wh, 1:D, l] = vgamma[,l,] # gamma v
        K[1:D, (D+G+1):dim_wh, l] = t(K[(D+G+1):dim_wh, 1:D, l])  # v gamma 
        
        K[(D+1):(G + D), (D+G+1):dim_wh, l] <- t(apply(array(apply(lambda_xi,2,function(x) sw*z*x),c(N,G,L)), # b gamma
                                                       c(2,3), function(xx) crossprod(xx, X_gamma))[,,l])
        K[(D+G+1):dim_wh, (D+1):(G + D), l] = t(K[(D+1):(G + D), (D+G+1):dim_wh, l]) # gamma b
        
        valori_gg <- apply(array(apply(array(apply(array(apply(lambda_xi,2,function(x) sw*z*x),c(N,G,L)),
                                                   c(2,3), function(xx) xx*X_gamma),c(N,P_gamma,G,L)),
                                       c(3,4), function(xxx) xxx*X_gamma),c(N,P_gamma,G,L)), # gamma gamma
                           c(2,4), sum)[,l]
        if(P_gamma == 1){
          K[(D+G+1):dim_wh, (D+G+1):dim_wh, l] <- valori_gg
        } else {
          for(jj in seq_len(P_gamma)){
            K[(D+G+jj), (D+G+jj), l] <- valori_gg[jj]
          }
        }
        
        wh[, l] <- -ginv(2 * K[, , l]) %*% gam[, l]
      }
      
      v <- wh[1:D, ]
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
      }
      ord.b=order(apply(b,1,mean))
      b <- b[ord.b,]
      
      gamma=wh[(D+G+1):dim_wh,]
      
      
      # log(P(Y|Z,X,xi))
      
      for(g in 1:G){
        
        detC <- C[1, , g] * C[4, , g] - C[3, , g] * C[2, , g]
        
        lxi[,g] <-
          (0.5 * log(detC)) + (0.5 * apply(rbind(C[4, , g] / detC, -C[2, , g] / detC, -C[3, , g] /
                                                   detC, C[1, , g] / detC, t(mu[, , g])), 2, function(x)
                                                     t((x[-(1:4)])) %*% matrix(x[1:4], nrow = D) %*% (x[-(1:4)]))) + 
          rowSums(
            (log(sigma_xi[, , g])) - (0.5 * xi[, , g]) - (lambda_xi[, , g] * xi[, , g]^2) + 
              (sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`)) + 
              (lambda_xi[, , g] * (X_gamma%*%gamma)^2) +
              (2 * sweep(lambda_xi[, , g] * (X_gamma%*%gamma), MARGIN=2, b[g, ], `*`)) + 
              (sweep(Y - 0.5, MARGIN = 2, b[g, ], `*`)) +
              ((Y-0.5) * (X_gamma%*%gamma))
          )
      } 
    }# end D=2
    
    ## D>2 ####
    if (D > 2) {
      for (g in 1:G)
      {
        # Latent Posterior Statistics
        
        C[, , g] <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(v), t(v))))
        
        mu[, , g] <-
          t(apply(rbind(C[, , g], v %*% t(
            Y - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`) +
              2 * lambda_xi[, , g] * (X_gamma%*%gamma)
          )), 2, function(x)
            matrix(x[1:(D * D)], nrow = D) %*% x[-(1:(D * D))]))
        
        # Variational Parameters (xi)
        
        YY[, , g] <- C[, , g] + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY[, , g], 2, function(x)
              rowSums((t(v) %*% matrix(x, ncol = D)) * t(v))) + 
              (2 * b[g, ] * t(v)) %*% t(mu[, , g]) + 
              2 * t(X_gamma%*%gamma) * tcrossprod(t(v),mu[,,g]) + 
              2 * t(X_gamma%*%gamma) * b[g,] + 
              + matrix(
                b[g, ] ^ 2,
                nrow = L,
                ncol = N,
                byrow = FALSE
              ) + 
              t((X_gamma%*%gamma)^2)
          )
        
      }
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      
      # Model Parameters 
      
      gam[(D + 1):(D + G), ] <- t(sw*z) %*% (Y - 0.5) # b
      
      gam[(D+G+1):dim_wh, ] <- t(apply(apply(array(apply(sw*z,2,function(x) x*X_gamma),c(N,P_gamma,G)), # gamma
                                             c(2,3), function(y) crossprod(y, Y-0.5)),
                                       c(1,2), function(xy) sum(xy)))
      
      aa <- aperm(array(sw*z, c(N, G, D)), c(1, 3, 2)) * mu
      bb <- aperm(array(sw*z, c(N, G, D * D)), c(3, 1, 2)) * YY
      
      K[(D + 1):(D + G), (D + 1):(D + G), ] <- # b b
        apply(apply(aperm(array(sw*z, c(
          N, G, L
        )), c(1, 3, 2)) * lambda_xi, c(2, 3), sum), 1, diag)
      
      kk <- 0
      vgamma <- 0
      for (g in 1:G) {
        K[D + g, 1:D, ] <- crossprod(aa[, , g], lambda_xi[, , g]) # b v
        K[1:D, D + g, ] <- K[D + g, 1:D, ] # v b 
        kk <- kk + bb[, , g] %*% lambda_xi[, , g]
        vgamma <- vgamma + apply(array(apply(array(apply(aa[,,g],2,function(x) x*lambda_xi[,,g]),c(N,L,D)),
                                             c(2,3), function(xx) xx*X_gamma),c(N,P_gamma,L,D)),
                                 c(2,3,4), sum)
      }
      
      K[1:D, 1:D, ] <- kk # v v
      
      for (l in 1:L)	{
        
        gam[1:D, l] <- apply(aa * (Y[, l] - 0.5), 2, sum) # v
        
        K[(D+G+1):dim_wh, 1:D, l] = vgamma[,l,] # gamma v
        
        K[1:D, (D+G+1):dim_wh, l] = t(K[(D+G+1):dim_wh, 1:D, l])  # v gamma 
        
        K[(D+1):(G + D), (D+G+1):dim_wh, l] <- t(apply(array(apply(lambda_xi,2,function(x) sw*z*x),c(N,G,L)), # b gamma
                                                       c(2,3), function(xx) crossprod(xx, X_gamma))[,,l])
        
        K[(D+G+1):dim_wh, (D+1):(G + D), l] = t(K[(D+1):(G + D), (D+G+1):dim_wh, l]) # gamma b
        
        valori_gg <- apply(array(apply(array(apply(array(apply(lambda_xi,2,function(x) sw*z*x),c(N,G,L)),
                                                   c(2,3), function(xx) xx*X_gamma),c(N,P_gamma,G,L)),
                                       c(3,4), function(xxx) xxx*X_gamma),c(N,P_gamma,G,L)), # gamma gamma
                           c(2,4), sum)[,l]
        if(P_gamma == 1){
          K[(D+G+1):dim_wh, (D+G+1):dim_wh, l] <- valori_gg
        } else {
          for(jj in seq_len(P_gamma)){
            K[(D+G+jj), (D+G+jj), l] <- valori_gg[jj]
          }
        }
        
        wh[, l] <- -ginv(2 * K[, , l]) %*% gam[, l]
      }
      
      
      v <- wh[1:D, ]
      v <- as.matrix(v)
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
      }
      ord.b=order(apply(b,1,mean))
      b <- b[ord.b,]
      
      gamma=wh[(D+G+1):dim_wh,]
      
      
      # log(P(Y|U,Z,X,xi))
      
      for(g in 1:G){
        
        detC <- apply(C[, , g], 2, function(x)
          det(matrix(x, D, D)))
        
        lxi[,g] <-
          (0.5 * log(detC)) + (0.5 * apply(rbind(C[, , g], t(mu[, , g])), 2, function(x)
            t((x[-(1:(D * D))])) %*% solve(matrix(x[1:(D * D)], nrow = D)) %*% (x[-(
              1:(D * D))]))) + rowSums(
                (log(sigma_xi[, , g])) - (0.5 * xi[, , g]) - (lambda_xi[, , g] *
                                                                xi[, , g]^2) + sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + 
                  (lambda_xi[, , g] * (X_gamma%*%gamma)^2) +
                  (2 * sweep(lambda_xi[, , g] * (X_gamma%*%gamma), MARGIN=2, b[g, ], `*`)) + 
                  (sweep(Y - 0.5, MARGIN = 2, b[g, ], `*`)) +
                  ((Y-0.5) * (X_gamma%*%gamma))
              )
      } 
      
    } # end D > 2
    
    
    # M-step
    
    lk = sum(sw*z*log(eta))
    it = 0; lko = lk
    XXdis = array(0,c(G,(G-1)*ncol(X),N))
    for(i in 1:N){
      XXdis[,,i] = diag(G)[,-1]%*%(diag(G-1)%x%t(X[i,]))
    }
    while((lk-lko>10^-6 & it<100) | it==0){
      it = it+1; lko = lk 
      sc = 0; Fi = 0
      for(i in 1:N){
        pdis = eta[i,]
        sc = sc + sw[i] * t(XXdis[,,i])%*%(z[i,]-pdis)
        Fi = Fi + sw[i] * t(XXdis[,,i])%*%(diag(pdis)-pdis%o%pdis)%*%XXdis[,,i]
      }
      
      dbe = as.vector(ginv(Fi)%*%sc)
      mdbe = max(abs(dbe))
      if(mdbe>0.5) dbe = dbe/mdbe*0.5
      be0 = c(beta)
      flag = TRUE
      while(flag){
        beta = be0+dbe
        Eta = matrix(0,N,G)
        for(i in 1:N){
          
          if(ncol(X)==1) Eta[i,] = XXdis[,,i]*beta
          else Eta[i,] = XXdis[,,i]%*%beta
        }	
        if(max(abs(Eta))>100){
          dbe = dbe/2
          flag = TRUE	
        }else{
          flag = FALSE
        }	        	
      }
      if(iter/10 == floor(iter/10))       print(beta)
      
      beta = matrix(beta, P, G-1)    
      exb <- exp(X %*% beta) # updfe priors
      eta <- cbind(1,exb)/(rowSums(exb)+1)
      
      lk = sum(sw*z*log(eta))
    }
    
    
    # E-step
    
    num <- eta * exp(lxi)
    if(any(is.nan(num))) browser()
    #v[is.nan(v)] <- 0
    den <- apply(num, 1, sum)
    z <- num / den
  
    
    # Log-likelihood
    
    ll <- sum(sw * log(den))
    
    
    # Stopping Criteria
    
    diff <- sum(abs(ll - ll_old))
    
    if(sum((ll - ll_old))<0) print(paste(iter, ll, sum((ll - ll_old)), "aaaaaa"))
    else print(paste(iter, ll, sum((ll - ll_old))))
    
  } # end while
  
  
  # Correction to the log-likelihood: Gauss-Hermite Quadrature
  
  Qd <- npoints ^ D       
  GaHer <- glmmML::ghq(npoints, FALSE)                     
  ugh <- as.matrix(expand.grid(rep(list(GaHer$zeros), D))) 
  ugh.star <- sqrt(2)*ugh                                 
  p.gh = as.matrix(expand.grid(rep(list(GaHer$weights), times = D)))	
  Ww =  (2)^(D/2) * exp(apply(ugh, 1, crossprod)) * apply(p.gh, 1, prod)
  W.prod =  apply(p.gh, 1, prod)
  Phi <- apply(ugh.star, 1, mvtnorm::dmvnorm)               
  
  fxy <- array(0, c(N, Qd, G))
  for (g in 1:G)
  {
    for(i in 1:N){
      if(D>1){
        Agh <- t(tcrossprod(t(v), ugh.star)) + b[g, ] + 
          (X_gamma%*%gamma)[i,]
      } else {
        Agh <- t(tcrossprod(t(v), ugh.star)) + b[g, ] + 
          (X_gamma%*%gamma)[i,]
      }
      
      pgh <- 1 / (1 + exp(-Agh))
      fxy[i, , g] <-
        exp(tcrossprod(Y[i,], log(pgh)) + tcrossprod(1 - Y[i,], log(1 - pgh)))
    }
  }
  
  LLva <- ll
  npar = (G-1)*P + (P_gamma * L) + (G * L) + ((L * D) - (D * (D - 1)/2))
  AICva <- -2 * LLva + 2 * npar
  
  llGH1 <- apply(aperm(apply(fxy,c(1,3),function(x) x*Ww*Phi),c(2,1,3)),c(1,3),sum)
  LL = sum(sw * log(rowSums(t(apply(llGH1,1,function(x) x*eta)))))
  AIC <- -2 * LL + 2 * npar
  
  rownames(b) <- NULL
  colnames(b) <- colnames(b, do.NULL = FALSE, prefix = "Item ")
  rownames(b) <- rownames(b, do.NULL = FALSE, prefix = "Cluster ")
  
  rownames(v) <- rownames(v, do.NULL = FALSE, prefix = "Dim ")
  colnames(v) <- colnames(v, do.NULL = FALSE, prefix = "Item ")
  
  gamma_out <- rbind(rep(0, L), gamma)
  colnames(gamma_out) <- colnames(Y)
  rownames(gamma_out) <- colnames(X)
  
  colnames(beta) <- 2:G
  rownames(beta) <- c(colnames(X))
  
  names(LLva) <- c("Log-Likelihood (variational approximation):")
  names(AICva) <- c("AIC (variational approximation):")
  names(LL) <- c("Log-Likelihood (G-H Quadrature correction):")
  names(AIC) <- c("AIC (G-H Quadrature correction):")
  
  out1 <-
    list(
      b = b,
      v = v,
      gamma = gamma_out,
      beta = beta,
      eta = eta,
      mu = mu,
      C = C,
      z = z,
      ugh.star = ugh.star,
      LL = LL,
      AIC = AIC,
      LLva = LLva,
      AICva = AICva, 
      npar = npar
    )
  
  return(out1)
  
}


f_MELTA_g1 <- function(Y, X, D, tol, maxiter, npoints, sw)
{
  
  # Function parameters:
  # Y = matrix of response variables
  # X = covariates
  # D = latent trait dimension
  # tol = tolerance level for convergence
  # maxiter = maximum number of iterations
  # npoints = quadrature points
  # sw = sampling weights
  
  
  N <- nrow(Y)         # units
  L <- ncol(Y)         # items
  P <- ncol(X)         # covariates (including intercept)
  
  # Remove intercept from X for gamma (identifiability: gamma_0k = 0)
  X_gamma <- as.matrix(X[, -1, drop=FALSE])
  P_gamma <- ncol(X_gamma)  # P - 1
  dim_wh  <- D + 1 + P_gamma
  
  # Parameter initialization
  
  gamma <- matrix(rep(0, P_gamma*L), nrow=P_gamma, ncol=L)
  
  v <- matrix(rnorm(L * D), D, L)
  
  b=matrix(rnorm(L,0, 0.000001), 1, L) # matrix(rnorm(G*L), G, L)
  
  
  # Variational approximation initialization
  
  xi <- matrix(20, N, L)
  sigma_xi <- 1 / (1 + exp(-xi))
  lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
  
  C <- matrix(0, D*D, N)
  mu <- matrix(0, N, D)
  
  
  # Iterative process
  
  ll <- -Inf
  diff <- 1
  iter <- 0
  cond <- TRUE
  
  lxi <- matrix(0, N, 1) # variational approx
  
  YY <- matrix(0, D * D, N)
  gam <- matrix(0, dim_wh, L)
  K <- array(diag(dim_wh), c(dim_wh, dim_wh, L))
  wh <- matrix(0, dim_wh, L) # matrix to store model parameters
  
  while (diff > tol & iter < maxiter)
  {
    iter <- iter + 1
    ll_old <- ll
    
    ## D=1 ####
    if (D == 1) {
      # Posterior Statistics
      
      C <-
        1 / (1 - 2 * rowSums(sweep(lambda_xi, MARGIN = 2, v ^ 2, `*`)))
      mu <-
        C * rowSums(sweep(
          Y - 0.5 + 2 * lambda_xi * (X_gamma%*%gamma) + 
            2 * sweep(lambda_xi, MARGIN = 2, b, `*`),
          MARGIN = 2, v, `*`))
      
      YY <- matrix(C + mu ^ 2, ncol = 1)
      
      # Variational Parameters (xi)
      
      xi <-  
        (YY %*% v^2) + (mu %*% (2 * b * v)) +
        (apply((t(apply(2 * (X_gamma%*%gamma),1,function(x) x* v))), 2, function(y) y*mu)) + 
        (2 * t(apply(X_gamma%*%gamma,1,function(xx)  xx * b))) + 
        (matrix(b^2, nrow = N, ncol = L, byrow = TRUE)) + 
        (matrix((X_gamma%*%gamma)^2, nrow = N, ncol = L, byrow = FALSE))
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      # Model Parameters
      
      gam[1, ] <- colSums(crossprod(sw * mu, Y - 0.5)) # v
      gam[2, ] <- crossprod(sw, Y - 0.5) # b
      gam[3:(P_gamma+2), ] <- t(apply(apply(X_gamma,2,function(x) x*sw), # gamma
                                      +           c(2), function(y) crossprod(y, Y-0.5)))
      
      
      K[1, 1, ] <- apply(sw * YY[,1] * lambda_xi, 2, sum) # vv
      
      K[2, 1, ] <- apply( # b v
        sw * mu * lambda_xi, c(2), sum)
      
      K[1, 2, ] <- K[2, 1, ] # v b
      
      K[3:(P_gamma+2), 1, ] <- apply(lambda_xi * sw * mu, # gamma v
                                     c(2), function(xx) crossprod(xx, X_gamma))
      
      K[1, 3:(P_gamma+2), ] <- K[3:(P_gamma+2), 1, ] # v gamma
      
      for(l in 1:L){
        
        K[2,2,l] <- apply(lambda_xi*sw,2,sum)[l]
        
        K[2, (3):(P_gamma+2), l] <- apply(lambda_xi * sw, # b gamma
                                          c(2), function(xx) crossprod(xx, X_gamma))[,l]
        
        K[(3):(P_gamma+2), 2, l] = t(K[2, (3):(P_gamma+2), l]) # gamma b
        
        valori_gg <- apply(array(apply(array(apply(lambda_xi * sw,
                                                   c(2), function(xx) xx*X_gamma),c(N,P_gamma,L)),
                                       c(3), function(xxx) xxx*X_gamma),c(N,P_gamma,L)), # gamma gamma
                           c(2,3), sum)[,l]
        if(P_gamma == 1){
          K[(D+2):dim_wh, (D+2):dim_wh, l] <- valori_gg
        } else {
          for(jj in seq_len(P_gamma)){
            K[(D+1+jj), (D+1+jj), l] <- valori_gg[jj]
          }
        }
        
        wh[, l] <- -ginv(2 * K[, , l]) %*% gam[, l] # parameter updates
      }
      
      
      v <- wh[1, ]
      v <- as.matrix(t(v))
      
      b <- wh[2, ]
      b <- as.matrix(t(b))
      
      gamma=wh[(D+2):dim_wh,]
      
      
      # log(P(Y|Z,X,xi))
      
      lxi <-
        (0.5 * log(C)) + (mu ^ 2 / (2 * C)) + rowSums(
          (log(sigma_xi)) - (0.5 * xi) - (lambda_xi * xi ^ 2) + 
            (sweep(lambda_xi, MARGIN = 2, b ^ 2, `*`)) +
            (lambda_xi * (X_gamma%*%gamma)^2) +
            (2 * sweep(lambda_xi * (X_gamma%*%gamma), MARGIN=2, b, `*`)) + 
            (sweep(Y - 0.5, MARGIN = 2, b, `*`)) + 
            ((Y-0.5) * (X_gamma%*%gamma)))
    } # end D=1
    
    ## D=2 ####
    if (D == 2) {
      
      # Latent Posterior Statistics
      
      C <-
        apply(lambda_xi, 1, function(x)
          solve(diag(D) - 2 * crossprod(x * t(v), t(v))))
      
      mu <-
        t(apply(rbind(C, 
                      tcrossprod(v, Y - 0.5 + 2 * sweep(lambda_xi, MARGIN = 2, b, `*`)+
                                   2 * (lambda_xi * (X_gamma%*%gamma)))), 
                2, function(x) matrix(x[1:4], nrow = D) %*% x[-(1:4)]))
      
      # Variational Parameters (xi)
      
      YY <- C + apply(mu, 1, tcrossprod)
      
      b=c(b)
      xi <-
        t(
          apply(YY, 2, function(x)
            rowSums(crossprod(v, matrix(x, ncol = D)) * t(v))) + 
            tcrossprod(2 * b * t(v), mu) + 
            2 * t(X_gamma%*%gamma) * tcrossprod(t(v),mu) + 
            2 * t(X_gamma%*%gamma) * b + 
            + matrix(
              b ^ 2,
              nrow = L,
              ncol = N,
              byrow = FALSE
            ) + 
            t((X_gamma%*%gamma)^2)
        )
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      
      # Model Parameters 
      
      gam[3, ] <- crossprod(sw, Y - 0.5) # b
      
      gam[4:(P_gamma+3), ] <- t(apply(sw*X_gamma, # gamma
                                      c(2), function(y) crossprod(y, Y-0.5)))
      
      aa <- matrix(sw, N, D) * mu
      bb <- t(matrix(sw, N, D * D)) * YY
      
      K[3, 3, ] <- # b b
        apply(sw * lambda_xi, c(2), sum)
      
      kk <- 0
      vgamma <- 0
      
      K[3, 1:2, ] <- crossprod(aa, lambda_xi) # b v
      K[1:2, 3, ] <- K[3, (1:2), ] # v b
      
      kk <- kk + bb %*% lambda_xi
      
      vgamma <- vgamma + apply(array(apply(array(apply(aa,2,function(x) x*lambda_xi),c(N,L,D)),
                                           c(2,3), function(xx) xx*X_gamma),c(N,P_gamma,L,D)),
                               c(2,3,4), sum)
      
      K[1:2, 1:2, ] <- kk # v v
      
      for (l in 1:L)
      {
        gam[1:D, l] <- apply(aa * (Y[, l] - 0.5), 2, sum) # v
        
        K[4:(P_gamma+3), 1:D, l] = vgamma[,l,] # gamma v
        K[1:D, 4:(P_gamma+3), l] = t(K[4:(P_gamma+3), 1:D, l])  # v gamma 
        
        K[3, 4:(P_gamma+3), l] <- apply(lambda_xi * sw, # b gamma
                                        c(2), function(xx) crossprod(xx, X_gamma))[,l]
        K[4:(P_gamma+3), 3, l] = t(K[3, 4:(P_gamma+3), l]) # gamma b
        
        valori_gg <- apply(array(apply(array(apply(lambda_xi * sw,
                                                   c(2), function(xx) xx*X_gamma),c(N,P_gamma,L)),
                                       c(3), function(xxx) xxx*X_gamma),c(N,P_gamma,L)), # gamma gamma
                           c(2,3), sum)[,l]
        if(P_gamma == 1){
          K[(D+2):dim_wh, (D+2):dim_wh, l] <- valori_gg
        } else {
          for(jj in seq_len(P_gamma)){
            K[(D+1+jj), (D+1+jj), l] <- valori_gg[jj]
          }
        }
        
        wh[, l] <- -ginv(2 * K[, , l]) %*% gam[, l]
      }
      
      v <- wh[1:D, ]
      
      b <- wh[3, ]
      
      gamma=wh[4:(P_gamma+3),]
      
      
      # log(P(Y|Z,X,xi))
      
      detC <- C[1, ] * C[4, ] - C[3, ] * C[2, ]
      
      lxi <-
        (0.5 * log(detC)) + (0.5 * apply(rbind(C[4, ] / detC, -C[2, ] / detC, -C[3, ] /
                                                 detC, C[1, ] / detC, t(mu)), 2, function(x)
                                                   t((x[-(1:4)])) %*% matrix(x[1:4], nrow = D) %*% (x[-(1:4)]))) + 
        rowSums(
          (log(sigma_xi)) - (0.5 * xi) - (lambda_xi * xi^2) + 
            (sweep(lambda_xi, MARGIN = 2, b ^ 2, `*`)) + 
            (lambda_xi * (X_gamma%*%gamma)^2) +
            (2 * sweep(lambda_xi * (X_gamma%*%gamma), MARGIN=2, b, `*`)) + 
            (sweep(Y - 0.5, MARGIN = 2, b, `*`)) +
            ((Y-0.5) * (X_gamma%*%gamma))
        )
    }# end D=2
    
    ## D>2 ####
    if (D > 2) {
      
      # Latent Posterior Statistics
      
      C <-
        apply(lambda_xi, 1, function(x)
          solve(diag(D) - 2 * crossprod(x * t(v), t(v))))
      
      mu <-
        t(apply(rbind(C, v %*% t(
          Y - 0.5 + 2 * sweep(lambda_xi, MARGIN = 2, b, `*`) +
            2 * lambda_xi * (X_gamma%*%gamma)
        )), 2, function(x)
          matrix(x[1:(D * D)], nrow = D) %*% x[-(1:(D * D))]))
      
      # Variational Parameters (xi)
      
      YY <- C + apply(mu, 1, tcrossprod)
      
      b=c(b)
      xi <-
        t(
          apply(YY, 2, function(x)
            rowSums((t(v) %*% matrix(x, ncol = D)) * t(v))) + 
            (2 * b * t(v)) %*% t(mu) + 
            2 * t(X_gamma%*%gamma) * tcrossprod(t(v),mu) + 
            2 * t(X_gamma%*%gamma) * b + 
            + matrix(
              b ^ 2,
              nrow = L,
              ncol = N,
              byrow = FALSE
            ) + 
            t((X_gamma%*%gamma)^2)
        )
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      
      # Model Parameters 
      
      gam[(D + 1), ] <- t(sw) %*% (Y - 0.5) # b
      
      gam[(D+2):dim_wh, ] <- t(apply(sw*X_gamma, # gamma
                                     c(2), function(y) crossprod(y, Y-0.5)))
      
      aa <- matrix(sw, N, D) * mu
      bb <- t(matrix(sw, N, D * D)) * YY
      
      K[(D + 1), (D + 1), ] <- # b b
        apply(sw * lambda_xi, c(2), sum)
      
      kk <- 0
      vgamma <- 0
      
      K[D + 1, 1:D, ] <- crossprod(aa, lambda_xi) # b v
      K[1:D, D + 1, ] <- K[D + 1, 1:D, ] # v b 
      kk <- kk + bb %*% lambda_xi
      vgamma <- vgamma + apply(array(apply(array(apply(aa,2,function(x) x*lambda_xi),c(N,L,D)),
                                           c(2,3), function(xx) xx*X_gamma),c(N,P_gamma,L,D)),
                               c(2,3,4), sum)
      
      K[1:D, 1:D, ] <- kk # v v
      
      for (l in 1:L)	{
        
        gam[1:D, l] <- apply(aa * (Y[, l] - 0.5), 2, sum) # v
        
        K[(D+2):dim_wh, 1:D, l] = vgamma[,l,] # gamma v
        
        K[1:D, (D+2):dim_wh, l] = t(K[(D+2):dim_wh, 1:D, l])  # v gamma 
        
        K[(D+1), (D+2):dim_wh, l] <- apply(lambda_xi * sw, # b gamma
                                           c(2), function(xx) crossprod(xx, X_gamma))[,l]
        
        K[(D+2):dim_wh, (D+1), l] = t(K[(D+1), (D+2):dim_wh, l]) # gamma b
        
        # gamma gamma diagonal G=1 (loop for P_gamma=1 compatibility)
        valori_gg <- apply(array(apply(array(apply(lambda_xi * sw,
                                                   c(2), function(xx) xx*X_gamma),c(N,P_gamma,L)),
                                       c(3), function(xxx) xxx*X_gamma),c(N,P_gamma,L)), # gamma gamma
                           c(2,3), sum)[,l]
        if(P_gamma == 1){
          K[(D+2):dim_wh, (D+2):dim_wh, l] <- valori_gg
        } else {
          for(jj in seq_len(P_gamma)){
            K[(D+1+jj), (D+1+jj), l] <- valori_gg[jj]
          }
        }
        
        wh[, l] <- -ginv(2 * K[, , l]) %*% gam[, l]
      }
      
      
      v <- wh[1:D, ]
      v <- as.matrix(v)
      
      b <- wh[D + 1, ]
      
      gamma=wh[(D+2):dim_wh,]
      
      
      # log(P(Y|U,Z,X,xi))
      
      detC <- apply(C, 2, function(x)
        det(matrix(x, D, D)))
      
      lxi <-
        (0.5 * log(detC)) + (0.5 * apply(rbind(C, t(mu)), 2, function(x)
          t((x[-(1:(D * D))])) %*% solve(matrix(x[1:(D * D)], nrow = D)) %*% (x[-(
            1:(D * D))]))) + rowSums(
              (log(sigma_xi)) - (0.5 * xi) - (lambda_xi *
                                                xi^2) + sweep(lambda_xi, MARGIN = 2, b ^ 2, `*`) + 
                (lambda_xi * (X_gamma%*%gamma)^2) +
                (2 * sweep(lambda_xi * (X_gamma%*%gamma), MARGIN=2, b, `*`)) + 
                (sweep(Y - 0.5, MARGIN = 2, b, `*`)) +
                ((Y-0.5) * (X_gamma%*%gamma))
            )
      
    } # end D > 2
    
    
    # Log-likelihood
    
    ll <- sum(sw * lxi)
    
    
    # Stopping Criteria
    
    diff <- sum(abs(ll - ll_old))
    
    if(sum((ll - ll_old))<0) print(paste(iter, ll, sum((ll - ll_old)), "aaaaaa"))
    else print(paste(iter, ll, sum((ll - ll_old))))
    
  } # end while
  
  
  # Correction to the log-likelihood: Gauss-Hermite Quadrature
  
  Qd <- npoints ^ D       
  GaHer <- glmmML::ghq(npoints, FALSE)                     
  ugh <- as.matrix(expand.grid(rep(list(GaHer$zeros), D))) 
  ugh.star <- sqrt(2)*ugh                                 
  p.gh = as.matrix(expand.grid(rep(list(GaHer$weights), times = D)))	
  Ww =  (2)^(D/2) * exp(apply(ugh, 1, crossprod)) * apply(p.gh, 1, prod)
  W.prod =  apply(p.gh, 1, prod)
  Phi <- apply(ugh.star, 1, mvtnorm::dmvnorm)               
  
  fxy <- matrix(0, N, Qd)
  for(i in 1:N){
    if(D>1){
      Agh <- t(tcrossprod(t(v), ugh.star)) + b + 
        (X_gamma%*%gamma)[i,]
    } else {
      Agh <- t(apply(tcrossprod(t(v), ugh.star), c(2), function(x) x + b + 
                       (X_gamma%*%gamma)[i,]))
    }
    pgh <- 1 / (1 + exp(-Agh))
    fxy[i, ] <-
      exp(tcrossprod(Y[i,], log(pgh)) + tcrossprod(1 - Y[i,], log(1 - pgh)))
  }
  
  LLva <- ll
  npar = (P_gamma * L) + (L) + ((L * D) - (D * (D - 1)/2))
  AICva <- -2 * LLva + 2 * npar
  
  fxy[!is.finite(fxy)] <- 0
  llGH1 <- rowSums(t(apply(fxy,c(1),function(x) x*Ww*Phi)))
  LL = sum(sw * log(llGH1))
  AIC <- -2 * LL + 2 * npar
  
  b=matrix(b,ncol=L)
  colnames(v) <- colnames(v, do.NULL = FALSE, prefix = "Item ")
  
  rownames(v) <- rownames(v, do.NULL = FALSE, prefix = "Dim ")
  colnames(v) <- colnames(v, do.NULL = FALSE, prefix = "Item ")
  
  gamma_out <- rbind(rep(0, L), gamma)
  colnames(gamma_out) <- colnames(Y)
  rownames(gamma_out) <- colnames(X)
  
  names(LLva) <- c("Log-Likelihood (variational approximation):")
  names(AICva) <- c("AIC (variational approximation):")
  names(LL) <- c("Log-Likelihood (G-H Quadrature correction):")
  names(AIC) <- c("AIC (G-H Quadrature correction):")
  
  out1 <-
    list(
      b = b,
      v = v,
      gamma = gamma_out,
      mu = mu,
      C = C,
      ugh.star = ugh.star,
      LL = LL,
      AIC = AIC,
      LLva = LLva,
      AICva = AICva,
      npar = npar
    )
  
  return(out1)
  
}



f_MELTA_g1_final <- function(Y, X, D, tol, maxiter, npoints, sw, 
                             b.init, v.init, gamma.init)
{
  
  # Function parameters:
  # Y = matrix of response variables
  # X = covariates
  # D = latent trait dimension
  # tol = tolerance level for convergence
  # maxiter = maximum number of iterations
  # npoints = quadrature points
  # sw = sampling weigths
  
  
  N <- nrow(Y)         # units
  L <- ncol(Y)         # items
  P <- ncol(X)         # covariates (including intercept)
  
  # Remove intercept from X for gamma (identifiability: gamma_0k = 0)
  X_gamma <- as.matrix(X[, -1, drop=FALSE])
  P_gamma <- ncol(X_gamma)
  dim_wh  <- D + 1 + P_gamma
  
  # Parameter initialization
  if(!is.null(gamma.init) && nrow(as.matrix(gamma.init)) == P)
    gamma.init <- as.matrix(gamma.init)[-1, , drop=FALSE]
  gamma <- gamma.init
  
  v <- v.init
  
  b=b.init
  
  
  # Variational approximation initialization
  
  xi <- matrix(20, N, L)
  sigma_xi <- 1 / (1 + exp(-xi))
  lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
  
  C <- matrix(0, D*D, N)
  mu <- matrix(0, N, D)
  
  
  # Iterative process
  
  ll <- -Inf
  diff <- 1
  iter <- 0
  cond <- TRUE
  
  # tol <- 10^-5
  
  lxi <- matrix(0, N, 1) # variational approx
  
  YY <- matrix(0, D * D, N)
  gam <- matrix(0, dim_wh, L)
  K <- array(diag(dim_wh), c(dim_wh, dim_wh, L))
  wh <- matrix(0, dim_wh, L) # matrix to store model parameters
  
  while (diff > tol & iter < maxiter)
  {
    iter <- iter + 1
    ll_old <- ll
    
    ## D=1 ####
    if (D == 1) {
      # Posterior Statistics
      
      C <-
        1 / (1 - 2 * rowSums(sweep(lambda_xi, MARGIN = 2, v ^ 2, `*`)))
      mu <-
        C * rowSums(sweep(
          Y - 0.5 + 2 * lambda_xi * (X_gamma%*%gamma) + 
            2 * sweep(lambda_xi, MARGIN = 2, b, `*`),
          MARGIN = 2, v, `*`))
      
      YY <- matrix(C + mu ^ 2, ncol = 1)
      
      # Variational Parameters (xi)
      
      xi <-  
        (YY %*% v^2) + (mu %*% (2 * b * v)) +
        (apply((t(apply(2 * (X_gamma%*%gamma),1,function(x) x* v))), 2, function(y) y*mu)) + 
        (2 * t(apply(X_gamma%*%gamma,1,function(xx)  xx * b))) + 
        (matrix(b^2, nrow = N, ncol = L, byrow = TRUE)) + 
        (matrix((X_gamma%*%gamma)^2, nrow = N, ncol = L, byrow = FALSE))
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      # Model Parameters
      
      gam[1, ] <- colSums(crossprod(sw * mu, Y - 0.5)) # v
      gam[2, ] <- crossprod(sw, Y - 0.5) # b
      gam[3:(P_gamma+2), ] <- t(apply(apply(X_gamma,2,function(x) x*sw), # gamma
                                      +           c(2), function(y) crossprod(y, Y-0.5)))
      
      
      K[1, 1, ] <- apply(sw * YY[,1] * lambda_xi, 2, sum) # vv
      
      K[2, 1, ] <- apply( # b v
        sw * mu * lambda_xi, c(2), sum)
      
      K[1, 2, ] <- K[2, 1, ] # v b
      
      K[3:(P_gamma+2), 1, ] <- apply(lambda_xi * sw * mu, # gamma v
                                     c(2), function(xx) crossprod(xx, X_gamma))
      
      K[1, 3:(P_gamma+2), ] <- K[3:(P_gamma+2), 1, ] # v gamma
      
      for(l in 1:L){
        
        K[2,2,l] <- apply(lambda_xi*sw,2,sum)[l]
        
        K[2, (3):(P_gamma+2), l] <- apply(lambda_xi * sw, # b gamma
                                          c(2), function(xx) crossprod(xx, X_gamma))[,l]
        
        K[(3):(P_gamma+2), 2, l] = t(K[2, (3):(P_gamma+2), l]) # gamma b
        
        valori_gg <- apply(array(apply(array(apply(lambda_xi * sw,
                                                   c(2), function(xx) xx*X_gamma),c(N,P_gamma,L)),
                                       c(3), function(xxx) xxx*X_gamma),c(N,P_gamma,L)), # gamma gamma
                           c(2,3), sum)[,l]
        if(P_gamma == 1){
          K[(D+2):dim_wh, (D+2):dim_wh, l] <- valori_gg
        } else {
          for(jj in seq_len(P_gamma)){
            K[(D+1+jj), (D+1+jj), l] <- valori_gg[jj]
          }
        }
        
        wh[, l] <- -ginv(2 * K[, , l]) %*% gam[, l] # parameter updates
      }
      
      
      v <- wh[1, ]
      v <- as.matrix(t(v))
      
      b <- wh[2, ]
      b <- as.matrix(t(b))
      
      gamma=wh[(D+2):dim_wh,]
      
      
      # log(P(Y|Z,X,xi))
      
      lxi <-
        (0.5 * log(C)) + (mu ^ 2 / (2 * C)) + rowSums(
          (log(sigma_xi)) - (0.5 * xi) - (lambda_xi * xi ^ 2) + 
            (sweep(lambda_xi, MARGIN = 2, b ^ 2, `*`)) +
            (lambda_xi * (X_gamma%*%gamma)^2) +
            (2 * sweep(lambda_xi * (X_gamma%*%gamma), MARGIN=2, b, `*`)) + 
            (sweep(Y - 0.5, MARGIN = 2, b, `*`)) + 
            ((Y-0.5) * (X_gamma%*%gamma)))
    } # end D=1
    
    ## D=2 ####
    if (D == 2) {
      
      # Latent Posterior Statistics
      
      C <-
        apply(lambda_xi, 1, function(x)
          solve(diag(D) - 2 * crossprod(x * t(v), t(v))))
      
      mu <-
        t(apply(rbind(C, 
                      tcrossprod(v, Y - 0.5 + 2 * sweep(lambda_xi, MARGIN = 2, b, `*`)+
                                   2 * (lambda_xi * (X_gamma%*%gamma)))), 
                2, function(x) matrix(x[1:4], nrow = D) %*% x[-(1:4)]))
      
      # Variational Parameters (xi)
      
      YY <- C + apply(mu, 1, tcrossprod)
      
      b=c(b)
      xi <-
        t(
          apply(YY, 2, function(x)
            rowSums(crossprod(v, matrix(x, ncol = D)) * t(v))) + 
            tcrossprod(2 * b * t(v), mu) + 
            2 * t(X_gamma%*%gamma) * tcrossprod(t(v),mu) + 
            2 * t(X_gamma%*%gamma) * b + 
            + matrix(
              b ^ 2,
              nrow = L,
              ncol = N,
              byrow = FALSE
            ) + 
            t((X_gamma%*%gamma)^2)
        )
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      
      # Model Parameters 
      
      gam[3, ] <- crossprod(sw, Y - 0.5) # b
      
      gam[4:(P_gamma+3), ] <- t(apply(sw*X_gamma, # gamma
                                      c(2), function(y) crossprod(y, Y-0.5)))
      
      aa <- matrix(sw, N, D) * mu
      bb <- t(matrix(sw, N, D * D)) * YY
      
      K[3, 3, ] <- # b b
        apply(sw * lambda_xi, c(2), sum)
      
      kk <- 0
      vgamma <- 0
      
      K[3, 1:2, ] <- crossprod(aa, lambda_xi) # b v
      K[1:2, 3, ] <- K[3, (1:2), ] # v b
      
      kk <- kk + bb %*% lambda_xi
      
      vgamma <- vgamma + apply(array(apply(array(apply(aa,2,function(x) x*lambda_xi),c(N,L,D)),
                                           c(2,3), function(xx) xx*X_gamma),c(N,P_gamma,L,D)),
                               c(2,3,4), sum)
      
      K[1:2, 1:2, ] <- kk # v v
      
      for (l in 1:L)
      {
        gam[1:D, l] <- apply(aa * (Y[, l] - 0.5), 2, sum) # v
        
        K[4:(P_gamma+3), 1:D, l] = vgamma[,l,] # gamma v
        K[1:D, 4:(P_gamma+3), l] = t(K[4:(P_gamma+3), 1:D, l])  # v gamma 
        
        K[3, 4:(P_gamma+3), l] <- apply(lambda_xi * sw, # b gamma
                                        c(2), function(xx) crossprod(xx, X_gamma))[,l]
        K[4:(P_gamma+3), 3, l] = t(K[3, 4:(P_gamma+3), l]) # gamma b
        
        valori_gg <- apply(array(apply(array(apply(lambda_xi * sw,
                                                   c(2), function(xx) xx*X_gamma),c(N,P_gamma,L)),
                                       c(3), function(xxx) xxx*X_gamma),c(N,P_gamma,L)), # gamma gamma
                           c(2,3), sum)[,l]
        if(P_gamma == 1){
          K[(D+2):dim_wh, (D+2):dim_wh, l] <- valori_gg
        } else {
          for(jj in seq_len(P_gamma)){
            K[(D+1+jj), (D+1+jj), l] <- valori_gg[jj]
          }
        }
        
        wh[, l] <- -ginv(2 * K[, , l]) %*% gam[, l]
      }
      
      v <- wh[1:D, ]
      
      b <- wh[3, ]
      
      gamma=wh[4:(P_gamma+3),]
      
      
      # log(P(Y|Z,X,xi))
      
      detC <- C[1, ] * C[4, ] - C[3, ] * C[2, ]
      
      lxi <-
        (0.5 * log(detC)) + (0.5 * apply(rbind(C[4, ] / detC, -C[2, ] / detC, -C[3, ] /
                                                 detC, C[1, ] / detC, t(mu)), 2, function(x)
                                                   t((x[-(1:4)])) %*% matrix(x[1:4], nrow = D) %*% (x[-(1:4)]))) + 
        rowSums(
          (log(sigma_xi)) - (0.5 * xi) - (lambda_xi * xi^2) + 
            (sweep(lambda_xi, MARGIN = 2, b ^ 2, `*`)) + 
            (lambda_xi * (X_gamma%*%gamma)^2) +
            (2 * sweep(lambda_xi * (X_gamma%*%gamma), MARGIN=2, b, `*`)) + 
            (sweep(Y - 0.5, MARGIN = 2, b, `*`)) +
            ((Y-0.5) * (X_gamma%*%gamma))
        )
    }# end D=2
    
    ## D>2 ####
    if (D > 2) {
      
      # Latent Posterior Statistics
      
      C <-
        apply(lambda_xi, 1, function(x)
          solve(diag(D) - 2 * crossprod(x * t(v), t(v))))
      
      mu <-
        t(apply(rbind(C, v %*% t(
          Y - 0.5 + 2 * sweep(lambda_xi, MARGIN = 2, b, `*`) +
            2 * lambda_xi * (X_gamma%*%gamma)
        )), 2, function(x)
          matrix(x[1:(D * D)], nrow = D) %*% x[-(1:(D * D))]))
      
      # Variational Parameters (xi)
      
      YY <- C + apply(mu, 1, tcrossprod)
      
      b=c(b)
      xi <-
        t(
          apply(YY, 2, function(x)
            rowSums((t(v) %*% matrix(x, ncol = D)) * t(v))) + 
            (2 * b * t(v)) %*% t(mu) + 
            2 * t(X_gamma%*%gamma) * tcrossprod(t(v),mu) + 
            2 * t(X_gamma%*%gamma) * b + 
            + matrix(
              b ^ 2,
              nrow = L,
              ncol = N,
              byrow = FALSE
            ) + 
            t((X_gamma%*%gamma)^2)
        )
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      
      # Model Parameters 
      
      gam[(D + 1), ] <- t(sw) %*% (Y - 0.5) # b
      
      gam[(D+2):dim_wh, ] <- t(apply(sw*X_gamma, # gamma
                                     c(2), function(y) crossprod(y, Y-0.5)))
      
      aa <- matrix(sw, N, D) * mu
      bb <- t(matrix(sw, N, D * D)) * YY
      
      K[(D + 1), (D + 1), ] <- # b b
        apply(sw * lambda_xi, c(2), sum)
      
      kk <- 0
      vgamma <- 0
      
      K[D + 1, 1:D, ] <- crossprod(aa, lambda_xi) # b v
      K[1:D, D + 1, ] <- K[D + 1, 1:D, ] # v b 
      kk <- kk + bb %*% lambda_xi
      vgamma <- vgamma + apply(array(apply(array(apply(aa,2,function(x) x*lambda_xi),c(N,L,D)),
                                           c(2,3), function(xx) xx*X_gamma),c(N,P_gamma,L,D)),
                               c(2,3,4), sum)
      
      K[1:D, 1:D, ] <- kk # v v
      
      for (l in 1:L)	{
        
        gam[1:D, l] <- apply(aa * (Y[, l] - 0.5), 2, sum) # v
        
        K[(D+2):dim_wh, 1:D, l] = vgamma[,l,] # gamma v
        
        K[1:D, (D+2):dim_wh, l] = t(K[(D+2):dim_wh, 1:D, l])  # v gamma 
        
        K[(D+1), (D+2):dim_wh, l] <- apply(lambda_xi * sw, # b gamma
                                           c(2), function(xx) crossprod(xx, X_gamma))[,l]
        
        K[(D+2):dim_wh, (D+1), l] = t(K[(D+1), (D+2):dim_wh, l]) # gamma b
        
        valori_gg <- apply(array(apply(array(apply(lambda_xi * sw,
                                                   c(2), function(xx) xx*X_gamma),c(N,P_gamma,L)),
                                       c(3), function(xxx) xxx*X_gamma),c(N,P_gamma,L)), # gamma gamma
                           c(2,3), sum)[,l]
        if(P_gamma == 1){
          K[(D+2):dim_wh, (D+2):dim_wh, l] <- valori_gg
        } else {
          for(jj in seq_len(P_gamma)){
            K[(D+1+jj), (D+1+jj), l] <- valori_gg[jj]
          }
        }
        
        wh[, l] <- -ginv(2 * K[, , l]) %*% gam[, l]
      }
      
      
      v <- wh[1:D, ]
      v <- as.matrix(v)
      
      b <- wh[D + 1, ]
      
      gamma=wh[(D+2):dim_wh,]
      
      
      # log(P(Y|U,Z,X,xi))
      
      detC <- apply(C, 2, function(x)
        det(matrix(x, D, D)))
      
      lxi <-
        (0.5 * log(detC)) + (0.5 * apply(rbind(C, t(mu)), 2, function(x)
          t((x[-(1:(D * D))])) %*% solve(matrix(x[1:(D * D)], nrow = D)) %*% (x[-(
            1:(D * D))]))) + rowSums(
              (log(sigma_xi)) - (0.5 * xi) - (lambda_xi *
                                                xi^2) + sweep(lambda_xi, MARGIN = 2, b ^ 2, `*`) + 
                (lambda_xi * (X_gamma%*%gamma)^2) +
                (2 * sweep(lambda_xi * (X_gamma%*%gamma), MARGIN=2, b, `*`)) + 
                (sweep(Y - 0.5, MARGIN = 2, b, `*`)) +
                ((Y-0.5) * (X_gamma%*%gamma))
            )
      
    } # end D > 2
    
    
    # Log-likelihood
    
    ll <- sum(sw * lxi)
    
    
    # Stopping Criteria
    
    diff <- sum(abs(ll - ll_old))
    
    if(sum((ll - ll_old))<0) print(paste(iter, ll, sum((ll - ll_old)), "aaaaaa"))
    else print(paste(iter, ll, sum((ll - ll_old))))
    
  } # end while
  
  
  # Correction to the log-likelihood: Gauss-Hermite Quadrature
  
  Qd <- npoints ^ D       
  GaHer <- glmmML::ghq(npoints, FALSE)                     
  ugh <- as.matrix(expand.grid(rep(list(GaHer$zeros), D))) 
  ugh.star <- sqrt(2)*ugh                                 
  p.gh = as.matrix(expand.grid(rep(list(GaHer$weights), times = D)))	
  Ww =  (2)^(D/2) * exp(apply(ugh, 1, crossprod)) * apply(p.gh, 1, prod)
  W.prod =  apply(p.gh, 1, prod)
  Phi <- apply(ugh.star, 1, mvtnorm::dmvnorm)               
  
  fxy <- matrix(0, N, Qd)
  for(i in 1:N){
    if(D>1){
      Agh <- t(tcrossprod(t(v), ugh.star)) + b + 
        (X_gamma%*%gamma)[i,]
    } else {
      Agh <- t(apply(tcrossprod(t(v), ugh.star), c(2), function(x) x + b + 
                       (X_gamma%*%gamma)[i,]))
    }
    pgh <- 1 / (1 + exp(-Agh))
    fxy[i, ] <-
      exp(tcrossprod(Y[i,], log(pgh)) + tcrossprod(1 - Y[i,], log(1 - pgh)))
  }
  
  LLva <- ll
  npar = (P_gamma * L) + (L) + ((L * D) - (D * (D - 1)/2))
  AICva <- -2 * LLva + 2 * npar
  
  fxy[!is.finite(fxy)] <- 0
  llGH1 <- rowSums(t(apply(fxy,c(1),function(x) x*Ww*Phi)))
  LL = sum(sw * log(llGH1))
  AIC <- -2 * LL + 2 * npar
  
  b=matrix(b,ncol=L)
  colnames(v) <- colnames(v, do.NULL = FALSE, prefix = "Item ")
  
  rownames(v) <- rownames(v, do.NULL = FALSE, prefix = "Dim ")
  colnames(v) <- colnames(v, do.NULL = FALSE, prefix = "Item ")
  
  gamma_out <- rbind(rep(0, L), gamma)
  colnames(gamma_out) <- colnames(Y)
  rownames(gamma_out) <- colnames(X)
  
  names(LLva) <- c("Log-Likelihood (variational approximation):")
  names(AICva) <- c("AIC (variational approximation):")
  names(LL) <- c("Log-Likelihood (G-H Quadrature correction):")
  names(AIC) <- c("AIC (G-H Quadrature correction):")
  
  out1 <-
    list(
      b = b,
      v = v,
      gamma = gamma_out,
      mu = mu,
      C = C,
      ugh.star = ugh.star,
      LL = LL,
      AIC = AIC,
      LLva = LLva,
      AICva = AICva, 
      npar = npar
    )
  
  return(out1)
  
}


f_MELTA_methods <- function(Y, X, G, D, nstarts, tol, maxiter, npoints, sw)
{
  if (D == 0) {
    if (any(G == 1)) {
      out <- f_D0_nstarts(Y, X, G, D, nstarts, tol, maxiter, sw)
    } else{
      if (length(G) == 1) {
        out <- f_D0_nstarts(Y, X, G, D, nstarts, tol, maxiter, sw)
      } else{
        out <- vector("list", length(G) + 1)
        names(out) <- c(paste('G', G, sep = '='), 'AIC')
        i <- 0
        for (g in G) {
          i <- i + 1
          out[[i]] <- f_D0_nstarts(Y, X, g, D, nstarts, tol, maxiter, sw)
        }
        out[[length(G) + 1]] <- tableAIC(out)
      }
    }
  } else {
    if (D > 0 && G == 1){
      if(length(D) == 1){ 
        out <- f_MELTA_g1_nstarts_vfix(Y, X, D, nstarts, tol, maxiter, npoints, sw)
        out$eta <- 1
      }else{
        out<-vector("list", length(D) + 1)
        names(out) <- c(paste('Dim', D, sep = '='), 'AIC')
        i<-0
        for(diy in D){
          i <- i + 1
          out[[i]] <- f_MELTA_g1_nstarts_vfix(Y, X, diy, nstarts, tol, maxiter, npoints, sw)
          out[[i]]$eta <- 1
        }
        cat('AIC results',"\n \n")
        out[[length(D) + 1]]<-tableAIC(out)
      }
    }
    if (D > 0 && G > 1) {
      out <- f_MELTA_nstarts_vfix(Y, X, G, D, nstarts, tol, maxiter, npoints, sw)
      class(out) <- c("MELTA")
    }
  }
  class(out) <- c("MELTA")
  return(out)
}

MELTA <- function(Y, X, G, D, nstarts = 10, tol = 10^-4, maxiter = 100, npoints = 7, sw){
  if (any(G < 1)) {
    print("Specify G > 0!")
    return("Specify G > 0!")
  }
  
  if (any(D < 0)) {
    print("Specify D >= 0!")
    return("Specify D >= 0!")
  }
  
  if (length(D) == 1 && length(G) == 1) {
    out <- f_MELTA_methods(Y, X, G, D, nstarts, tol, maxiter, npoints, sw)
  } else{
    out <- vector("list", length(D) * length(G) + 3)
    names(out) <- c(t(outer(
      paste('G=', G, sep = ''),
      paste('D=', D, sep = ''),
      paste,
      sep = ','
    )),
    'AIC', 'LL', 'LLva')
    
    aictab <- matrix(0, length(D), length(G))
    lltab <- matrix(0, length(D), length(G))
    
    rownames(aictab) <- paste('D=', D, sep = '')
    colnames(aictab) <- paste('G=', G, sep = '')
    
    rownames(lltab) <- paste('D=', D, sep = '')
    colnames(lltab) <- paste('G=', G, sep = '')
    
    llvatab <- lltab
    
    i <- 0
    for (g in G){
      for (diy in D) {
        i <- i + 1
        
        out[[i]] <- f_MELTA_methods(Y, X, g, diy, nstarts, tol, maxiter, npoints, sw)
        print("ciao")
        aictab[i] <- out[[i]]$AIC
        lltab[i] <- out[[i]]$LL
        
        if (diy == 0) {
          llvatab[i] <- out[[i]]$LL
        } else{
          llvatab[i] <- out[[i]]$LLva
        }
        
        out[[length(G) * length(D) + 1]] <- ResTable(aictab, restype = 'AIC')
        out[[length(G) * length(D) + 2]] <- ResTable(lltab, restype = 'll')
        out[[length(G) * length(D) + 3]] <- ResTable(llvatab, restype = 'llva')
      }
    }
    class(out) <- "mMELTA"
  }
  
  out
}


ResTable <- function(aicG, restype)
{
  if (restype == 'll') {
    resAIC <- vector('list', 1)
    names(resAIC) <- 'Table of LL (G-H Quadrature correction)'
    resAIC[[1]] <- aicG
  }
  
  if (restype == 'llva') {
    resAIC <- vector('list', 1)
    names(resAIC) <- 'Table of LL (variational approximation)'
    resAIC[[1]] <- aicG
  }
  
  if (restype == 'AIC') {
    resAIC <- vector('list', 2)
    names(resAIC) <- c('Table of BIC Results', 'Model Selection')
    resAIC[[1]] <- aicG
    resAIC[[2]] <-
      paste(colnames(aicG)[t(aicG == min(aicG)) * seq(1:ncol(aicG))], rownames(aicG)[(aicG ==
                                                                                        min(aicG)) * seq(1:nrow(aicG))], sep = ', ')
    names(resAIC[[2]]) <- 'Model with lower AIC:'
  }
  
  resAIC
}


f_MELTA_nstarts_vfix <- function(Y, X, G, D, nstarts, tol, maxiter, npoints, sw)
{
  out <- f_MELTA_vfix(Y, X, G, D, tol, maxiter, npoints, sw)
  
  if(nstarts > 1){ 
    for(i in 2:nstarts){
      out1 <- f_MELTA_vfix(Y, X, G, D, tol, maxiter, npoints, sw)
      if(out1$LL > out$LL) out <- out1
    }
  }
  
  return(out)
}


f_MELTA_g1_nstarts_vfix <- function(Y, X, D, nstarts, tol, maxiter, npoints, sw)
{
  out <- f_MELTA_g1(Y, X, D, tol, maxiter, npoints, sw)
  
  if(nstarts > 1){ 
    for(i in 2:nstarts){
      out1 <- f_MELTA_g1(Y, X, D, tol, maxiter, npoints, sw)
      if(out1$LL > out$LL) out <- out1
    }
  }
  
  return(out)
}


tableAIC <- function(out){
  
  lout <- length(out) - 1
  aicG <- numeric(lout)
  names(aicG) <- names(out[- (lout + 1)])
  
  for(i in 1:lout) aicG[i] <- out[[i]]$AIC
  
  resAIC <- vector('list', 2)
  names(resAIC) <- c('Model Selection', 'Table of AIC Results')
  resAIC[[1]] <- names(aicG)[aicG == min(aicG)]
  names(resAIC[[1]]) <- 'Model with lower AIC:'
  resAIC[[2]] <- aicG
  resAIC
}



print.MELTA <- function(x){
  stopifnot(inherits(x, 'MELTA'))
  cat("beta:\n")
  print(x$beta)
  cat("b:\n")
  print(x$b)
  cat("\nv:\n")
  print(x$v)
  cat("\neta:\n")
  print(x$eta)
  cat("\nz:\n")
  print(x$z)
  cat("\nz1:\n")
  print(x$z1)
  cat("\n")
  print(x$LL)
  print(x$AIC)
}



print.mMELTA <- function(x){
  cat("Log-Likelihood:\n")
  print(x$LL$`Table of LL (G-H Quadrature correction)`)
  cat("AIC:\n")
  print(x$AIC$`Table of AIC Results`)
  print(x$AIC$`Model Selection`)
}
