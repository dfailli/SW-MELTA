cls<- function(lab) 
{
  gr <- sort(unique(lab))
  z <- matrix(0, length(lab), length(gr))
  for (i in 1:length(gr)) z[lab == gr[i], i] <- 1
  z
}

f_MELTA_vfix <- function(Y, X, G, D, tol, maxiter, npoints, sw)
{
  
  N <- nrow(Y)         # units
  L <- ncol(Y)         # items
  P <- ncol(X)         # covariates (including intercept)
  
  # X without intercept for gamma (identifiability)
  X_gamma  <- as.matrix(X[, -1, drop=FALSE])
  P_gamma  <- ncol(X_gamma)  # P - 1
  
  # Total parameter dimension per item: D (v) + G (b) + P_gamma (gamma)
  dim_wh <- D + G + P_gamma
  
  # Parameter initialization
  beta  <- matrix(rep(0, P*(G-1)), ncol=G-1)
  gamma <- matrix(rep(0, P_gamma*L), nrow=P_gamma, ncol=L)
  
  exb <- exp(X %*% beta)
  eta <- cbind(1,exb)/(rowSums(exb)+1)
  
  z <- matrix(NA, nrow=N, ncol=G)
  for(i in 1:N) z[i,] <- t(rmultinom(1, size=1, prob=eta[i,]))
  
  v <- matrix(rnorm(L * D), D, L)
  b <- matrix(rnorm(L*G, 0, 0.000001), G, L)
  ord.b <- order(apply(b,1,mean))
  b <- b[ord.b,]
  
  # Variational approximation initialization
  xi       <- array(20, c(N, L, G))
  sigma_xi <- 1 / (1 + exp(-xi))
  lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
  
  C  <- array(0, c(D*D, N, G))
  mu <- array(0, c(N, D, G))
  
  ll   <- -Inf
  diff <- 1
  iter <- 0
  tol  <- 10^-4
  
  print(c(beta))
  
  lxi <- matrix(0, N, G)
  YY  <- array(0, c(D * D, N, G))
  
  # gam and K now use dim_wh instead of D+G+P
  gam <- matrix(0, dim_wh, L)
  K   <- array(diag(dim_wh), c(dim_wh, dim_wh, L))
  wh  <- matrix(0, dim_wh, L)
  
  while (diff > tol & iter < maxiter)
  {
    iter   <- iter + 1
    ll_old <- ll
    
    ## D=1 ####
    if (D == 1) {
      for (g in 1:G)
      {
        C[,,g] <- 1 / (1 - 2 * rowSums(sweep(lambda_xi[,,g], MARGIN=2, v^2, `*`)))
        mu[,,g] <- C[,,g] * rowSums(sweep(
          Y - 0.5 + 2 * lambda_xi[,,g] * (X_gamma %*% gamma) +
            2 * sweep(lambda_xi[,,g], MARGIN=2, b[g,], `*`),
          MARGIN=2, v, `*`))
        
        YY[,,g] <- matrix(C[,,g] + mu[,,g]^2, ncol=1)
        
        xi[,,g] <-
          (YY[,,g] %*% v^2) + (mu[,,g] %*% (2 * b[g,] * v)) +
          (apply((t(apply(2 * (X_gamma %*% gamma), 1, function(x) x * v))), 2, function(y) y * mu[,,g])) +
          (2 * t(apply(X_gamma %*% gamma, 1, function(xx) xx * b[g,]))) +
          (matrix(b[g,]^2, nrow=N, ncol=L, byrow=TRUE)) +
          (matrix((X_gamma %*% gamma)^2, nrow=N, ncol=L, byrow=FALSE))
      }
      
      xi        <- sqrt(xi)
      sigma_xi  <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      # Model Parameters
      gam[1, ]                    <- colSums(crossprod(sw * z * mu[,1,], Y - 0.5))  # v
      gam[2:(1+G), ]              <- crossprod(sw * z, Y - 0.5)                      # b
      gam[(D+G+1):(dim_wh), ] <- t(apply(apply(array(apply(sw*z, 2, function(x) x*X_gamma), c(N,P_gamma,G)),
                                               c(2,3), function(y) crossprod(y, Y-0.5)),
                                         c(1,2), function(xy) sum(xy)))           # gamma
      
      K[1, 1, ] <- apply(aperm(array(sw * z * YY[1,,], c(N,G,L)), c(1,3,2)) * lambda_xi, 2, sum) # vv
      
      K[2:(G+1), 1, ] <- t(apply(aperm(array(sw * z * mu[,1,], c(N,G,L)), c(1,3,2)) * lambda_xi, c(2,3), sum)) # bv
      K[1, 2:(G+1), ] <- K[2:(G+1), 1, ]  # vb
      
      # gamma v
      tmp_gv <- apply(array(apply(array(apply(lambda_xi, 2, function(x) x*(sw * z * mu[,1,])), c(N,G,L)),
                                  c(2,3), function(xx) crossprod(xx, X_gamma)),c(G,L,P_gamma)),
                      c(2,3), sum)
      if(P_gamma == 1) tmp_gv <- t(tmp_gv)
      K[(D+G+1):(dim_wh), 1, ] <- tmp_gv
      K[1, (D+G+1):(dim_wh), ] <- K[(D+G+1):(dim_wh), 1, ]  # v gamma
      
      for(l in 1:L){
        diag(K[2:(G+1), 2:(G+1), l]) <- apply(array(apply(lambda_xi, 2, function(x) sw*z*x), c(N,G,L)),
                                              c(2,3), sum)[,l]  # bb
        
        # b gamma
        tmp_bg <- t(aperm(array(apply(array(apply(lambda_xi, 2, function(x) sw*z*x), c(N,G,L)),
                                      c(2,3), function(xx) crossprod(xx, X_gamma)), c(G,L,P_gamma)),c(1,3,2))[,,l])
        K[2:(G+1), (D+G+1):(dim_wh), l] <- tmp_bg
        K[(D+G+1):(dim_wh), 2:(G+1), l] <- t(tmp_bg)  # gamma b
        
        # gamma gamma
        idx_gamma <- (D + G + 1):dim_wh
        valori_diag <- apply(array(apply(array(apply(array(apply(lambda_xi, 2, function(x) sw*z*x), c(N,G,L)),
                                                     c(2,3), function(xx) xx*X_gamma), c(N,P_gamma,G,L)),
                                         c(3,4), function(xxx) xxx*X_gamma), c(N,P_gamma,G,L)),
                             c(2,4), sum)[, l]
        if (length(idx_gamma) == 1) {
          # Se c'è solo 1 covariata (matrice 1x1), è una cella singola. Assegnazione diretta.
          K[idx_gamma, idx_gamma, l] <- valori_diag
        } else {
          # Se ci sono più covariate, creiamo una matrice di indici bidimensionali 
          # per toccare *solo* gli elementi della diagonale di quel blocco.
          indici_diagonale <- cbind(idx_gamma, idx_gamma)
          K[, , l][indici_diagonale] <- valori_diag
        }
        
        # diag(K[(D+G+1):(dim_wh), (D+G+1):(dim_wh), l]) <-
        #   apply(array(apply(array(apply(array(apply(lambda_xi, 2, function(x) sw*z*x), c(N,G,L)),
        #                                 c(2,3), function(xx) xx*X_gamma), c(N,P_gamma,G,L)),
        #                     c(3,4), function(xxx) xxx*X_gamma), c(N,P_gamma,G,L)),
        #         c(2,4), sum)[,l]
        
        wh[, l] <- -ginv(2 * K[,,l]) %*% gam[,l]
      }
      
      v     <- as.matrix(t(wh[1:D, ]))
      for(g in 1:G) b[g,] <- wh[D+g, ]
      ord.b <- order(apply(b,1,mean))
      b     <- b[ord.b,]
      gamma <- wh[(D+G+1):(dim_wh), , drop=FALSE]
      
      # log(P(Y|Z,X,xi))
      for(g in 1:G){
        lxi[,g] <-
          (0.5 * log(C[1,,g])) + (mu[,1,g]^2 / (2 * C[1,,g])) + rowSums(
            (log(sigma_xi[,,g])) - (0.5 * xi[,,g]) - (lambda_xi[,,g] * xi[,,g]^2) +
              (sweep(lambda_xi[,,g], MARGIN=2, b[g,]^2, `*`)) +
              (lambda_xi[,,g] * (X_gamma %*% gamma)^2) +
              (2 * sweep(lambda_xi[,,g] * (X_gamma %*% gamma), MARGIN=2, b[g,], `*`)) +
              (sweep(Y - 0.5, MARGIN=2, b[g,], `*`)) +
              ((Y-0.5) * (X_gamma %*% gamma)))
      }
    } # end D=1
    
    ## D=2 ####
    if (D == 2) {
      for (g in 1:G)
      {
        C[,,g] <- apply(lambda_xi[,,g], 1, function(x)
          solve(diag(D) - 2 * crossprod(x * t(v), t(v))))
        
        mu[,,g] <- t(apply(rbind(C[,,g],
                                 tcrossprod(v, Y - 0.5 + 2 * sweep(lambda_xi[,,g], MARGIN=2, b[g,], `*`) +
                                              2 * (lambda_xi[,,g] * (X_gamma %*% gamma)))),
                           2, function(x) matrix(x[1:4], nrow=D) %*% x[-(1:4)]))
        
        YY[,,g] <- C[,,g] + apply(mu[,,g], 1, tcrossprod)
        
        xi[,,g] <- t(
          apply(YY[,,g], 2, function(x) rowSums(crossprod(v, matrix(x, ncol=D)) * t(v))) +
            tcrossprod(2 * b[g,] * t(v), mu[,,g]) +
            2 * t(X_gamma %*% gamma) * tcrossprod(t(v), mu[,,g]) +
            2 * t(X_gamma %*% gamma) * b[g,] +
            matrix(b[g,]^2, nrow=L, ncol=N, byrow=FALSE) +
            t((X_gamma %*% gamma)^2)
        )
      }
      
      xi        <- sqrt(xi)
      sigma_xi  <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      gam[3:(2+G), ]              <- crossprod(sw*z, Y - 0.5)  # b
      gam[(D+G+1):(dim_wh), ] <- t(apply(apply(array(apply(sw*z, 2, function(x) x*X_gamma), c(N,P_gamma,G)),
                                               c(2,3), function(y) crossprod(y, Y-0.5)),
                                         c(1,2), function(xy) sum(xy)))  # gamma
      
      aa <- aperm(array(sw*z, c(N,G,D)), c(1,3,2)) * mu
      bb_arr <- aperm(array(sw*z, c(N,G,D*D)), c(3,1,2)) * YY
      
      K[3:(2+G), 3:(2+G), ] <- apply(apply(aperm(array(sw*z, c(N,G,L)), c(1,3,2)) * lambda_xi,
                                           c(2,3), sum), 1, diag)  # bb
      
      kk     <- 0
      vgamma <- 0
      for(g in 1:G){
        K[2+g, 1:2, ] <- crossprod(aa[,,g], lambda_xi[,,g])  # bv
        K[1:2, 2+g, ] <- K[2+g, 1:2, ]                       # vb
        kk     <- kk + bb_arr[,,g] %*% lambda_xi[,,g]
        vgamma <- vgamma + apply(array(apply(array(apply(aa[,,g], 2, function(x) x*lambda_xi[,,g]), c(N,L,D)),
                                             c(2,3), function(xx) xx*X_gamma), c(N,P_gamma,L,D)),
                                 c(2,3,4), sum)
      }
      K[1:2, 1:2, ] <- kk  # vv
      
      for(l in 1:L){
        gam[1:D, l] <- apply(aa * (Y[,l] - 0.5), 2, sum)  # v
        
        K[(D+G+1):(dim_wh), 1:D, l] <- vgamma[,l,]          # gamma v
        K[1:D, (D+G+1):(dim_wh), l] <- t(K[(D+G+1):(dim_wh), 1:D, l])  # v gamma
        
        # b gamma
        tmp_bg <- t(aperm(array(apply(array(apply(lambda_xi, 2, function(x) sw*z*x), c(N,G,L)),
                                      c(2,3), function(xx) crossprod(xx, X_gamma)),c(G, L, P_gamma)),c(1,3,2))[,,l])
        K[(D+1):(G+D), (D+G+1):(dim_wh), l] <- tmp_bg
        K[(D+G+1):(dim_wh), (D+1):(G+D), l] <- t(tmp_bg)  # gamma b
        
        # gamma gamma
        idx_gamma <- (D + G + 1):dim_wh
        valori_diag <- apply(array(apply(array(apply(array(apply(lambda_xi, 2, function(x) sw*z*x), c(N,G,L)),
                                                     c(2,3), function(xx) xx*X_gamma), c(N,P_gamma,G,L)),
                                         c(3,4), function(xxx) xxx*X_gamma), c(N,P_gamma,G,L)),
                             c(2,4), sum)[, l]
        if (length(idx_gamma) == 1) {
          # Se la matrice è 1x1, sovrascriviamo direttamente l'unica cella interessata
          K[idx_gamma, idx_gamma, l] <- valori_diag
        } else {
          # Se ci sono più covariate, usiamo l'indicizzazione lineare coordinata per coordinata
          indici_diagonale <- cbind(idx_gamma, idx_gamma)
          K[, , l][indici_diagonale] <- valori_diag
        }
        
        wh[,l] <- -ginv(2 * K[,,l]) %*% gam[,l]
      }
      
      v <- wh[1:D, ]
      for(g in 1:G) b[g,] <- wh[D+g, ]
      ord.b <- order(apply(b,1,mean))
      b     <- b[ord.b,]
      gamma <- wh[(D+G+1):(dim_wh), , drop=FALSE]
      
      for(g in 1:G){
        detC <- C[1,,g] * C[4,,g] - C[3,,g] * C[2,,g]
        lxi[,g] <-
          (0.5 * log(detC)) +
          (0.5 * apply(rbind(C[4,,g]/detC, -C[2,,g]/detC, -C[3,,g]/detC, C[1,,g]/detC, t(mu[,,g])),
                       2, function(x) t(x[-(1:4)]) %*% matrix(x[1:4], nrow=D) %*% x[-(1:4)])) +
          rowSums(
            (log(sigma_xi[,,g])) - (0.5 * xi[,,g]) - (lambda_xi[,,g] * xi[,,g]^2) +
              (sweep(lambda_xi[,,g], MARGIN=2, b[g,]^2, `*`)) +
              (lambda_xi[,,g] * (X_gamma %*% gamma)^2) +
              (2 * sweep(lambda_xi[,,g] * (X_gamma %*% gamma), MARGIN=2, b[g,], `*`)) +
              (sweep(Y - 0.5, MARGIN=2, b[g,], `*`)) +
              ((Y-0.5) * (X_gamma %*% gamma)))
      }
    } # end D=2
    
    
    ## D>2 ####
    if (D > 2) {
      for(g in 1:G){
        C[,,g] <- apply(lambda_xi[,,g], 1, function(x)
          solve(diag(D) - 2 * crossprod(x * t(v), t(v))))
        
        mu[,,g] <- t(apply(rbind(C[,,g], v %*% t(
          Y - 0.5 + 2 * sweep(lambda_xi[,,g], MARGIN=2, b[g,], `*`) +
            2 * lambda_xi[,,g] * (X_gamma %*% gamma))),
          2, function(x) matrix(x[1:(D*D)], nrow=D) %*% x[-(1:(D*D))]))
        
        YY[,,g] <- C[,,g] + apply(mu[,,g], 1, tcrossprod)
        
        xi[,,g] <- t(
          apply(YY[,,g], 2, function(x) rowSums((t(v) %*% matrix(x, ncol=D)) * t(v))) +
            (2 * b[g,] * t(v)) %*% t(mu[,,g]) +
            2 * t(X_gamma %*% gamma) * tcrossprod(t(v), mu[,,g]) +
            2 * t(X_gamma %*% gamma) * b[g,] +
            matrix(b[g,]^2, nrow=L, ncol=N, byrow=FALSE) +
            t((X_gamma %*% gamma)^2)
        )
      }
      
      xi        <- sqrt(xi)
      sigma_xi  <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      gam[(D+1):(D+G), ]          <- t(sw*z) %*% (Y - 0.5)  # b
      gam[(D+G+1):(dim_wh), ] <- t(apply(apply(array(apply(sw*z, 2, function(x) x*X_gamma), c(N,P_gamma,G)),
                                               c(2,3), function(y) crossprod(y, Y-0.5)),
                                         c(1,2), function(xy) sum(xy)))  # gamma
      
      aa     <- aperm(array(sw*z, c(N,G,D)), c(1,3,2)) * mu
      bb_arr <- aperm(array(sw*z, c(N,G,D*D)), c(3,1,2)) * YY
      
      K[(D+1):(D+G), (D+1):(D+G), ] <- apply(apply(aperm(array(sw*z, c(N,G,L)), c(1,3,2)) * lambda_xi,
                                                   c(2,3), sum), 1, diag)  # bb
      
      kk     <- 0
      vgamma <- 0
      for(g in 1:G){
        K[D+g, 1:D, ] <- crossprod(aa[,,g], lambda_xi[,,g])  # bv
        K[1:D, D+g, ] <- K[D+g, 1:D, ]                       # vb
        kk     <- kk + bb_arr[,,g] %*% lambda_xi[,,g]
        vgamma <- vgamma + apply(array(apply(array(apply(aa[,,g], 2, function(x) x*lambda_xi[,,g]), c(N,L,D)),
                                             c(2,3), function(xx) xx*X_gamma), c(N,P_gamma,L,D)),
                                 c(2,3,4), sum)
      }
      K[1:D, 1:D, ] <- kk  # vv
      
      for(l in 1:L){
        gam[1:D, l] <- apply(aa * (Y[,l] - 0.5), 2, sum)  # v
        
        K[(D+G+1):(dim_wh), 1:D, l] <- vgamma[,l,]          # gamma v
        K[1:D, (D+G+1):(dim_wh), l] <- t(K[(D+G+1):(dim_wh), 1:D, l])  # v gamma
        
        # b gamma
        tmp_bg <- t(aperm(array(apply(array(apply(lambda_xi, 2, function(x) sw*z*x), c(N,G,L)),
                                      c(2,3), function(xx) crossprod(xx, X_gamma)),c(G,L,P_gamma)),c(1,3,2))[,,l])
        K[(D+1):(G+D), (D+G+1):(dim_wh), l] <- tmp_bg
        K[(D+G+1):(dim_wh), (D+1):(G+D), l] <- t(tmp_bg)  # gamma b
        
        # gamma gamma
        idx_gamma <- (D + G + 1):dim_wh
        valori_diag <- apply(array(apply(array(apply(array(apply(lambda_xi, 2, function(x) sw*z*x), c(N,G,L)),
                                                     c(2,3), function(xx) xx*X_gamma), c(N,P_gamma,G,L)),
                                         c(3,4), function(xxx) xxx*X_gamma), c(N,P_gamma,G,L)),
                             c(2,4), sum)[, l]
        if (length(idx_gamma) == 1) {
          K[idx_gamma, idx_gamma, l] <- valori_diag
        } else {
          indici_diagonale <- cbind(idx_gamma, idx_gamma)
          K[, , l][indici_diagonale] <- valori_diag
        }
        
        wh[,l] <- -ginv(2 * K[,,l]) %*% gam[,l]
      }
      
      v <- as.matrix(wh[1:D, ])
      for(g in 1:G) b[g,] <- wh[D+g, ]
      ord.b <- order(apply(b,1,mean))
      b     <- b[ord.b,]
      gamma <- wh[(D+G+1):(dim_wh), , drop=FALSE]
      
      for(g in 1:G){
        detC <- apply(C[,,g], 2, function(x) det(matrix(x, D, D)))
        lxi[,g] <-
          (0.5 * log(detC)) +
          (0.5 * apply(rbind(C[,,g], t(mu[,,g])), 2, function(x)
            t(x[-(1:(D*D))]) %*% solve(matrix(x[1:(D*D)], nrow=D)) %*% x[-(1:(D*D))])) +
          rowSums(
            (log(sigma_xi[,,g])) - (0.5 * xi[,,g]) - (lambda_xi[,,g] * xi[,,g]^2) +
              sweep(lambda_xi[,,g], MARGIN=2, b[g,]^2, `*`) +
              (lambda_xi[,,g] * (X_gamma %*% gamma)^2) +
              (2 * sweep(lambda_xi[,,g] * (X_gamma %*% gamma), MARGIN=2, b[g,], `*`)) +
              (sweep(Y - 0.5, MARGIN=2, b[g,], `*`)) +
              ((Y-0.5) * (X_gamma %*% gamma)))
      }
    } # end D>2
    
    
    # M-step (beta) — uses full X with intercept
    lk  <- sum(sw*z*log(eta))
    it  <- 0
    lko <- lk
    XXdis <- array(0, c(G, (G-1)*P, N))
    for(i in 1:N) XXdis[,,i] <- diag(G)[,-1] %*% (diag(G-1) %x% t(X[i,]))
    
    while((lk-lko > 10^-6 & it < 100) | it == 0){
      it  <- it + 1
      lko <- lk
      sc  <- 0
      Fi  <- 0
      for(i in 1:N){
        pdis <- eta[i,]
        sc   <- sc + sw[i] * t(XXdis[,,i]) %*% (z[i,] - pdis)
        Fi   <- Fi + sw[i] * t(XXdis[,,i]) %*% (diag(pdis) - pdis %o% pdis) %*% XXdis[,,i]
      }
      dbe  <- as.vector(ginv(Fi) %*% sc)
      mdbe <- max(abs(dbe))
      if(mdbe > 0.5) dbe <- dbe/mdbe*0.5
      be0  <- c(beta)
      flag <- TRUE
      while(flag){
        beta <- be0 + dbe
        Eta  <- matrix(0, N, G)
        for(i in 1:N){
          if(ncol(X)==1) Eta[i,] <- XXdis[,,i]*beta
          else           Eta[i,] <- XXdis[,,i] %*% beta
        }
        if(max(abs(Eta)) > 100){ dbe <- dbe/2; flag <- TRUE } else flag <- FALSE
      }
      if(iter/10 == floor(iter/10)) print(beta)
      beta <- matrix(beta, P, G-1)
      exb  <- exp(X %*% beta)
      eta  <- cbind(1,exb) / (rowSums(exb)+1)
      lk   <- sum(sw*z*log(eta))
    }
    
    # E-step
    num <- eta * exp(lxi)
    if(any(is.nan(num))) browser()
    den <- apply(num, 1, sum)
    z   <- num / den
    
    ll   <- sum(sw * log(den))
    diff <- sum(abs(ll - ll_old))
    
    if(sum((ll - ll_old)) < 0) print(paste(iter, ll, sum((ll - ll_old)), "aaaaaa"))
    else                        print(paste(iter, ll, sum((ll - ll_old))))
    
  } # end while
  
  
  # Gauss-Hermite Quadrature correction
  Qd    <- npoints^D
  GaHer <- glmmML::ghq(npoints, FALSE)
  ugh   <- as.matrix(expand.grid(rep(list(GaHer$zeros), D)))
  ugh.star <- sqrt(2)*ugh
  p.gh  <- as.matrix(expand.grid(rep(list(GaHer$weights), times=D)))
  Ww    <- (2)^(D/2) * exp(apply(ugh, 1, crossprod)) * apply(p.gh, 1, prod)
  Phi   <- apply(ugh.star, 1, mvtnorm::dmvnorm)
  
  fxy <- array(0, c(N, Qd, G))
  for(g in 1:G){
    for(i in 1:N){
      Agh <- t(tcrossprod(t(v), ugh.star)) + b[g,] + (X_gamma %*% gamma)[i,]
      pgh <- 1 / (1 + exp(-Agh))
      fxy[i,,g] <- exp(tcrossprod(Y[i,], log(pgh)) + tcrossprod(1-Y[i,], log(1-pgh)))
    }
  }
  
  LLva  <- ll
  npar  <- (G-1)*P + (P_gamma * L) + (G * L) + ((L * D) - (D * (D-1)/2))
  AICva <- -2 * LLva + 2 * npar
  
  llGH1 <- apply(aperm(apply(fxy, c(1,3), function(x) x*Ww*Phi), c(2,1,3)), c(1,3), sum)
  LL    <- sum(sw * log(rowSums(t(apply(llGH1, 1, function(x) x*eta)))))
  AIC   <- -2 * LL + 2 * npar
  
  # Add intercept row (zeros) back to gamma for output consistency
  gamma_out <- rbind(rep(0, L), gamma)
  
  rownames(b) <- NULL
  colnames(b) <- colnames(b, do.NULL=FALSE, prefix="Item ")
  rownames(b) <- rownames(b, do.NULL=FALSE, prefix="Cluster ")
  rownames(v) <- rownames(v, do.NULL=FALSE, prefix="Dim ")
  colnames(v) <- colnames(v, do.NULL=FALSE, prefix="Item ")
  colnames(gamma_out) <- colnames(Y)
  rownames(gamma_out) <- colnames(X)
  colnames(beta) <- 2:G
  rownames(beta) <- colnames(X)
  
  names(LLva) <- "Log-Likelihood (variational approximation):"
  names(AICva) <- "AIC (variational approximation):"
  names(LL)   <- "Log-Likelihood (G-H Quadrature correction):"
  names(AIC)  <- "AIC (G-H Quadrature correction):"
  
  list(b=b, v=v, gamma=gamma_out, beta=beta, eta=eta, mu=mu, C=C,
       z=z, ugh.star=ugh.star, LL=LL, AIC=AIC, LLva=LLva, AICva=AICva, lxi=lxi)
}


f_mlta_wfix <- function(S, counts, G, D, tol, maxiter, pdGH)
{
  Ns <- nrow(S)
  N <- sum(counts)
  M <- ncol(S) 
  rownames(S) <- NULL
  
  # Initialize EM Algorithm
  
  eta <- rep(1 / G, G)
  p <- matrix(1 / M, G, M)
  
  lab <- sample(1:G,
                prob = eta,
                size = Ns,
                replace = TRUE)
  z <- counts * cls(lab)
  
  # Preliminary M-step
  
  eta <- colSums(z) / N
  
  ###
  
  xi <- array(20, c(Ns, M, G))
  sigma_xi <- 1 / (1 + exp(-xi))
  lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
  
  w <- matrix(rnorm(M * D), D, M)
  # b <- matrix(rnorm(M, G), G, M)
  b=matrix(rnorm(M*G,0, 0.000001), G, M) 
  ord.b=order(apply(b,1,mean))
  b <- b[ord.b,]
  wh <- matrix(0, D + G, M)
  
  C <- array(0, c(D * D, Ns, G))
  
  mu <- array(0, c(Ns, D, G))
  
  YY <- array(0, c(D * D, Ns, G))
  
  gam <- matrix(0, D + G, M)
  K <- array(diag(D + G), c(D + G, D + G, M))
  
  lxi <- matrix(0, G, Ns)
  
  # Iterative process
  
  v <- matrix(0, Ns, G)
  ll <- -Inf
  diff <- 1
  iter <- 0
  
  tol <- 10^-4	#0.05
  
  while (diff > tol & iter < maxiter)
  {
    iter <- iter + 1
    ll_old <- ll
    
    if (D == 1) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C[, , g] <-
          1 / (1 - 2 * rowSums(sweep(lambda_xi[, , g], MARGIN = 2, w ^ 2, `*`)))
        mu[, , g] <-
          C[, , g] * rowSums(sweep(
            S - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`),
            MARGIN = 2, w, `*`))
        
        YY[, , g] <- matrix(C[, , g] + mu[, , g] ^ 2, ncol = 1)
        
        xi[, , g] <- YY[, , g] %*% w^2 + mu[, , g] %*% (2 * b[g, ] * w) + 
          matrix(b[g, ]^2, nrow = Ns, ncol = M, byrow = TRUE)
      }
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      gam[1, ] <- colSums(crossprod(z * mu[, 1, ], S - 0.5))
      gam[2:(1 + G), ] <- crossprod(z, S - 0.5)
      
      K[1, 1, ] <- apply(aperm(array(z * YY[1, , ], c(Ns, G, M)), c(1, 3, 2)) * lambda_xi, 2, sum)
      
      K[2:(G + 1), 1, ] <- t(apply(aperm(array(
        z * mu[, 1, ], c(Ns, G, M)
      ), c(1, 3, 2)) * lambda_xi, c(2, 3), sum))
      K[1, 2:(G + 1), ] <- K[2:(G + 1), 1, ]
      
      for (m in 1:M)
      {
        for (g in 1:G)
          K[1 + g, 1 + g, m] <- sum(z[, g] * lambda_xi[, m, g])
        
        wh[, m] <- -solve(2 * K[, , m]) %*% gam[, m]
      }
      
      w <- wh[1:D, ]
      
      w <- as.matrix(t(w))
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
        
        lxi[g, ] <-
          0.5 * log(C[1, , g]) + mu[, 1, g] ^ 2 / (2 * C[1, , g]) + rowSums(
            log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] * xi[, , g] ^ 2 + 
              sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + 
              sweep(S - 0.5, MARGIN = 2, b[g, ], `*`))
      } # end for (g in 1:G)
    }
    
    if (D == 2) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C[, , g] <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(w), t(w))))
        
        mu[, , g] <-
          t(apply(rbind(C[, , g], 
                        tcrossprod(w, S - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`))), 
                  2, function(x) matrix(x[1:4], nrow = D) %*% x[-(1:4)]))
        
        # # STEP 2: Optimising the Variational Parameters (xi)
        
        YY[, , g] <- C[, , g] + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY[, , g], 2, function(x)
              rowSums(crossprod(w, matrix(x, ncol = D)) * t(w))) + 
              tcrossprod(2 * b[g, ] * t(w), mu[, , g]) + matrix(
                b[g, ] ^ 2,
                nrow = M,
                ncol = Ns,
                byrow = FALSE
              )
          )
      }
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      # wh<-rbind(w,b)
      gam[3:(2 + G), ] <- crossprod(z, S - 0.5)
      
      aa <- aperm(array(z, c(Ns, G, D)), c(1, 3, 2)) * mu
      bb <- aperm(array(z, c(Ns, G, D * D)), c(3, 1, 2)) * YY
      
      K[3:(2 + G), 3:(2 + G), ] <-
        apply(apply(aperm(array(z, c(
          Ns, G, M
        )), c(1, 3, 2)) * lambda_xi, c(2, 3), sum), 1, diag)
      
      kk <- 0
      for (g in 1:G) {
        K[2 + g, 1:2, ] <- crossprod(aa[, , g], lambda_xi[, , g])
        K[1:2, 2 + g, ] <- K[2 + g, (1:2), ]
        kk <- kk + bb[, , g] %*% lambda_xi[, , g]
      }
      
      K[1:2, 1:2, ] <- kk
      
      for (m in 1:M)
      {
        gam[1:D, m] <- apply(aa * (S[, m] - 0.5), 2, sum)
        wh[, m] <- -solve(2 * K[, , m]) %*% gam[, m]
      }
      
      w <- wh[1:D, ]
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
        
        # Approximation of log(p(x|z))
        
        detC <- C[1, , g] * C[4, , g] - C[3, , g] * C[2, , g]
        
        lxi[g, ] <-
          0.5 * log(detC) + 0.5 * apply(rbind(C[4, , g] / detC, -C[2, , g] / detC, -C[3, , g] /
                                                detC, C[1, , g] / detC, t(mu[, , g])), 2, function(x)
                                                  t((x[-(1:4)])) %*% matrix(x[1:4], nrow = D) %*% (x[-(1:4)])) + 
          rowSums(
            log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] * xi[, , g]^2 + 
              sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + 
              sweep(S - 0.5, MARGIN = 2, b[g, ], `*`))
      } # end for (g in 1:G)
      
    }# end D=2
    
    if (D > 2) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C[, , g] <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(w), t(w))))
        
        mu[, , g] <-
          t(apply(rbind(C[, , g], w %*% t(
            S - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`)
          )), 2, function(x)
            matrix(x[1:(D * D)], nrow = D) %*% x[-(1:(D * D))]))
        
        # # STEP 2: Optimising the Variational Parameters (xi)
        
        YY[, , g] <- C[, , g] + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY[, , g], 2, function(x)
              rowSums((t(w) %*% matrix(x, ncol = D)) * t(w))) + 
              (2 * b[g, ] * t(w)) %*% t(mu[, , g]) + 
              matrix(
                b[g, ] ^ 2,
                nrow = M,
                ncol = Ns,
                byrow = FALSE
              )
          )
        
      }
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      gam[(D + 1):(D + G), ] <- t(z) %*% (S - 0.5)
      
      aa <- aperm(array(z, c(Ns, G, D)), c(1, 3, 2)) * mu
      bb <- aperm(array(z, c(Ns, G, D * D)), c(3, 1, 2)) * YY
      
      K[(D + 1):(D + G), (D + 1):(D + G), ] <-
        apply(apply(aperm(array(z, c(
          Ns, G, M
        )), c(1, 3, 2)) * lambda_xi, c(2, 3), sum), 1, diag)
      
      kk <- 0
      for (g in 1:G) {
        K[D + g, 1:D, ] <- crossprod(aa[, , g], lambda_xi[, , g])
        K[1:D, D + g, ] <- K[D + g, 1:D, ]
        kk <- kk + bb[, , g] %*% lambda_xi[, , g]
      }
      
      K[1:D, 1:D, ] <- kk
      
      for (m in 1:M)	{
        gam[1:D, m] <- apply(aa * (S[, m] - 0.5), 2, sum)
        wh[, m] <- -solve(2 * K[, , m]) %*% gam[, m]
      }
      
      w <- wh[1:D, ]
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
        
        # Approximation of log(p(x|z))
        
        detC <- apply(C[, , g], 2, function(x)
          det(matrix(x, D, D)))
        
        lxi[g, ] <-
          0.5 * log(detC) + 0.5 * apply(rbind(C[, , g], t(mu[, , g])), 2, function(x)
            t((x[-(1:(D * D))])) %*% solve(matrix(x[1:(D * D)], nrow = D)) %*% (x[-(
              1:(D * D))])) + rowSums(
                log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] *
                  xi[, , g]^2 + sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + 
                  sweep(S - 0.5, MARGIN = 2, b[g, ], `*`)
              )
      } # end for (g in 1:G)
      
    } # end if D > 2
    
    # E-step
    
    v <- t(eta * exp(lxi))
    v[is.nan(v)] <- 0
    vsum <- apply(v, 1, sum)
    z <- counts * v / vsum
    ll <- sum(counts * log(vsum))
    
    # M-step
    
    eta <- apply(z, 2, sum) / N
    
    # Stopping Criteria
    
    diff <- sum(abs(ll - ll_old))
    
  } # end while(diff>tol)
  
  # Correction to the log-likelihood 
  
  # Gauss-Hermite Quadrature
  
  npoints <- round(pdGH ^ (1 / D))
  ny <- npoints ^ D
  GaHer <- glmmML::ghq(npoints, FALSE)
  Ygh <- expand.grid(rep(list(GaHer$zeros), D))
  Ygh <- as.matrix(Ygh)
  Wgh <-
    apply(as.matrix(expand.grid(rep(
      list(GaHer$weights), D
    ))), 1, prod) * apply(exp(Ygh ^ 2), 1, prod)
  
  Hy <- apply(Ygh, 1, mvtnorm::dmvnorm)
  Beta <- Hy * Wgh / sum(Hy * Wgh)
  
  fxy <- array(0, c(Ns, ny, G))
  
  for (g in 1:G)
  {
    Agh <- t(tcrossprod(t(w), Ygh) + b[g, ])
    pgh <- 1 / (1 + exp(-Agh))
    fxy[, , g] <-
      exp(tcrossprod(S, log(pgh)) + tcrossprod(1 - S, log(1 - pgh)))
  }
  
  eta <- as.vector(eta)
  
  LLva <- ll
  BICva <- -2 * LLva + {
    G * M + M * D - D * {
      D - 1
    } / 2 + G - 1
  } * log(N)
  
  llGH1 <- apply(rep(Beta, each = Ns) * fxy, c(1, 3), sum)
  LL <- sum(counts * log(colSums(eta * t(llGH1))))
  
  BIC <- -2 * LL + {
    G * M + M * D - D * {
      D - 1
    } / 2 + G - 1
  } * log(N)
  
  expected <- colSums(eta * t(llGH1)) * N
  
  rownames(b) <- NULL
  colnames(b) <- colnames(b, do.NULL = FALSE, prefix = "Item ")
  rownames(b) <- rownames(b, do.NULL = FALSE, prefix = "Group ")
  
  rownames(w) <- rownames(w, do.NULL = FALSE, prefix = "Dim ")
  colnames(w) <- colnames(w, do.NULL = FALSE, prefix = "Item ")
  
  eta <- matrix(eta, byrow = T, G)
  rownames(eta) <- rownames(eta, do.NULL = FALSE, prefix = "Group ")
  colnames(eta) <- "eta"
  names(LLva) <- c("Log-Likelihood (variational approximation):")
  names(BICva) <- c("BIC (variational approximation):")
  names(LL) <- c("Log-Likelihood (G-H Quadrature correction):")
  names(BIC) <- c("BIC (G-H Quadrature correction):")
  
  out1 <-
    list(
      b = b,
      w = w,
      eta = eta,
      LL = LL,
      BIC = BIC,
      LLva = LLva,
      BICva = BICva,
      expected = expected,
      mu = mu,
      C = C,
      z = z
    )
  
  return(out1)
  
}


f_mlta_wfix.conc <- function(X, DM, G, D, tol, maxiter, pdGH, beta0)
{
  
  N <- nrow(X)
  M <- ncol(X) 
  J <- ncol(DM)
  
  # Initialize EM Algorithm
  
  if(!is.null(beta0)) beta = beta0 else{
    beta <- rep(0,J*(G-1))
    beta <- matrix(beta,ncol=G-1)
  }
  
  exb <- exp(DM %*% beta)
  eta <- cbind(1,exb)/(rowSums(exb)+1)  # priors
  
  
  z <- matrix(NA, nrow=N, ncol=G)
  for(i in 1:N)   z[i,] <- t(rmultinom(1, size = 1, prob = eta[i,]))
  
  # Set up
  
  p <- (t(z) %*% X) / colSums(z)
  
  
  ###
  
  xi <- array(20, c(N, M, G))
  sigma_xi <- 1 / (1 + exp(-xi))
  lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
  
  w <- matrix(rnorm(M * D), D, M)
  # b <- matrix(rnorm(M, G), G, M)
  b=matrix(rnorm(M*G,0, 0.000001), G, M) 
  ord.b=order(apply(b,1,mean))
  b <- b[ord.b,]
  
  wh <- matrix(0, D + G, M)
  
  C <- array(0, c(D * D, N, G))
  
  mu <- array(0, c(N, D, G))
  
  YY <- array(0, c(D * D, N, G))
  
  gam <- matrix(0, D + G, M)
  K <- array(diag(D + G), c(D + G, D + G, M))
  
  lxi <- matrix(0, G, N)
  
  # Iterative process
  
  v <- matrix(0, N, G)
  W <- list()
  beta_change <- rep(NA, (J)*(G-1))
  beta_change <- matrix(beta_change, ncol=(G-1), nrow=J)
  se.beta <- rep(NA, (J)*(G-1))
  se.beta <- matrix(se.beta, ncol=(G-1), nrow=J)
  ll <- -Inf
  diff <- 1
  iter <- 0
  cond <- TRUE
  
  tol <- 10^-4
  print(c(beta))
  while (diff > tol & iter < maxiter)
  {
    iter <- iter + 1
    beta.old=beta
    ll_old <- ll
    
    if (D == 1) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C[, , g] <-
          1 / (1 - 2 * rowSums(sweep(lambda_xi[, , g], MARGIN = 2, w ^ 2, `*`)))
        mu[, , g] <-
          C[, , g] * rowSums(sweep(
            X - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`),
            MARGIN = 2, w, `*`))
        
        YY[, , g] <- matrix(C[, , g] + mu[, , g] ^ 2, ncol = 1)
        
        xi[, , g] <- YY[, , g] %*% w^2 + mu[, , g] %*% (2 * b[g, ] * w) + 
          matrix(b[g, ]^2, nrow = N, ncol = M, byrow = TRUE)
      }
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      gam[1, ] <- colSums(crossprod(z * mu[, 1, ], X - 0.5))
      gam[2:(1 + G), ] <- crossprod(z, X - 0.5)
      
      K[1, 1, ] <- apply(aperm(array(z * YY[1, , ], c(N, G, M)), c(1, 3, 2)) * lambda_xi, 2, sum)
      
      K[2:(G + 1), 1, ] <- t(apply(aperm(array(
        z * mu[, 1, ], c(N, G, M)
      ), c(1, 3, 2)) * lambda_xi, c(2, 3), sum))
      K[1, 2:(G + 1), ] <- K[2:(G + 1), 1, ]
      
      for (m in 1:M)
      {
        for (g in 1:G)
          K[1 + g, 1 + g, m] <- sum(z[, g] * lambda_xi[, m, g])
        
        wh[, m] <- -solve(2 * K[, , m]) %*% gam[, m]
      }
      
      w <- wh[1:D, ]
      
      w <- as.matrix(t(w))
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
        
        lxi[g, ] <-
          0.5 * log(C[1, , g]) + mu[, 1, g] ^ 2 / (2 * C[1, , g]) + rowSums(
            log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] * xi[, , g] ^ 2 + 
              sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + 
              sweep(X - 0.5, MARGIN = 2, b[g, ], `*`))
      } # end for (g in 1:G)
    }
    
    if (D == 2) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C[, , g] <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(w), t(w))))
        
        mu[, , g] <-
          t(apply(rbind(C[, , g], 
                        tcrossprod(w, X - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`))), 
                  2, function(x) matrix(x[1:4], nrow = D) %*% x[-(1:4)]))
        
        # # STEP 2: Optimising the Variational Parameters (xi)
        
        YY[, , g] <- C[, , g] + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY[, , g], 2, function(x)
              rowSums(crossprod(w, matrix(x, ncol = D)) * t(w))) + 
              tcrossprod(2 * b[g, ] * t(w), mu[, , g]) + matrix(
                b[g, ] ^ 2,
                nrow = M,
                ncol = N,
                byrow = FALSE
              )
          )
      }
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      gam[3:(2 + G), ] <- crossprod(z, X - 0.5)
      
      aa <- aperm(array(z, c(N, G, D)), c(1, 3, 2)) * mu
      bb <- aperm(array(z, c(N, G, D * D)), c(3, 1, 2)) * YY
      
      K[3:(2 + G), 3:(2 + G), ] <-
        apply(apply(aperm(array(z, c(
          N, G, M
        )), c(1, 3, 2)) * lambda_xi, c(2, 3), sum), 1, diag)
      
      kk <- 0
      for (g in 1:G) {
        K[2 + g, 1:2, ] <- crossprod(aa[, , g], lambda_xi[, , g])
        K[1:2, 2 + g, ] <- K[2 + g, (1:2), ]
        kk <- kk + bb[, , g] %*% lambda_xi[, , g]
      }
      
      K[1:2, 1:2, ] <- kk
      
      for (m in 1:M)
      {
        gam[1:D, m] <- apply(aa * (X[, m] - 0.5), 2, sum)
        wh[, m] <- -solve(2 * K[, , m]) %*% gam[, m]
      }
      
      w <- wh[1:D, ]
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
        
        # Approximation of log(p(x|z))
        
        detC <- C[1, , g] * C[4, , g] - C[3, , g] * C[2, , g]
        
        lxi[g, ] <-
          0.5 * log(detC) + 0.5 * apply(rbind(C[4, , g] / detC, -C[2, , g] / detC, -C[3, , g] /
                                                detC, C[1, , g] / detC, t(mu[, , g])), 2, function(x)
                                                  t((x[-(1:4)])) %*% matrix(x[1:4], nrow = D) %*% (x[-(1:4)])) + 
          rowSums(
            log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] * xi[, , g]^2 + 
              sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + 
              sweep(X - 0.5, MARGIN = 2, b[g, ], `*`))
      } # end for (g in 1:G)
      
    }# end D=2
    
    if (D > 2) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C[, , g] <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(w), t(w))))
        
        mu[, , g] <-
          t(apply(rbind(C[, , g], w %*% t(
            X - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`)
          )), 2, function(x)
            matrix(x[1:(D * D)], nrow = D) %*% x[-(1:(D * D))]))
        
        # # STEP 2: Optimising the Variational Parameters (xi)
        
        YY[, , g] <- C[, , g] + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY[, , g], 2, function(x)
              rowSums((t(w) %*% matrix(x, ncol = D)) * t(w))) + 
              (2 * b[g, ] * t(w)) %*% t(mu[, , g]) + 
              matrix(
                b[g, ] ^ 2,
                nrow = M,
                ncol = N,
                byrow = FALSE
              )
          )
        
      }
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      gam[(D + 1):(D + G), ] <- t(z) %*% (X - 0.5)
      
      aa <- aperm(array(z, c(N, G, D)), c(1, 3, 2)) * mu
      bb <- aperm(array(z, c(N, G, D * D)), c(3, 1, 2)) * YY
      
      K[(D + 1):(D + G), (D + 1):(D + G), ] <-
        apply(apply(aperm(array(z, c(
          N, G, M
        )), c(1, 3, 2)) * lambda_xi, c(2, 3), sum), 1, diag)
      
      kk <- 0
      for (g in 1:G) {
        K[D + g, 1:D, ] <- crossprod(aa[, , g], lambda_xi[, , g])
        K[1:D, D + g, ] <- K[D + g, 1:D, ]
        kk <- kk + bb[, , g] %*% lambda_xi[, , g]
      }
      
      K[1:D, 1:D, ] <- kk
      
      for (m in 1:M)	{
        gam[1:D, m] <- apply(aa * (X[, m] - 0.5), 2, sum)
        wh[, m] <- -solve(2 * K[, , m]) %*% gam[, m]
      }
      
      w <- wh[1:D, ]
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
        
        # Approximation of log(p(x|z))
        
        detC <- apply(C[, , g], 2, function(x)
          det(matrix(x, D, D)))
        
        lxi[g, ] <-
          0.5 * log(detC) + 0.5 * apply(rbind(C[, , g], t(mu[, , g])), 2, function(x)
            t((x[-(1:(D * D))])) %*% solve(matrix(x[1:(D * D)], nrow = D)) %*% (x[-(
              1:(D * D))])) + rowSums(
                log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] *
                  xi[, , g]^2 + sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + 
                  sweep(X - 0.5, MARGIN = 2, b[g, ], `*`)
              )
      } # end for (g in 1:G)
      
    } # end if D > 2
    
    # M-step 
    
    lk = sum(z*log(eta))
    it = 0; lko = lk
    XXdis = array(0,c(G,(G-1)*ncol(DM),N))
    for(i in 1:N){
      XXdis[,,i] = diag(G)[,-1]%*%(diag(G-1)%x%t(DM[i,]))
    }
    while((lk-lko>10^-6 & it<100) | it==0){
      it = it+1; lko = lk 
      sc = 0; Fi = 0
      for(i in 1:N){
        pdis = eta[i,]
        sc = sc+t(XXdis[,,i])%*%(z[i,]-pdis)
        Fi = Fi+t(XXdis[,,i])%*%(diag(pdis)-pdis%o%pdis)%*%XXdis[,,i]
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
          
          if(ncol(DM)==1) Eta[i,] = XXdis[,,i]*beta
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
      
      beta = matrix(beta, J, G-1)    
      exb <- exp(DM %*% beta) # updfe priors
      eta <- cbind(1,exb)/(rowSums(exb)+1)
      
      lk = sum(z*log(eta))
    }
    
    
    # E-step
    v <- eta * t(exp(lxi))
    if(any(is.nan(v))) browser()
    #v[is.nan(v)] <- 0
    vsum <- apply(v, 1, sum)
    z <- v / vsum
    ll <- sum(rep(1, N) * log(vsum))
    
    # Stopping Criteria
    
    diff <- sum(abs(ll - ll_old))
    
  } # end while(diff>tol)
  
  # Correction to the log-likelihood 
  
  # Gauss-Hermite Quadrature
  
  npoints <- round(pdGH ^ (1 / D))
  ny <- npoints ^ D
  GaHer <- glmmML::ghq(npoints, FALSE)
  Ygh <- expand.grid(rep(list(GaHer$zeros), D))
  Ygh <- as.matrix(Ygh)
  Wgh <-
    apply(as.matrix(expand.grid(rep(
      list(GaHer$weights), D
    ))), 1, prod) * apply(exp(Ygh ^ 2), 1, prod)
  
  Hy <- apply(Ygh, 1, mvtnorm::dmvnorm)
  Beta <- Hy * Wgh / sum(Hy * Wgh)
  
  fxy <- array(0, c(N, ny, G))
  
  for (g in 1:G)
  {
    Agh <- t(tcrossprod(t(w), Ygh) + b[g, ])
    pgh <- 1 / (1 + exp(-Agh))
    fxy[, , g] <-
      exp(tcrossprod(X, log(pgh)) + tcrossprod(1 - X, log(1 - pgh)))
  }
  
  
  LLva <- ll
  BICva <- -2 * LLva + {
    G * M + M * D - D * {
      D - 1
    } / 2 + J*(G - 1)
  } * log(N)
  
  llGH1 <- apply(rep(Beta, each = N) * fxy, c(1, 3), sum)
  LL <- sum(log(rowSums(eta * (llGH1))))
  
  BIC <- -2 * LL + {
    G * M + M * D - D * {
      D - 1
    } / 2 + J*(G - 1)
  } * log(N)
  
  expected <- colSums(eta * (llGH1)) * N
  
  q <- qnorm(0.975) # CI
  u.beta <- beta + q*se.beta
  l.beta <- beta - q*se.beta
  
  rownames(b) <- NULL
  colnames(b) <- colnames(b, do.NULL = FALSE, prefix = "Item ")
  rownames(b) <- rownames(b, do.NULL = FALSE, prefix = "Group ")
  
  rownames(w) <- rownames(w, do.NULL = FALSE, prefix = "Dim ")
  colnames(w) <- colnames(w, do.NULL = FALSE, prefix = "Item ")
  
  colnames(beta) <- 2:G
  rownames(beta) <- c(colnames(DM))
  
  names(LLva) <- c("Log-Likelihood (variational approximation):")
  names(BICva) <- c("BIC (variational approximation):")
  names(LL) <- c("Log-Likelihood (G-H Quadrature correction):")
  names(BIC) <- c("BIC (G-H Quadrature correction):")
  
  out1 <-
    list(
      b = b,
      w = w,
      eta = eta,
      LL = LL,
      BIC = BIC,
      LLva = LLva,
      BICva = BICva,
      expected = expected,
      mu = mu,
      C = C,
      z = z,
      beta = beta, 
      se.beta = se.beta, 
      rrr=exp(beta), 
      u.beta=u.beta, 
      l.beta=l.beta
    )
  
  out1
  
}
