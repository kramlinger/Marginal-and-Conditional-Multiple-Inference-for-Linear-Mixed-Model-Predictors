
# in order to create table 4 in the supplement, run the parametric bootstrap 
# "bt" for all tests. 


bt <- function(TEST) {
  
  library(dplyr)
  library(purrr)
  
  ni <- tally(prisons)$n
  prisons.original <- prisons %>% left_join(tally(prisons), by = "facility_state")
  
  lambda.bt <- matrix(NA, 1000, 6)
  
  set.seed(12345)
  for (oo in 1:1000) { # number of bootstrap replications 
    
    # sample according to original n_i
    asdf <- mapply(function(n, data) data[sample(1:n, n, TRUE),], 
                   ni, group_split(prisons.original), SIMPLIFY = FALSE)
    
    # bring data frame of bootstrap replication in original format
    prisons <- do.call("rbind", asdf)  %>% group_by(facility_state)

    # calculate states data frame
    states <- summarise(prisons, 
                        county_log_mortality = mean(county_log_mortality),
                        prison_log_mortality = mean(prison_log_mortality),
                        weights = sum(weights))
    
    # repeat the calculation of covariance matrix as in study.R
    
    ### FIT MODEL 
    library("lme4")
    fit <- lmer(prison_log_mortality ~ county_log_mortality + (1 | facility_state), 
                data = prisons) # log transformed weights give better fit

    ### ESTIMATE COVARIANCE MATRICES
    ni <- tally(prisons)$n
    m <- length(ni)
    n <- sum(ni)
    aij <- prisons$weights
    ai <- states$weights
    
    delta <- as.data.frame(VarCorr(fit))$vcov
    beta <- fixef(fit)
    
    Z <- getME(fit, "Z")
    X <- getME(fit, "X")
    y <- getME(fit, "y")
    l <- cbind(1, states$county_log_mortality)
    V <- delta[1] * tcrossprod(Z) + diag(delta[2] / aij, n)
    Vi <- solve(V)
    
    gamma <- delta[1] / (delta[1] + delta[2] / ai)
    b <- delta[1] * t(Z) %*% Vi
    d <- l - b  %*% X
    
    mu <- l %*% beta + b %*% (y - X %*% beta) # = l %*% beta + getME(fit, "b")
    
    ### SIGMA 
    g1 <- delta[1] * diag(1, m) - delta[1]^2 * t(Z) %*% Vi %*% Z
    
    a <- ai * delta[1] + delta[2]
    Ive <- sum(ai / a^2) / 2
    Ivv <- sum(ai^2 / a^2) / 2
    Iee <- sum( (ni - 1) / delta[2]^2 + 1 / a^2) / 2
    I <- matrix(c(Ivv, Ive, Ive, Iee), 2, 2)
    Vbar <- solve(I)
    
    
    Vv <- tcrossprod(Z)
    Ve <- diag(1 / aij)
    bv <- t(Z) %*% Vi - delta[1] * t(Z) %*% Vi %*% Vv %*% Vi
    be <- -delta[1] * t(Z) %*% Vi %*% Ve %*% Vi
    # Db <- list(bv, be)
    
    b <- split(colSums(b), rep(1:m, ni))
    bv <- split(colSums(bv), rep(1:m, ni))
    be <- split(colSums(be), rep(1:m, ni))
    Db <- mapply(cbind, bv, be)
    
    Vj <- list(); NI <- c(0, cumsum(ni))
    for (i in 1:m) Vj[[i]] <- V[(NI[i] + 1):NI[i + 1], (NI[i] + 1):NI[i + 1]] 
    
    
    g3 <- c()
    for (i in 1:m) g3[i] <- sum(diag(t(Db[[i]]) %*% Vj[[i]] %*% Db[[i]] %*% Vbar))
    
    Sigma <- g1 + diag(2 * g3) + d %*% vcov(fit) %*% t(d)
    
    
    
    ### cSigma
    l1 <- c()
    R <- diag(delta[2] / aij)
    Rj <- mapply(diag, split(delta[2] / aij, rep(1:m, ni)), ni)
    for (i in 1:m) l1[i] <- b[[i]] %*% Rj[[i]] %*% b[[i]]
    
    K <- diag(delta[2] /aij) %*% Vi %*% X %*% vcov(fit)
    Kj <- list() 
    for (i in 1:m) Kj[[i]] <- K[(NI[i] + 1):NI[i + 1], ] 
    KRK <- mapply(function(a,b){
      if(length(b) != 1) crossprod(a, solve(b)) %*% a #t(a) %*% diag(1/b) %*% a
      else{ as.matrix(a); a %*% solve(b) %*% t(a)}
    }, Kj, Rj, SIMPLIFY = F)
    KRK <- as.matrix(Reduce("+", KRK))
    
    L2 <- matrix(NA, m, m)
    for (i in 1:m) {
      for (j in 1:m) {
        L2[i,j] <- as.numeric(d[i,] %*% t(b[[j]] %*% Kj[[j]])) + 
          as.numeric(d[j,] %*% t(b[[i]] %*% Kj[[i]])) + 
          d[j,] %*% KRK %*% d[i,]
      }
    }
    
    l4 <- c() 
    for (i in 1:m) l4[i] <- sum(diag(t(Db[[i]]) %*% Rj[[i]] %*% Db[[i]] %*% Vbar)) / ni[i]
    
    bvv <- -t(Z) %*% Vi %*% Vv %*% Vi - t(Z) %*% Vi %*% Vv %*% Vi + 
      2 * delta[1] * t(Z) %*% Vi %*% Vv %*% Vi %*% Vv %*% Vi
    bee <- 2 * delta[1] * t(Z) %*% Vi %*% Ve %*% Vi %*% Ve %*% Vi
    bve <- -t(Z) %*% Vi %*% Ve %*% Vi + 
      delta[1] * t(Z) %*% Vi %*% Ve %*% Vi %*% Vv %*% Vi + 
      delta[1] * t(Z) %*% Vi %*% Vv %*% Vi %*% Ve %*% Vi
    bvv <- split(colSums(bvv), rep(1:m, ni))
    bee <- split(colSums(bee), rep(1:m, ni))
    bve <- split(colSums(bve), rep(1:m, ni))
    D2b <- mapply(function(a,b,c){
      D2b <- list(list(), list())
      D2b[[1]][[1]] <- a
      D2b[[1]][[2]] <- b
      D2b[[2]][[1]] <- b
      D2b[[2]][[2]] <- c
      D2b
    }, bvv, bve, bee,SIMPLIFY = F)
    
    D2L <- list() 
    for (i in 1:m) {
      DLvv <- bv[[i]] %*% Rj[[i]] %*% bv[[i]] + b[[i]] %*% Rj[[i]] %*% bvv[[i]]
      DLve <- be[[i]] %*% Rj[[i]] %*% bv[[i]] + b[[i]] %*% Rj[[i]] %*% bve[[i]] + 
        b[[i]] %*% Rj[[i]] %*% bv[[i]] / delta[2]
      DLee <-  2 * b[[i]] %*% Rj[[i]] %*% be[[i]] / delta[2] + 
        b[[i]] %*% Rj[[i]] %*% bee[[i]] + be[[i]] %*% Rj[[i]] %*% be[[i]]
      D2L[[i]] <- matrix(c(DLvv, DLve, DLve, DLee), 2, 2)
    }
    
    l5 <- c() 
    for (i in 1:m) {
      l5[i] <- sum(diag(D2L[[i]] %*% Vbar))
    }
    
    
    l3.1 <- 0
    dV <- list(Vv, Ve)
    w <- sapply(b, sum) %*% t(Z); w <- bdiag(split(w, rep(1:m, ni)))
    uu <- list(); for(i in 1:m) uu[[i]] <- w[,i]; w <- uu; rm(uu)
    
    wv <- sapply(bv, sum) %*% t(Z); wv <- bdiag(split(wv, rep(1:m, ni)))
    uu <- list(); for(i in 1:m) uu[[i]] <- wv[,i]; wv <- uu; rm(uu)
    we <- sapply(be, sum) %*% t(Z); we <- bdiag(split(we, rep(1:m, ni)))
    uu <- list(); for(i in 1:m) uu[[i]] <- we[,i]; we <- uu; rm(uu)
    Dw <- mapply(rbind, wv, we, SIMPLIFY = F)
    
    wvv <- sapply(bvv, sum) %*% t(Z); wvv <- bdiag(split(wvv, rep(1:m, ni)))
    uu <- list(); for(i in 1:m) uu[[i]] <- wvv[,i]; wvv <- uu; rm(uu)
    wve <- sapply(bve, sum) %*% t(Z); wve <- bdiag(split(wve, rep(1:m, ni)))
    uu <- list(); for(i in 1:m) uu[[i]] <- wve[,i]; wve <- uu; rm(uu)
    wee <- sapply(bee, sum) %*% t(Z); wee <- bdiag(split(wee, rep(1:m, ni)))
    uu <- list(); for(i in 1:m) uu[[i]] <- wee[,i]; wee <- uu; rm(uu)
    D2w <- mapply(function(a,b,c){
      D2b <- list(list(), list())
      D2b[[1]][[1]] <- a
      D2b[[1]][[2]] <- b
      D2b[[2]][[1]] <- b
      D2b[[2]][[2]] <- c
      D2b
    }, wvv, wve, wee,SIMPLIFY = F)
    
    
    
    l3.1 <- rep(0, m)
    for (i in 1:m) {
      for (e in 1:2) {
        l3.1[i] <- l3.1[i] + Vbar[e,] %*% Dw[[i]] %*% R %*% Vi %*% dV[[e]] %*% Vi %*% R %*% w[[i]] / ni[i]^2
      }
    }
    
    l3.2 <- rep(0, m)
    for (i in 1:m) {
      for (e in 1:2) {
        for (d in 1:2) {
          for (f in 1:2) {
            for (g in 1:2) {
              l3.2[i] <- l3.2[i] + 4 * Vbar[e,f] * Vbar[f,g] * I[e,d] * D2w[[i]][[e]][[d]] %*% R %*% w[[i]] / ni[i]^3
            }
          }
        }
      }
    }
    
    DvI <- -matrix(c(sum(ai^3 / a^3), sum(ai^2 / a^3), sum(ai^2 / a^3), 
                     sum(ai / a^3)), 2, 2)
    DeI <- -matrix(c(sum(ai^2 / a^3), sum(ai / a^3), sum(ai / a^3), 
                     sum((ni - 1) / delta[2]^3 - 1 / a^3)), 2, 2)
    DI <- list(DvI, DeI)
    
    l3.3 <- rep(0, m)
    for (i in 1:m) {
      for (e in 1:2) {
        for (d in 1:2) {
          for (f in 1:2) {
            l3.3[i] <- l3.3[i] - 2 * Vbar[e,f] * I[e,d] * Vbar[d,] %*% DI[[e]] %*% Vbar %*% Dw[[i]] %*% R %*% w[[i]] / ni[i]^3
          }
        }
      }
    }
    
    l3.4 <- rep(0, m)
    for (i in 1:m) {
      for (e in 1:2) {
        for (d in 1:2) {
          for (g in 1:2) {
            l3.4[i] <- l3.4[i] + 2 * DI[[g]][e,d] * Vbar[e,d] * Vbar[e,] %*% Dw[[i]] %*% R %*% w[[i]] / ni[i]^3
          }
        }
      }
    }
    
    l3 <- c(l3.1 + l3.2 + l3.3 + l3.4)
    
    SigmaC <- diag(l1 + l4 + 2 * l3 - l5) + L2
    
    
    ### SPECIFY NONCENTRALITY ESTIMATOR 
    Xj <- mapply(function(a) t(matrix(a, 2, byrow = T)), split(X, rep(1:m, ni)))
    d <- l - t(mapply("%*%", b, Xj))
    A <- (sapply(b, sum) - 1) / ni * t(Z) + d %*% vcov(fit) %*% t(X) %*% Vi
    lambda <- function(L) {
      A <- L %*% A
      as.numeric(t(A %*% y) %*% solve(L %*% SigmaC %*% t(L)) %*% A %*% y +
                   -sum(diag(solve(R) %*% t(A) %*% solve(L %*% SigmaC %*% t(L)) %*% A)) + 
                   -t(A %*% X %*% beta) %*% solve(L %*% SigmaC %*% t(L)) %*% A %*% X %*% beta)
    }
    
    
    lamfun <- function(test) {
      test <- states$facility_state[states$facility_state %in% test]
      
      
      mstar <- length(test) - 1
      Lstar <- cbind(diag(mstar), 0) - matrix(1/(mstar + 1), mstar, mstar + 1)
      L <- matrix(0, mstar, 45)
      L[,(1:45)[states$facility_state %in% test]] <- Lstar
      lam <- lambda(L)
    }
    
    
    lambda.bt[oo,] <- sapply(TEST, lamfun)
  }
  return(lambda.bt)
}
