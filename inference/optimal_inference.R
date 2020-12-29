# this file contains the necessary routines to execute data carving after model selection
# either for single variable testing or group testing
# the sampler is defined elsewhere


carve.lasso <- function(X, y, ind, beta, tol.beta, lambda, sigma = NULL, family = "gaussian",
                        intercept = TRUE, which.check = NULL, selected = TRUE, ndraw = 8000, burnin = 2000,
                        sig.level = 0.05, FWER = TRUE, aggregation = 0.05,  time.constant = 1e-6, verbose = FALSE) {
  # to be applied after Lasso Selection
  # X: full X matrix
  # y: full y vector
  # ind: indices used for selection, i.e. Lasso was applied to X[ind, ] and y[ind]
  # beta: coefficients obtained from Lasso selection
  # tol.beta tolerance, to assume variable to be active
  # lambda: penalty parameter used for Lasso 
  # (with objective 1/2||X*\beta-y||^2+\lambda*||\beta||_1, i.e different normalization than in glmnet)
  # sigma: standard deviation, assumed to be known for Gaussian
  # family: Gaussian or binomial
  # intercept: was the model fit using an intercept
  # which.check: which variables shall be checked (might skip "unimportant" ones)
  # selected: whether to use the selected viewpoint for aggregation
  # otherwise no sampling is of need
  # ndraw: number of points to draw, when using hit-and-run sampler
  # burnin: number of initial points to burn
  # ndraw and burning might be increased for high degrees of freedom or due to required significance
  # the following three values are used to defined the minimally required sample size
  # sig.level: level for the hypothesis test
  # FWER: whether a FWER correction will be applied
  # aggregation: aggregation parameter \gamma.min
  # time.constant: the Hamiltonian sampler is assumed to be stuck, if it does not
  # finish after (ndraw + burnin) * number of constraints * dimensionality * time.constant seconds
  # verbose: whether to print key steps
  
  knownfamilies <- c("gaussian", "binomial")
  if (!(family %in% knownfamilies)) {
    stop(paste("Unknown family, family should be one of", paste(knownfamilies, collapse = ", ")))
  }
  n <- dim(X)[1]
  p <- dim(X)[2]
  splitn <- length(ind)
  if (family == "gaussian") {
    if (is.null(sigma)) stop("Sigma needs to be provided for Gaussian family")
    if (length(beta) == p + 1) {
      if (!intercept) {
        if (verbose) print("expected p coefficients, received p+1, ignoring first")
        beta <- beta[-1]
      }
    } else if (length(beta) == p) {
      if (intercept) {
        if (verbose) print("expected p+1 coefficients, received p, reestimating intercept")
        intc <- mean(y[ind] - X[ind, ] %*% beta)
        beta <- c(intc, beta)
      }
    } else {
      stop("uninterpretable coefficients")
    }
    if (intercept) {
      X <- cbind(rep(mean(abs(X)), n), X)
      beta[1] <- beta[1] / mean(abs(X))
    } 
  } else if (family == "binomial") {
    if (!intercept) stop("Binomial family should use an intercept")
    if (length(beta) == p) {
      stop("Need to reestimate the intercept, no shortcut, dual work")
    } else if (length(beta) < p || length(beta) > p + 1) {
      stop("uninterpretable coefficients")
    }
    # transformed data
    xbh <- beta[1] + X %*% beta[-1]
    ph <- exp(xbh) / (1 + exp(xbh))
    w <- ph * (1 - ph)
    W <- diag(n)
    diag(W) <- w
    Yadj <- xbh + solve(W) %*% (y - ph)
    beta[1] <- beta[1] / mean(abs(X))
    X <- cbind(rep(mean(abs(X)), n), X)
    y <- sqrt(W) %*% Yadj
    X <- sqrt(W) %*% X
    sigma <- 1
  }
  
  chosen <-  which(abs(beta) > tol.beta) # selected variables
  s <- length(chosen)
  if (s == 0) {
    return(NULL)
  }
  if (intercept && s == 1) {
    return(NULL)
  }
  
  y1 <- y[ind]
  X1 <- X[ind, ]
  X.E <- X[, chosen] # active variables on  full data set 
  X.Ei <- ginv(X.E)
  X.E1 <- X1[ ,chosen] # active variables on selection data set 
  X.Ei1 <- ginv(X.E1)
  inv.info.E <- tcrossprod(X.Ei, X.Ei)
  inv.info.E1 <- tcrossprod(X.Ei1, X.Ei1)
  beta.E <- X.Ei %*% y # beta_hat according to OLS on full data set
  beta.E1 <- X.Ei1 %*% y1 # beta_hat according to OLS on selection data set
  z.E <- sign(beta[chosen])
  
  C <- solve(crossprod(X.E1, X.E1))

  if (intercept) {
    lam.vec <- c(0, rep(lambda, s-1))
    b1 <- -diag(x = (z.E), nrow = s) %*% C %*% (lam.vec * z.E) # linear constraint as derived in Lee et al. 2016
  } else {
    b1 <- -lambda * diag(x = (z.E), nrow = s) %*% C %*% z.E # linear constraint as derived in Lee et al. 2016 
  }
  if (selected) {
    if (n - splitn > s) {
      # in this case we sample from c(beta_hat_full, beta_hat_selection)
      # linear constraint as derived in Lee et al. 2016,
      # adjusted for the fact, that sampling is not from y1, but from a linear transformation thereof
      linear.part <- matrix(0,nrow = s,ncol = 2*s)
      linear.part[, (s + 1):(2 * s)] <- -diag(z.E)
      b <- b1 
      
      # covariance of c(beta_hat_full, beta_hat_selection)
      cov <- matrix(0, 2 * s, 2 * s)
      cov[1:s, 1:s] <- inv.info.E
      cov[(s + 1):(2 * s), 1:s] <- inv.info.E
      cov[1:s, (s + 1):(2 * s)] <- inv.info.E
      cov[(s + 1):(2 * s), (s + 1):(2 * s)] <- inv.info.E1
      con.covariance = cov * sigma**2
      
      # for the conditional law
      # we will change the linear function for each coefficient
      selector <- matrix(0, nrow = s, ncol = 2 * s)
      selector[, 1:s]  <- diag(s)
      # LHS of the additional equality constraint needed to sample for a specific coefficients
      # as explained in chapter 'Selective Inference for Linear Regression' of Fithian et al.
      conditional.linear <- crossprod(X.E, X.E) %*% selector
      # a valid initial condition (assuming glmnet fulfills KKT, not always the case)
      initial <- c(beta.E, beta.E1) 
      OLS.func <- selector 
    } else {
      # in this case we sample from c(beta_hat_selection, y2)
      # linear constraint as derived in Lee et al. 2016,
      # adjusted for the fact, that sampling is not from y1, but from a linear transformation thereof
      linear.part <- matrix(0, nrow = s, ncol = s + n - splitn)
      linear.part[, 1:s] <- -diag(x = (z.E), nrow = s)
      b <- b1
      
      # specify covariance of Gaussian vector c(beta_hat_selection, y2)
      # beta_hat_selection and y2 are independent, so only diagonal blocks are non-empty
      cov <- matrix(0, nrow = (s + n - splitn), ncol = (s + n - splitn))
      cov[1:s, 1:s] <- inv.info.E1
      cov[(s + 1):(s + n - splitn), (s + 1):(s + n - splitn)] <- diag(n - splitn) 
      con.covariance <- cov * sigma ** 2
      
      # LHS of the additional equality constraint needed to sample for a specific coefficients
      # as explained in chapter 'Selective Inference for Linear Regression' of Fithian et al.
      conditional.linear <- matrix(0, nrow = s, ncol = (s + n - splitn))
      conditional.linear[, 1:s]  <- ginv(inv.info.E1)
      conditional.linear[,(s + 1):(s + n - splitn)] <- t(X[-ind, chosen])
      OLS.func <- inv.info.E %*% conditional.linear
      
      # a valid initial condition (assuming glmnet fulfills KKT, not always the case)
      initial = c(beta.E1, y[-ind]) 
    }
    
    pvalues <- numeric(s)
    ndraw0 <- ndraw
    burnin0 <- burnin
    if (FWER) {
      if (intercept) {
        sig.level <- sig.level/(s-1)
      } else {
        sig.level <- sig.level/s
      }
    }
    correction.factor <- (1 - log(aggregation)) / aggregation
    min.size <- ceiling(1 / (sig.level / correction.factor))
    nconstraints <- s - 1
    DOF <- length(initial) - nconstraints
    if (DOF > 15) {
      ndraw <- ceiling(ndraw0 * DOF / 15)
      burnin <- ceiling(burnin0 * DOF / 15)
      if (verbose) {
        print(paste("Increasing ndraw to ", ndraw, "and burnin to ", burnin, "due to degrees of freedom"))
      }
    }
    
    nskipped <- 0 
    for (j in 1:s) {
      if (intercept && j == 1) next ()
      if (!is.null(which.check)){
        if ((intercept && !((chosen[j]-1) %in% which.check)) || (!intercept && !((chosen[j]) %in% which.check))){
          pvalues[j] <- 1
          next()
        }
      }
      if (intercept) {
        if (verbose) print(paste("Inference for", chosen[j]-1))
      } else {
        if (verbose) print(paste("Inference for", chosen[j])) 
      }
      
      
      eta <- OLS.func[j,]
      # direction to project samples on to afterwards

      # new mean and covariance adjusting for additional equality constraint
      if (s > 1) {
        conditional.law <- conditional(con.covariance, rep(0, dim(conditional.linear)[2]), conditional.linear[-j, ], 
                                      crossprod(X.E, y)[-j])
        conditional.law.covariance <- con.covariance - conditional.law$delta.cov
        conditional.law.mean <- -conditional.law$delta.mean
      } else {
        conditional.law.covariance <- con.covariance
        conditional.law.mean <- rep(0, dim(conditional.linear)[2])
      }
      
      skip <<- FALSE # indicator whether Hamiltonian sampler shall be skipped right away
      ft <<- TRUE # indicator whether it is the first chain for the given covariate
      # both indicators are to be shared with other functions
      
      white.out <- whiten(conditional.law.covariance, linear.part, b, conditional.law.mean, rank = DOF)
      forward.map <- white.out$forward.map
      inverse.map <- white.out$inverse.map
      new.A <- white.out$new.A
      new.b <- white.out$new.b
      white.Y <- forward.map(initial)
      white.eta <- forward.map(conditional.law.covariance %*% eta)
      if (max (new.A %*% white.Y - new.b) > 0) stop("Constraints not fulfilled after whitening")
      
      # get a sample of points fulfilling all constraints
      Z <- sample.from.constraints(new.A, new.b,
                                   white.Y, white.eta, ndraw = ndraw, burnin = burnin,
                                   how.often = 10, verbose = verbose, time.constant = time.constant)
      Z <- t(inverse.map(t(Z)))
      continue <- FALSE
      i <- 0
      while (!continue) {
        i <- i + 1
        continue <- TRUE
        # project sampled points on to "direction of interest" to calculate p-value
        null.statistics <- Z %*% eta
        if (i > 5) break()
        observed <- sum(initial * eta)
        lennull <- length(null.statistics)
        # calculate p-values on two halves of samples to estimate convergence
        if (beta[chosen[j]] > 0) {
          pval1 <- (sum(null.statistics[1:floor(lennull / 2)] >= observed) + 1) / (floor(lennull / 2) + 1)
          pval2 <- (sum(null.statistics[(floor(lennull / 2) + 1):lennull] >= observed) + 1) / (ceiling(lennull / 2) + 1)
        } else {
          pval1 <- (sum(null.statistics[1:floor(lennull / 2)] <= observed) + 1) / (floor(lennull / 2) + 1)
          pval2 <- (sum(null.statistics[(floor(lennull / 2) + 1):lennull] <= observed) + 1) / (ceiling(lennull / 2) + 1)
        }
        if (min(pval1, pval2) < sig.level && ft) {
          # if potentially significant
          if (verbose) print("Checking significance with longer chain")
          ft <<- FALSE
          if (verbose) {
            print(paste("Increasing ndraw to ", 2 * max(min.size, ndraw), "and burnin to ", ndraw))
          }
          Z <- sample.from.constraints(new.A, new.b, forward.map(Z[lennull, ]), white.eta,
                                       ndraw = 2 * max(min.size, ndraw), burnin = 0,
                                       how.often = 10, verbose = verbose, time.constant = time.constant)
          Z <- t(inverse.map(t(Z)))
          continue <- FALSE
        } else if (min(pval1, pval2) < sig.level && max(pval1, pval2) > 1.5 * sig.level) {
          # if significance is contradicting on the two parts
          if (verbose) print("Appending the chain")
          ft <<- FALSE
          if (verbose) {
            print(paste("Adding", 2 * max(min.size, ndraw), "draws"))
          }
          Zn <- sample.from.constraints(new.A, new.b, forward.map(Z[lennull, ]), white.eta,
                                        ndraw = 2 * max(min.size, ndraw), burnin = 0,
                                        how.often = 10, verbose = verbose, time.constant = time.constant)
          
          Zn <- t(inverse.map(t(Zn)))
          Z <- rbind(Z, Zn)
          continue <- FALSE
        }
      }
      
      add <- 1 # add <- 0 would allow for 0 pvalues
      if (beta[chosen[j]] > 0) {
        pvalues[j] <- (sum(null.statistics >= observed) + add) / (length(null.statistics) + add)
      } else {
        pvalues[j] <- (sum(null.statistics <= observed) + add) / (length(null.statistics) + add)
      }
      if (skip) nskipped <- nskipped + 1
      skip <<- FALSE
    }
    if (intercept){
      warning(paste("Hamiltonian sampler failed for", nskipped, "out of", s-1, "variables"))
    } else {
      warning(paste("Hamiltonian sampler failed for", nskipped, "out of", s, "variables"))
    }
    
    if (intercept) pvalues <- pvalues[-1]
    return(list(pv = pvalues))
  } else {
    # saturated model
    X.o <- rbind(X.E[ind, ], X.E[-ind, ])
    y.o <- c(y[ind], y[-ind])
    A1 <- -diag(z.E) %*% C %*% t(X.E1)
    A <- cbind(A1, matrix(0, nrow = dim(A1)[1], ncol = n - length(ind)))
    if (max(A %*% y.o - b1) > 0) stop("constraint error")
    Xinv <- ginv(X.o)
    pvalues <-  numeric(s)
    vlo <- numeric(s)
    vup <- numeric(s)
    estimates <- numeric(s)
    ses <- numeric(s)
    for (j in 1:s) {
      if (intercept && j==1) next ()
      if (intercept) {
        if (verbose) print(paste("Inference for", chosen[j]-1))
      } else {
        if (verbose) print(paste("Inference for", chosen[j])) 
      }
      
      
      eta <- Xinv[j, ]
      Pnoteta <- diag(length(y.o)) - eta %*% t(eta) / (norm(eta, "2") ^ 2)
      z <- Pnoteta %*% y.o
      c <- eta / (norm(eta, "2") ^ 2)
      posind <- which(A %*% c > 0)
      negind <- which(A %*% c < 0)
      V <- (b1 - A %*% z) / (A %*% c)
      Vplus <- min(V[posind])
      Vmin <- max(V[negind])
      if (Vplus < Vmin) stop("Boundary error")
      sigeta <- sigma * norm(eta, "2")
      etay <- sum(eta * y.o)
      if (etay < Vmin || etay > Vplus) stop("Projection not in range")
      vlo[j] <- Vmin
      vup[j] <- Vplus
      estimates[j] <- etay
      ses[j] <- sigeta
      
      pv <- selectiveInference:::tnorm.surv(etay, 0, sigeta, Vmin, Vplus)
      if (pv == 0 || pv == 1 || is.na(pv)) {
        pv <- selectiveInference:::tnorm.surv(etay, 0, sigeta, Vmin, Vplus, bits = 2)
        if (is.na(pv)) pv <- selectiveInference:::tnorm.surv(etay, 0, sigeta, Vmin, Vplus, bits = 100)
      }
      if (beta[chosen[j]] > 0) {
        pvalues[j] <- pv
      } else {
        pvalues[j] <- 1 - pv
      }
      if (is.na(pvalues[j])){
        if (intercept){
          warning(paste("p-value for", chosen[j] - 1, "was set to NA"))
        } else {
          warning(paste("p-value for", chosen[j], "was set to NA"))
        }
      } 
    }
    if (intercept) {
      pvalues <- pvalues[-1]
      vlo <- vlo[-1]
      vup <- vup[-1]
      estimates <- estimates[-1]
      ses <- ses[-1]
    } 
    # return extra information for saturated model in order to determine CI
    return(list(pv = pvalues, vlo = vlo, vup = vup, estimates = estimates, ses = ses))
  }
}

carve.lasso.group <- function(X, y, ind, groups, beta, tol.beta, lambda, sigma = NULL, family = "gaussian",
                              intercept = TRUE, which.check = NULL, ndraw = 8000, burnin = 2000,
                              sig.level = 0.05, FWER = TRUE, aggregation = 0.05, time.constant = 1e-6, verbose = FALSE) {
  # to be applied after Lasso Selection
  # X: full X matrix
  # y: full y vector
  # ind: indices used for selection, i.e. Lasso was applied to X[ind,] and y[ind]
  # groups: groups to be tested
  # beta: coefficients obtained from Lasso selection
  # tol.beta tolerance, to assume variable to be active
  # lambda: penalty parameter used for Lasso 
  # (with objective 1/2||X*\beta-y||^2+\lambda*||\beta||_1, i.e different normalization than in glmnet)
  # sigma: standard deviation, assumed to be known for Gaussian
  # family: gaussian or binomial
  # intercept: was the model fit using an intercept
  # which.check: which groups shall be checked (might skip "unimportant" ones)
  # ndraw: number of points to draw
  # burnin: number of initial points to burn
  # ndraw and burning might be increased for high degrees of freedom or due to required significance
  # the following three values are used to defined the minimally required sample size
  # sig.level: level for the hypothesis test
  # FWER: whether a FWER correction will be applied
  # aggregation: aggregation parameter \gamma.min
  # time.constant: the Hamiltonian sampler is assumed to be stuck, if it does not
  # finish after (ndraw + burnin) * number of constraints * dimensionality * time.constant seconds
  # verbose: whether to print key steps
  
  knownfamilies <- c("gaussian", "binomial")
  if (!(family %in% knownfamilies)) {
    stop(paste("Unknown family, family should be one of", paste(knownfamilies, collapse = ", ")))
  }
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (family == "gaussian") {
    if (is.null(sigma)) stop("Sigma needs to be provided for Gaussian family")
    if (length(beta) == p + 1) {
      if (!intercept) {
        if (verbose) print("expected p coefficients, received p+1, ignoring first")
        beta <- beta[-1]
      }
    } else if (length(beta) == p) {
      if (intercept) {
        if (verbose) print("expected p+1 coefficients, received p, reestimating intercept")
        intc <- mean(y[ind] - X[ind, ] %*% beta)
        beta <- c(intc, beta)
      }
    } else {
      stop("uninterpretable coefficients")
    }
    if (intercept) {
      X <- cbind(rep(mean(abs(X)), n), X)
      beta[1] <- beta[1] / mean(abs(X))
    } 
  } else if (family == "binomial") {
    if (!intercept) stop("Binomial family should use an intercept")
    if (length(beta) == p) {
      stop("Need to reestimate the intercept, no shortcut, dual work")
    } else if (length(beta) < p || length(beta) > p + 1) {
      stop("uninterpretable coefficients")
    }
    # transformed data
    xbh <- beta[1] + X %*% beta[-1]
    ph <- exp(xbh) / (1 + exp(xbh))
    w <- ph * (1 - ph)
    W <- diag(n)
    diag(W) <- w
    Yadj <- xbh + solve(W) %*% (y - ph)
    beta[1] <- beta[1] / mean(abs(X))
    X <- cbind(rep(mean(abs(X)), n), X)
    y <- sqrt(W) %*% Yadj
    X <- sqrt(W) %*% X
    sigma <- 1
  }
  
  splitn <- length(ind)
  chosen <- which(abs(beta) > tol.beta) # selected variables
  s <- length(chosen)
  if (s == 0) {
    warning("Empty model selected, returning 1")
    return(list(pv = rep(1, length(groups))))
  }
  if (intercept && s == 1 && !(0 %in% unlist(groups))) {
    warning("Only intercept selected, not checking for that, returning 1")
    return(list(pv = rep(1, length(groups))))
  }
  
  y1 <- y[ind]
  X1 <- X[ind, ]
  X.E <- X[, chosen] # active variables on  full data set 
  X.Ei <- ginv(X.E)
  X.E1 <- X1[, chosen] # active variables on selection data set 
  X.Ei1 <- ginv(X.E1)
  inv.info.E <- tcrossprod(X.Ei, X.Ei)
  inv.info.E1 <- tcrossprod(X.Ei1, X.Ei1)
  beta.E <- X.Ei %*% y # beta_hat according to OLS on full data set
  beta.E1 <- X.Ei1 %*% y1 # beta_hat according to OLS on selection data set
  z.E <- sign(beta[chosen])
  
  C <- solve(crossprod(X.E1, X.E1))
  if (intercept) {
    lam.vec <- c(0, rep(lambda, s-1))
    b1 <- -diag(x = (z.E), nrow = s) %*% C %*% (lam.vec * z.E) # linear constraint as derived in Lee et al. 2016
  } else {
    b1 <- -lambda * diag(x = (z.E), nrow = s) %*% C %*% z.E # linear constraint as derived in Lee et al. 2016 
  }
  if (n == splitn) {
    # pure post-selection inference
    # in this case we sample from c(beta_hat_selection)
    # linear constraint as derived in Lee et al. 2016,
    # adjusted for the fact, that sampling is not from y1, but from a linear transformation thereof
    linear.part <- matrix(0, nrow = s, ncol = s)
    linear.part[, 1:s] <- -diag(z.E)
    b <- b1
    
    # specify covariance of Gaussian vector c(beta_hat_selection)
    cov <- inv.info.E1
    con.covariance = cov * sigma**2
    
    # LHS of the additional equality constraint needed to sample for a specific coefficients
    # as explained in chapter 'Selective Inference for Linear Regression' of Fithian et al.
    conditional.linear = matrix(0,nrow = s, ncol = s)
    conditional.linear[, 1:s]  = ginv(inv.info.E1)
    
    OLS.func <- inv.info.E %*% conditional.linear
    # a valid initial condition
    initial <- c(beta.E1) 
  } else if (n - splitn > s) {
    # in this case we sample from c(beta_hat_full, beta_hat_selection)
    # linear constraint as derived in Lee et al. 2016,
    # adjusted for the fact, that sampling is not from y1, but from a linear transformation thereof
    linear.part <- matrix(0, nrow = s,ncol = 2 * s)
    linear.part[, (s + 1):(2 * s)] <- -diag(z.E)
    b <- b1 
    
    # specify covariance of 2s Gaussian vector
    # covariance of c(beta_hat_full, beta_hat_selection)
    cov <- matrix(0,2 * s, 2 * s)
    cov[1:s, 1:s] <- inv.info.E
    cov[(s + 1):(2 * s), 1:s] <- inv.info.E
    cov[1:s, (s + 1):(2 * s)] <- inv.info.E
    cov[(s + 1):(2 * s), (s + 1):(2 * s)] <- inv.info.E1
    con.covariance <- cov * sigma ** 2
    
    # for the conditional law
    # we will change the linear function for each coefficient
    selector <- matrix(0, nrow = s, ncol = 2 * s)
    selector[, 1:s] <- diag(s)
    # LHS of the additional equality constraint needed to sample for a specific coefficients
    # as explained in chapter 'Selective Inference for Linear Regression' of Fithian et al.
    conditional.linear <- crossprod(X.E, X.E) %*% selector
    
    # a valid initial condition (assuming glmnet fulfills KKT, not always the case)
    initial <- c(beta.E, beta.E1) 
    OLS.func <- selector 
  } else {
    # in this case we sample from c(beta_hat_selection, y2)
    # linear constraint as derived in Lee et al. 2016,
    # adjusted for the fact, that sampling is not from y1, but from a linear transformation thereof
    linear.part <- matrix(0, nrow = s, ncol = s + n - splitn)
    linear.part[, 1:s] <- -diag(z.E)
    b <- b1
    
    # specify covariance of Gaussian vector c(beta_hat_selection, y2)
    # beta_hat_selection and y2 are independent, so only diagonal blocks are non-empty
    cov <- matrix(0, nrow = (s + n - splitn), ncol = (s + n - splitn))
    cov[1:s, 1:s] <- inv.info.E1
    cov[(s + 1):(s + n - splitn), (s + 1):(s + n - splitn)] <- diag(n - splitn) 
    con.covariance <- cov * sigma ** 2
    # LHS of the additional equality constraint needed to sample for a specific coefficients
    # as explained in chapter 'Selective Inference for Linear Regression' of Fithian et al.
    conditional.linear <- matrix(0, nrow = s, ncol = (s + n - splitn))
    conditional.linear[, 1:s]  <- ginv(inv.info.E1)
    conditional.linear[,(s + 1):(s + n - splitn)] <- t(X[-ind, chosen])
    
    # write the OLS estimates of full model in terms of X.E1^{dagger}y1, y2
    OLS.func = inv.info.E%*%conditional.linear
    
    # a valid initial condition (assuming glmnet fulfills KKT, not always the case)
    initial = c(beta.E1, y[-ind]) 
  }
  pvaluessum <- numeric(length(groups))
  j <- 0
  ndraw0 <- ndraw
  burnin0 <- burnin
  
  if (intercept) {
    ngroup.tested <- sum(unlist(lapply(lapply(groups, intersect, chosen - 1), length)) > 0)
  } else {
    ngroup.tested <- sum(unlist(lapply(lapply(groups, intersect, chosen), length)) > 0)
  }
  if (FWER) {
    sig.level <- sig.level / ngroup.tested
  }
  correction.factor <- (1 - log(aggregation)) / aggregation
  min.size <- ceiling(1/(sig.level/correction.factor))
  
  nskipped <- 0
  for (group in groups) {
    j <- j + 1
    if (!is.null(which.check) && !(j %in% which.check)){
      if (verbose) print(paste("Not checking group", j))
      pvaluessum[j] <- 1
      next()
    } 
    if (intercept) {
      group.vars <- which((chosen-1) %in% group)
    } else {
      group.vars <- which(chosen %in% group)
    }
    
    if (length(group.vars) > 0) {
      if (verbose) print(paste("Inference for group", j))
      eta <- matrix(t(OLS.func[group.vars, ]), ncol = length(group.vars))
      # directions to project samples on to afterwards
      
      # new mean and covariance adjusting for additional equality constraint
      if (length(group.vars) < s) {
        conditional.law <- conditional(con.covariance, rep(0, dim(conditional.linear)[2]),
                                       matrix(conditional.linear[-group.vars, ], nrow = s -length(group.vars)),
                                       crossprod(X.E, y)[-group.vars])
        conditional.law.covariance <- con.covariance - conditional.law$delta.cov
        conditional.law.mean <- -conditional.law$delta.mean
      } else {
        # there is no additional equality constraint if all the selected variables are in the group
        conditional.law.covariance <- con.covariance
        conditional.law.mean <- rep(0,length(initial))
        
      }
      
      nconstraints <- s - length(group.vars)
      DOF <- length(initial) - nconstraints
      if (DOF > 15) {
        ndraw <- ceiling(ndraw0 * DOF / 15)
        burnin <- ceiling(burnin0 * DOF / 15)
      }
      
      skip <<- FALSE # indicator whether Hamiltonian sampler shall be skipped right away
      ft <<- TRUE # indicator whether it is the first chain for the given covariate
      # both indicators are to be shared with other functions
      
      white.out <- whiten(conditional.law.covariance, linear.part, b, conditional.law.mean, rank = DOF)
      forward.map <- white.out$forward.map
      inverse.map <- white.out$inverse.map
      new.A <- white.out$new.A
      new.b <- white.out$new.b
      white.Y <- forward.map(initial)
      white.eta <- forward.map(conditional.law.covariance %*% eta)
      if (max (new.A %*% white.Y - new.b) > 0) stop("Constraints not fulfilled after whitening")
      
      # get a sample of points fulfilling all constraints
      Z <- sample.from.constraints(new.A, new.b, white.Y, white.eta, ndraw = ndraw,
                                   burnin = burnin, how.often = 10, verbose = verbose)
      Z <- t(inverse.map(t(Z)))
      
      etasign <- t(t(eta)*sign(beta[chosen[group.vars]]))
      continue <- FALSE
      
      i <- 0
      while (!continue) {
        i <- i + 1
        continue <- TRUE
        # project sampled points on to "directions of interest" to calculate p-value
        all.null.statistics <- Z %*% etasign
        allobserved <- initial %*% etasign
        if (i > 5) break()
        sum.null.statistics <- apply(all.null.statistics, 1, sum)
        sumobserved <- sum(allobserved)
        lennull <- length(sum.null.statistics)
        # calculate p-values on two halves of samples to estimate convergence
        pval1 <- (sum(sum.null.statistics[1:floor(lennull / 2)] >= sumobserved) + 1) / (floor(lennull / 2) + 1)
        pval2 <- (sum(sum.null.statistics[(floor(lennull / 2)+1):lennull] >= sumobserved) + 1) / (ceiling(lennull / 2) + 1)
        if (min(pval1, pval2) < sig.level && ft) {
          # if potentially significant
          if (verbose) print("Checking significance with longer chain")
          ft <<- FALSE
          Z <- sample.from.constraints(new.A, new.b, forward.map(Z[lennull, ]), white.eta,
                                       ndraw = 2 * max(min.size, ndraw), burnin=0,
                                       how.often = 10, verbose = verbose)
          Z <- t(inverse.map(t(Z)))
          continue <- FALSE
        } else if (min(pval1, pval2) < sig.level && max(pval1, pval2) > 1.5 * sig.level) {
          # if significance is contradicting on the two parts
          if (verbose) print("Appending the chain")
          if (verbose) {
            print(paste("Adding", 2 * max(min.size, ndraw), "draws"))
          }
          ft <<- FALSE
          Zn <- sample.from.constraints(new.A, new.b, forward.map(Z[lennull, ]), white.eta,
                                        ndraw = 2 * max(min.size, ndraw), burnin = 0,
                                        how.often = 10, verbose = verbose)
          Zn <- t(inverse.map(t(Zn)))
          Z <- rbind(Z, Zn)
          continue <- FALSE
        }
      }
      
      
      # project sampled points on to "directions of interest" to calculate p-value
      all.null.statistics <- Z %*% etasign
      allobserved <- initial %*% etasign
      sum.null.statistics <- apply(all.null.statistics, 1, sum)
      sumobserved <- sum(allobserved)
      add <- 1 # add <- 0 would allow for 0 pvalues
      pvaluessum[j] <- (sum(sum.null.statistics >= sumobserved) + add) / (length(sum.null.statistics) + add)
    } else {
      if (verbose) print(paste("Not checking group", j))
      pvaluessum[j] <- 1
    }
    if (skip) nskipped <- nskipped + 1
    skip <<- FALSE
  }
  warning(paste("Hamiltonian sampler failed for", nskipped, "out of", ngroup.tested, "groups"))
  return(list(pv = pvaluessum))
}

conditional<-function(S, mean, C, d) {
  # Return an equivalent constraint 
  # after having conditioned on a linear equality.
  # Let the inequality constraints be specified by
  # `(A,b)` and the equality constraints be specified
  # by `(C,d)`. We form equivalent inequality constraints by 
  # considering the residual
  # AY - E(AY|CY=d)
  # delta.cov and delta.mean are reduction in covariance and mean, due to equality constraint
  # formula can be generally derived for any such equality constraint
  C <- matrix(C, ncol = nrow(S))
  M1 <- tcrossprod(S, C)
  M2 <- C %*% M1
  if (is.matrix(M2)) {
    M2i <- ginv(M2)
    delta.cov <- M1 %*% tcrossprod(M2i, M1)
    delta.mean <- M1 %*% M2i %*% (C %*% mean - d)
  } else {
    M2i <- 1 / M2
    delta.cov <- M1 %o% M1 / M2i
    delta.mean <- M1 * d  / M2i 
  }
  return(list(delta.cov = delta.cov, delta.mean = delta.mean))
}


whiten <- function(cov, linear.part, b, mmean, rank = NULL) {
  #   Return a whitened version of constraints in a different
  #   basis, and a change of basis matrix.

  # calculate a root of the covariance matrix using EVD
  if (is.null(rank)) rank <- rankMatrix(cov)[1]
  ev <- eigen(cov)
  D1 <- ev$values
  U <- ev$vectors
  Dtry <- tryCatch_W_E(sqrt(D1[rank:1]))
  while (!is.null(Dtry$warning)) {
    warning("Reducing Rank")
    rank <- rank - 1
    Dtry <- tryCatch_W_E(sqrt(D1[rank:1]))
  }
  D <- Dtry$value
  U <- U[, rank:1]
  
  sqrt.cov <- Re(t(t(U)*D))
  sqrt.inv <- Re(t(U) / D)
  # get equivalent constraint for whitened vector i.e linear.part%*%y<b <=> new.A%*%y_white<new.b
  new.A <- linear.part %*% sqrt.cov
  den <- sqrt(apply(new.A ^ 2, 1, sum))
  new.b <- b - linear.part %*% mmean
  new.A <- new.A / den
  new.b <- new.b / den
  mu <- mmean
  # colour a "white" point
  inverse.map <- function(Z) {
    (sqrt.cov %*% Z) + as.vector(mu)
  } 
  # whiten a "coloured" point
  forward.map <- function(W) sqrt.inv %*% (W - as.vector(mu))
  return(list(inverse.map = inverse.map, forward.map = forward.map, new.A = new.A, new.b = new.b))
}


sample.from.constraints <- function(new.A, new.b, white.Y, white.direction.of.interest,
                            how.often = -1, ndraw = 1000, burnin = 1000, white = FALSE,
                            use.constraint.directions = TRUE, verbose = FALSE, time.constant = 1e-6) {
  # routine to do whitening, activate the sampler, and recolour the samples
  if (how.often < 0) {
    how.often <- ndraw + burnin
  }
  
  if (skip) {
    # if Hamiltonian sampler got stuck before for same set-up do not try again
    #  sample from whitened points with new constraints
    Z <- sample.truncnorm.white(new.A, new.b, white.Y, white.direction.of.interest,
                                           how.often = how.often, ndraw = ndraw, burnin = burnin,
                                           sigma = 1, use.A = use.constraint.directions)
  } else {
    # this sampler seems to work better for most cases, though, it sometime takes "forever" => not usable
    # current wrap around is not windows supported
    nw <- length(white.Y)
    nconstraint <- dim(new.A)[1]
    nsample <- ndraw / 2 + burnin
    time.factor <- nw * nconstraint * nsample
    time.limit <- time.factor * time.constant
    tic()
    trywhite <- tryCatch_W_E(eval_with_timeout({rtmg(ndraw / 2, diag(nw), rep(0,nw), white.Y,
                                                     -new.A, as.vector(new.b), burn.in = burnin)},
                                               timeout = time.limit, on_timeout = "error"), 0)
    time <- toc(quiet = TRUE)
    time.diff <- round(time$toc - time$tic, 4)
    if (!is.null(trywhite$error) || !is.matrix(trywhite$value)) {
      skip <<- TRUE
      if (ft){
        first.text <- "this variable was tested for the first time;"
      } else {
        first.text <- "this variable was not tested for the first time;"
      }
      warning(paste("Hamiltonian not successful after", time.diff, "for", nsample, "samples,", nw, "dimensions and", nconstraint, "constraints"))
      warning(paste("Evaluation of Hamiltonian sampler not successful:", trywhite$error, first.text, "using hit-and-run sampler"))
      #  sample from whitened points with new constraints
      Z <- sample.truncnorm.white(new.A, new.b, white.Y, white.direction.of.interest,
                                              how.often = how.often, ndraw = ndraw, burnin = burnin,
                                              sigma = 1, use.A = use.constraint.directions)
    } else {
      warning(paste("Hamiltonian successful after", time.diff, "for", nsample, "samples,", nw, "dimensions and", nconstraint, "constraints"))
      Z <- trywhite$value
    }
    
  }
  return(Z)
}


constraint.checker<-function(X1, y1, beta, tol.beta, lambda,
                            family = "gaussian", intercept = TRUE, verbose = FALSE) {
  # to be applied after Lasso Selection, checks if Lasso convergence is sufficient to run the carving routine
  # X1: selection part of X
  # y1: selection part of y
  # beta: coefficients obtained from Lasso selection
  # sigma: standard deviation, assumed to be known for Gaussian
  # tol.beta tolerance, to assume variable to be active
  # lambda: penalty parameter used for Lasso 
  # (with objective 1/2||X*\beta-y||^2+\lambda*||\beta||_1, i.e different normalization than in glmnet)
  # family: Gaussian or binomial
  knownfamilies <- c("gaussian", "binomial")
  if (!family %in% knownfamilies) {
    stop(paste("Unknown family, family should be one of", paste(knownfamilies, collapse = ", ")))
  }
  n <- dim(X1)[1]
  p <- dim(X1)[2]
  if (family=="gaussian") {
    if (length(beta) == p + 1) {
      if (!intercept) {
        if (verbose) print("expected p coefficients, received p+1, ignoring first")
        beta <- beta[-1]
      }
    } else if (length(beta) == p) {
      if (intercept) {
        if (verbose) print("expected p+1 coefficients, received p, reestimating intercept")
        intc <- mean(y1 - X1 %*% beta)
        beta <- c(intc, beta)
      }
    } else {
      stop("uninterpretable coefficients")
    }
    
    if (intercept) {
      beta[1] <- beta[1] / mean(abs(X1))
      X1 <- cbind(rep(mean(abs(X1)), n), X1)
    } 
  } else if (family == "binomial") {
    if (!intercept) stop("Binomial family should use an intercept")
    if (length(beta) == p) {
      stop("Need to reestimate the intercept, no shortcut, dual work")
    } else if (length(beta) < p || length(beta) > p + 1) {
      stop("uninterpretable coefficients")
    }
    # transformed data
    xbh <- beta[1] + X1 %*% beta[-1]
    ph <- exp(xbh) / (1 + exp(xbh))
    w <- ph * (1 - ph)
    W <- diag(n)
    diag(W) <- w
    Yadj <- xbh + solve(W) %*% (y1 - ph)
    beta[1] <- beta[1] / mean(abs(X1))
    X1 <- cbind(rep(mean(abs(X1)), n), X1)
    y1 <- sqrt(W) %*% Yadj
    X1 <- sqrt(W) %*% X1
  }
 
  chosen <- which(abs(beta) > tol.beta) # selected variables
  s <- length(chosen)
  if ((s == 0 && !intercept) || (s == 1 && intercept)) {
    warning("Checking empty model, by default TRUE since only active constraints are checked")
    return(TRUE)
  }

  
  X.E1 <- X1[, chosen] # active variables on selection data set 
  
  X.Ei1 <- ginv(X.E1)
  
  z.E <- sign(beta[chosen])
  C <- solve(crossprod(X.E1, X.E1))
  A1 <- -diag(x = (z.E), nrow = s) %*% C %*% t(X.E1)
  if (intercept) {
    lam.vec <- c(0, rep(lambda, s - 1))
    b1 <- -diag(x = (z.E), nrow = s) %*% C %*% (lam.vec * z.E) # linear constraint as derived in Lee et al. 2016
  } else {
    b1 <- -lambda * diag(x = (z.E), nrow = s) %*% C %*% z.E # linear constraint as derived in Lee et al. 2016 
  }
  if (max(A1 %*% y1 - b1) < 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

eval_with_timeout <- function(expr, envir = parent.frame(), timeout, on_timeout = c("error", "warning", "silent")) {
  # substitute expression so it is not executed as soon it is used
  expr <- substitute(expr)
  
  # match on_timeout
  on_timeout <- match.arg(on_timeout)
  
  # execute expr in separate fork
  myfork <- parallel::mcparallel({
    eval(expr, envir = envir)
  }, silent = FALSE)
  
  # wait max n seconds for a result.
  myresult <- parallel::mccollect(myfork, wait = FALSE, timeout = timeout)
  # kill fork after collect has returned
  tools::pskill(myfork$pid, tools::SIGKILL)
  tools::pskill(-1 * myfork$pid, tools::SIGKILL)
  
  # clean up:
  parallel::mccollect(myfork, wait = FALSE)
  
  # timeout?
  if (is.null(myresult)) {
    if (on_timeout == "error") {
      stop("reached elapsed time limit")
    } else if (on_timeout == "warning") {
      warning("reached elapsed time limit")
    } else if (on_timeout == "silent") {
      myresult <- NA
    }
  }
  
  # move this to distinguish between timeout and NULL returns
  myresult <- myresult[[1]]
  
  if ("try-error" %in% class(myresult)) {
    stop(attr(myresult, "condition"))
  }
  
  # send the buffered response
  return(myresult)
}