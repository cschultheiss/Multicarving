# this file contains the necessary routines to execute data carving after model selection
# either for single variable testing or group testing
# the sampler is defined elsewhere


OptimalFixedLasso<-function(X, y, ind, beta, sigma = NULL, tol.beta, lambda, family = "gaussian",
                            intercept = TRUE, ndraw = 8000, burnin = 2000, sig_Level = 0.05,
                            FWER = TRUE, aggregation = 0.05, selected = TRUE, verbose = FALSE, which.check = NULL) {
  # to be applied after Lasso Selection
  # X: full X matrix
  # y: full y vector
  # ind: indices used for selection, i.e. Lasso was applied to X[ind, ] and y[ind]
  # beta: coefficients obtained from Lasso selection
  # sigma: standard deviation, assumed to be known for Gaussian
  # tol.beta tolerance, to assume variable to be active
  # lambda: penalty parameter used for Lasso 
  # (with objective 1/2||X*\beta-y||^2+\lambda*||\beta||_1, i.e different normalization than in glmnet)
  # family: Gaussian or binomial
  # ndraw: number of points to draw
  # burnin: number of initial points to burn
  # ndraw and burning might be increased for high the degrees of freedom or due to required significance
  # the following three values are used to defined the minimally required sample size
  # sig_level: level for the hypothesis test
  # FWER: whether a FWER correction will be applied
  # aggregation: aggregation parameter \lambda_min
  # selected: whether to use the selected viewpoint for aggregation
  # verbose: print some key steps
  
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
  X_E <- X[, chosen] # active variables on  full data set 
  X_Ei <- ginv(X_E)
  X_E1 <- X1[ ,chosen] # active variables on selection data set 
  X_Ei1 <- ginv(X_E1)
  inv_info_E <- tcrossprod(X_Ei, X_Ei)
  inv_info_E1 <- tcrossprod(X_Ei1, X_Ei1)
  beta_E <- X_Ei %*% y # beta_hat according to OLS on full data set
  beta_E1 <- X_Ei1 %*% y1 # beta_hat according to OLS on selection data set
  z_E <- sign(beta[chosen])
  
  C <- solve(crossprod(X_E1, X_E1))

  if (intercept) {
    lamvec <- c(0, rep(lambda, s-1))
    b1 <- -diag(x = (z_E), nrow = s) %*% C %*% (lamvec * z_E) # linear constraint as derived in Lee et al. 2016
  } else {
    b1 <- -lambda * diag(x = (z_E), nrow = s) %*% C %*% z_E # linear constraint as derived in Lee et al. 2016 
  }
  if (selected) {
    if (n - splitn > s) {
      # in this case we sample from c(beta_hat_full, beta_hat_selection)
      # linear constraint as derived in Lee et al. 2016,
      # adjusted for the fact, that sampling is not from y1, but from a linear transformation thereof
      linear_part <- matrix(0,nrow = s,ncol = 2*s)
      linear_part[, (s + 1):(2 * s)] <- -diag(z_E)
      b <- b1 
      
      # covariance of c(beta_hat_full, beta_hat_selection)
      cov <- matrix(0, 2 * s, 2 * s)
      cov[1:s, 1:s] <- inv_info_E
      cov[(s + 1):(2 * s), 1:s] <- inv_info_E
      cov[1:s, (s + 1):(2 * s)] <- inv_info_E
      cov[(s + 1):(2 * s), (s + 1):(2 * s)] <- inv_info_E1
      con.covariance = cov * sigma**2
      
      # for the conditional law
      # we will change the linear function for each coefficient
      selector <- matrix(0, nrow = s, ncol = 2 * s)
      selector[, 1:s]  <- diag(s)
      # LHS of the additional equality constraint needed to sample for a specific coefficients
      # as explained in chapter 'Selective Inference for Linear Regression' of Fithian et al.
      conditional_linear <- crossprod(X_E, X_E) %*% selector
      # a valid initial condition (assuming glmnet fulfills KKT, not always the case)
      initial <- c(beta_E, beta_E1) 
      OLS_func <- selector 
    } else {
      # in this case we sample from c(beta_hat_selection, y2)
      # linear constraint as derived in Lee et al. 2016,
      # adjusted for the fact, that sampling is not from y1, but from a linear transformation thereof
      linear_part <- matrix(0, nrow = s, ncol = s + n - splitn)
      linear_part[, 1:s] <- -diag(x = (z_E), nrow = s)
      b <- b1
      
      # specify covariance of Gaussian vector c(beta_hat_selection, y2)
      # beta_hat_selection and y2 are independent, so only diagonal blocks are non-empty
      cov <- matrix(0, nrow = (s + n - splitn), ncol = (s + n - splitn))
      cov[1:s, 1:s] <- inv_info_E1
      cov[(s + 1):(s + n - splitn), (s + 1):(s + n - splitn)] <- diag(n - splitn) 
      con.covariance <- cov * sigma ** 2
      
      # LHS of the additional equality constraint needed to sample for a specific coefficients
      # as explained in chapter 'Selective Inference for Linear Regression' of Fithian et al.
      conditional_linear <- matrix(0, nrow = s, ncol = (s + n - splitn))
      conditional_linear[, 1:s]  <- ginv(inv_info_E1)
      conditional_linear[,(s + 1):(s + n - splitn)] <- t(X[-ind, chosen])
      OLS_func <- inv_info_E %*% conditional_linear
      
      # a valid initial condition (assuming glmnet fulfills KKT, not always the case)
      initial = c(beta_E1, y[-ind]) 
    }
    
    pvalues <- numeric(s)
    ndraw0 <- ndraw
    burnin0 <- burnin
    if (FWER) {
      if (intercept) {
        sig_Level <- sig_Level/(s-1)
      } else {
        sig_Level <- sig_Level/s
      }
    }
    correction_factor <- (1 - log(aggregation)) / aggregation
    min_size <- ceiling(1 / (sig_Level / correction_factor))
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
      
      
      eta <- OLS_func[j,]
      # direction to project samples on to afterwards

      # new mean and covariance adjusting for additional equality constraint
      if (s > 1) {
        conditional_law <- conditional(con.covariance, rep(0, dim(conditional_linear)[2]), conditional_linear[-j, ], 
                                      crossprod(X_E, y)[-j])
        conditional_law.covariance <- con.covariance - conditional_law$delta_cov
        conditional_law.mean <- -conditional_law$delta_mean
      } else {
        conditional_law.covariance <- con.covariance
        conditional_law.mean <- rep(0, dim(conditional_linear)[2])
      }
      
      skip <<- FALSE # indicator whether Hamiltonian sampler shall be skipped right away
      ft <<- TRUE # indicator whether it is the first chain for the given covariate
      # both indicators are to be shared with other functions
      
      white_out <- whiten(conditional_law.covariance, linear_part, b, conditional_law.mean)
      forward_map <- white_out$forward_map
      inverse_map <- white_out$inverse_map
      new_A <- white_out$new_A
      new_b <- white_out$new_b
      white_Y <- forward_map(initial)
      white_eta <- forward_map(conditional_law.covariance %*% eta)
      if (max (new_A %*% white_Y - new_b) > 0) stop("Constraints not fulfilled after whitening")
      
      # get a sample of points fulfilling all constraints
      Z <- sample_from_constraints(new_A, new_b,
                                   white_Y, white_eta, ndraw = ndraw, burnin = burnin,
                                   how_often = 10, verbose = verbose)
      Z <- t(inverse_map(t(Z)))
      continue <- FALSE
      i <- 0
      while (!continue) {
        i <- i + 1
        continue <- TRUE
        # project sampled points on to "direction of interest" to calculate p-value
        null_statistics <- Z %*% eta
        if (i > 5) break()
        observed <- sum(initial * eta)
        lennull <- length(null_statistics)
        # calculate p-values on two halves of samples to estimate convergence
        if (beta[chosen[j]] > 0) {
          pval1 <- (sum(null_statistics[1:floor(lennull / 2)] >= observed) + 1) / (floor(lennull / 2) + 1)
          pval2 <- (sum(null_statistics[(floor(lennull / 2) + 1):lennull] >= observed) + 1) / (ceiling(lennull / 2) + 1)
        } else {
          pval1 <- (sum(null_statistics[1:floor(lennull / 2)] <= observed) + 1) / (floor(lennull / 2) + 1)
          pval2 <- (sum(null_statistics[(floor(lennull / 2) + 1):lennull] <= observed) + 1) / (ceiling(lennull / 2) + 1)
        }
        if (min(pval1, pval2) < sig_Level && ft) {
          # if potentially significant
          if (verbose) print("Checking significance with longer chain")
          ft <<- FALSE
          if (verbose) {
            print(paste("Increasing ndraw to ", 2 * max(min_size, ndraw), "and burnin to ", ndraw))
          }
          Z <- sample_from_constraints(new_A, new_b, forward_map(Z[lennull, ]), white_eta,
                                       ndraw = 2 * max(min_size, ndraw), burnin = 0,
                                       how_often = 10, verbose = verbose)
          Z <- t(inverse_map(t(Z)))
          continue <- FALSE
        } else if (min(pval1, pval2) < sig_Level && max(pval1, pval2) > 1.5 * sig_Level) {
          # if significance is contradicting on the two parts
          if (verbose) print("Appending the chain")
          ft <<- FALSE
          if (verbose) {
            print(paste("Adding", 2 * max(min_size, ndraw), "draws"))
          }
          Zn <- sample_from_constraints(new_A, new_b, forward_map(Z[lennull, ]), white_eta,
                                        ndraw = 2 * max(min_size, ndraw), burnin = 0,
                                        how_often = 10, verbose = verbose)
          
          Zn <- t(inverse_map(t(Zn)))
          Z <- rbind(Z, Zn)
          continue <- FALSE
        }
      }
      
      add <- 1 # add <- 0 would allow for 0 pvalues
      if (beta[chosen[j]] > 0) {
        pvalues[j] <- (sum(null_statistics >= observed) + add) / (length(null_statistics) + add)
      } else {
        pvalues[j] <- (sum(null_statistics <= observed) + add) / (length(null_statistics) + add)
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
    X_o <- rbind(X_E[ind, ], X_E[-ind, ])
    y_o <- c(y[ind], y[-ind])
    A1 <- -diag(z_E) %*% C %*% t(X_E1)
    A <- cbind(A1, matrix(0, nrow = dim(A1)[1], ncol = n - length(ind)))
    if (max(A %*% y_o - b1) > 0) stop("constraint error")
    Xinv <- ginv(X_o)
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
      Pnoteta <- diag(length(y_o)) - eta %*% t(eta) / (norm(eta, "2") ^ 2)
      z <- Pnoteta %*% y_o
      c <- eta / (norm(eta, "2") ^ 2)
      posind <- which(A %*% c > 0)
      negind <- which(A %*% c < 0)
      V <- (b1 - A %*% z) / (A %*% c)
      Vplus <- min(V[posind])
      Vmin <- max(V[negind])
      if (Vplus < Vmin) stop("Boundary error")
      sigeta <- sigma * norm(eta, "2")
      etay <- sum(eta * y_o)
      if (etay < Vmin || etay > Vplus) stop("Projection not in range")
      vlo[j] <- Vmin
      vup[j] <- Vplus
      estimates[j] <- etay
      ses[j] <- sigeta
      if (beta[chosen[j]] > 0) {
        pvalues[j] <- selectiveInference:::tnorm.surv(etay, 0, sigeta, Vmin, Vplus)
      } else {
        pvalues[j] <- 1-selectiveInference:::tnorm.surv(etay, 0, sigeta, Vmin, Vplus)
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


conditional<-function(S, mean, C, d) {
  # Return an equivalent constraint 
  # after having conditioned on a linear equality.
  # Let the inequality constraints be specified by
  # `(A,b)` and the equality constraints be specified
  # by `(C,d)`. We form equivalent inequality constraints by 
  # considering the residual
  # AY - E(AY|CY=d)
  # delta_cov and delta_mean are reduction in covariance and mean, due to equality constraint
  # formula can be generally derived for any such equality constraint
  C <- matrix(C, ncol = nrow(S))
  M1 <- tcrossprod(S, C)
  M2 <- C %*% M1
  if (is.matrix(M2)) {
    M2i <- ginv(M2)
    delta_cov <- M1 %*% tcrossprod(M2i, M1)
    delta_mean <- M1 %*% M2i %*% (C %*% mean - d)
  } else {
    M2i <- 1 / M2
    delta_cov <- M1 %o% M1 / M2i
    delta_mean <- M1 * d  / M2i 
  }
  return(list(delta_cov = delta_cov, delta_mean = delta_mean))
}


whiten <- function(cov, linear_part, b, mmean) {
  #   Return a whitened version of constraints in a different
  #   basis, and a change of basis matrix.

  # calculate a root of the covariance matrix using EVD
  rank <- rankMatrix(cov)[1]
  ev <- eigen(cov)
  D1 <- ev$values
  # rank <- sum(abs(D1)> 1e-5 * abs(D1[1]))
  U <- ev$vectors
  Dtry <- tryCatch_W_E(sqrt(D1[rank:1]))
  while (!is.null(Dtry$warning)) {
    warning("Reducing Rank")
    rank <- rank - 1
    Dtry <- tryCatch_W_E(sqrt(D1[rank:1]))
  }
  D <- Dtry$value
  U <- U[, rank:1]
  
  sqrt_cov <- Re(t(t(U)*D))
  sqrt_inv <- Re(t(U) / D)
  # get equivalent constraint for whitened vector i.e linear_part%*%y<b <=> new_A%*%y_white<new_b
  new_A <- linear_part %*% sqrt_cov
  den <- sqrt(apply(new_A ^ 2, 1, sum))
  new_b <- b - linear_part %*% mmean
  new_A <- new_A / den
  new_b <- new_b / den
  mu <- mmean
  # colour a "white" point
  inverse_map <- function(Z) {
    (sqrt_cov %*% Z) + as.vector(mu)
  } 
  # whiten a "coloured" point
  forward_map <- function(W) sqrt_inv %*% (W - as.vector(mu))
  return(list(inverse_map = inverse_map, forward_map = forward_map, new_A = new_A, new_b = new_b))
}


sample_from_constraints <- function(new_A, new_b, white_Y, white_direction_of_interest,
                            how_often = -1, ndraw = 1000, burnin = 1000, white = FALSE,
                            use_constraint_directions = TRUE, verbose = FALSE) {
  # routine to do whitening, activate the sampler, and recolour the samples
  if (how_often < 0) {
    how_often <- ndraw + burnin
  }
  
  if (skip) {
    # if Hamiltonian sampler got stuck before for same set-up do not try again
    #  sample from whitened points with new constraints
    Z <- sample_truncnorm_white(new_A, new_b, white_Y, white_direction_of_interest,
                                           how_often = how_often, ndraw = ndraw, burnin = burnin,
                                           sigma = 1, use_A = use_constraint_directions)
  } else {
    # this sampler seems to work better for most cases, though, it sometime takes "forever" => not usable
    # current wrap around is not windows supported
    nw <- length(white_Y)
    trywhite <- tryCatch_W_E(eval_with_timeout({rtmg(ndraw / 2, diag(nw), rep(0,nw), white_Y,
                                                     -new_A, as.vector(new_b), burn.in = burnin)},
                                               timeout = 6, on_timeout = "error"), 0)
    if (!is.null(trywhite$error) || !is.matrix(trywhite$value)) {
      skip <<- TRUE
      if (ft){
        first_text <- "this variable was tested for the first time"
      } else {
        first_text <- "this variable was not tested for the first time"
      }
      warning(paste("Evaluation of Hamiltonian sampler not successful:", trywhite$error, first_text, "using hit-and-run sampler"))
      #  sample from whitened points with new constraints
      Z <- sample_truncnorm_white(new_A, new_b, white_Y, white_direction_of_interest,
                                              how_often = how_often, ndraw = ndraw, burnin = burnin,
                                              sigma = 1, use_A = use_constraint_directions)
    } else {
      Z <- trywhite$value
    }
    
  }
  return(Z)
}


OptimalFixedLassoGroup <- function(X, y, ind, beta, sigma = NULL, tol.beta, lambda, family = "gaussian", groups,
                                   intercept = TRUE, ndraw = 8000, burnin = 2000,
                                   sig_Level = 0.05, aggregation = 0.05, FWER = TRUE, verbose = FALSE, which.check = NULL) {
  # to be applied after Lasso Selection
  # X: full X matrix
  # y: full y vector
  # ind: indices used for selection, i.e. Lasso was applied to X[ind,] and y[ind]
  # beta: coefficients obtained from Lasso selection
  # sigma: standard deviation, assumed to be known
  # tol.beta tolerance, to assume variable to be active
  # lambda: penalty parameter used for Lasso 
  # (with objective 1/2||X*\beta-y||^2+\lambda*||\beta||_1, i.e different normalization than in glmnet)
  # family: gaussian or binomial
  # groups: groups to be tested
  # ndraw: number of points to draw
  # burnin: number of initial points to burn
  # ndraw and burning might be increased for high the degrees of freedom or due to required significance
  # the following three values are used to defined the minimally required sample size
  # sig_level: level for the hypothesis test
  # FWER: whether a FWER correction will be applied
  # aggregation: aggregation parameter \lambda_min
  # selected: whether to use the selected viewpoint for aggregation
  # verbose: print some key steps
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
  X_E <- X[, chosen] # active variables on  full data set 
  X_Ei <- ginv(X_E)
  X_E1 <- X1[, chosen] # active variables on selection data set 
  X_Ei1 <- ginv(X_E1)
  inv_info_E <- tcrossprod(X_Ei, X_Ei)
  inv_info_E1 <- tcrossprod(X_Ei1, X_Ei1)
  beta_E <- X_Ei %*% y # beta_hat according to OLS on full data set
  beta_E1 <- X_Ei1 %*% y1 # beta_hat according to OLS on selection data set
  z_E <- sign(beta[chosen])
  
  C <- solve(crossprod(X_E1, X_E1))
  if (intercept) {
    lamvec <- c(0, rep(lambda, s-1))
    b1 <- -diag(x = (z_E), nrow = s) %*% C %*% (lamvec * z_E) # linear constraint as derived in Lee et al. 2016
  } else {
    b1 <- -lambda * diag(x = (z_E), nrow = s) %*% C %*% z_E # linear constraint as derived in Lee et al. 2016 
  }
  if (n == splitn) {
    # pure post-selection inference
    # in this case we sample from c(beta_hat_selection)
    # linear constraint as derived in Lee et al. 2016,
    # adjusted for the fact, that sampling is not from y1, but from a linear transformation thereof
    linear_part <- matrix(0, nrow = s, ncol = s)
    linear_part[, 1:s] <- -diag(z_E)
    b <- b1
    
    # specify covariance of Gaussian vector c(beta_hat_selection)
    cov <- inv_info_E1
    con.covariance = cov * sigma**2
    
    # LHS of the additional equality constraint needed to sample for a specific coefficients
    # as explained in chapter 'Selective Inference for Linear Regression' of Fithian et al.
    conditional_linear = matrix(0,nrow = s, ncol = s)
    conditional_linear[, 1:s]  = ginv(inv_info_E1)
    
    OLS_func <- inv_info_E %*% conditional_linear
    # a valid initial condition
    initial <- c(beta_E1) 
  } else if (n - splitn > s) {
    # in this case we sample from c(beta_hat_full, beta_hat_selection)
    # linear constraint as derived in Lee et al. 2016,
    # adjusted for the fact, that sampling is not from y1, but from a linear transformation thereof
    linear_part <- matrix(0, nrow = s,ncol = 2 * s)
    linear_part[, (s + 1):(2 * s)] <- -diag(z_E)
    b <- b1 
    
    # specify covariance of 2s Gaussian vector
    # covariance of c(beta_hat_full, beta_hat_selection)
    cov <- matrix(0,2 * s, 2 * s)
    cov[1:s, 1:s] <- inv_info_E
    cov[(s + 1):(2 * s), 1:s] <- inv_info_E
    cov[1:s, (s + 1):(2 * s)] <- inv_info_E
    cov[(s + 1):(2 * s), (s + 1):(2 * s)] <- inv_info_E1
    con.covariance <- cov * sigma ** 2
    
    # for the conditional law
    # we will change the linear function for each coefficient
    selector <- matrix(0, nrow = s, ncol = 2 * s)
    selector[, 1:s] <- diag(s)
    # LHS of the additional equality constraint needed to sample for a specific coefficients
    # as explained in chapter 'Selective Inference for Linear Regression' of Fithian et al.
    conditional_linear <- crossprod(X_E, X_E) %*% selector
    
    # a valid initial condition (assuming glmnet fulfills KKT, not always the case)
    initial <- c(beta_E, beta_E1) 
    OLS_func <- selector 
  } else {
    # in this case we sample from c(beta_hat_selection, y2)
    # linear constraint as derived in Lee et al. 2016,
    # adjusted for the fact, that sampling is not from y1, but from a linear transformation thereof
    linear_part <- matrix(0, nrow = s, ncol = s + n - splitn)
    linear_part[, 1:s] <- -diag(z_E)
    b <- b1
    
    # specify covariance of Gaussian vector c(beta_hat_selection, y2)
    # beta_hat_selection and y2 are independent, so only diagonal blocks are non-empty
    cov <- matrix(0, nrow = (s + n - splitn), ncol = (s + n - splitn))
    cov[1:s, 1:s] <- inv_info_E1
    cov[(s + 1):(s + n - splitn), (s + 1):(s + n - splitn)] <- diag(n - splitn) 
    con.covariance <- cov * sigma ** 2
    # LHS of the additional equality constraint needed to sample for a specific coefficients
    # as explained in chapter 'Selective Inference for Linear Regression' of Fithian et al.
    conditional_linear <- matrix(0, nrow = s, ncol = (s + n - splitn))
    conditional_linear[, 1:s]  <- ginv(inv_info_E1)
    conditional_linear[,(s + 1):(s + n - splitn)] <- t(X[-ind, chosen])
    
    # write the OLS estimates of full model in terms of X_E1^{dagger}y_1, y2
    OLS_func = inv_info_E%*%conditional_linear
    
    # a valid initial condition (assuming glmnet fulfills KKT, not always the case)
    initial = c(beta_E1, y[-ind]) 
  }
  pvaluessum <- numeric(length(groups))
  j <- 0
  ndraw0 <- ndraw
  burnin0 <- burnin
  
  if (intercept) {
    ngrouptested <- sum(unlist(lapply(lapply(groups, intersect, chosen - 1), length)) > 0)
  } else {
    ngrouptested <- sum(unlist(lapply(lapply(groups, intersect, chosen), length)) > 0)
  }
  if (FWER) {
    sig_Level <- sig_Level / ngrouptested
  }
  correction_factor <- (1 - log(aggregation)) / aggregation
  min_size <- ceiling(1/(sig_Level/correction_factor))

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
      eta <- matrix(t(OLS_func[group.vars, ]), ncol = length(group.vars))
      # directions to project samples on to afterwards
      
      # new mean and covariance adjusting for additional equality constraint
      if (length(group.vars) < s) {
        conditional_law <- conditional(con.covariance, rep(0, dim(conditional_linear)[2]),
                                       matrix(conditional_linear[-group.vars, ], nrow = s -length(group.vars)),
                                       crossprod(X_E, y)[-group.vars])
        conditional_law.covariance <- con.covariance - conditional_law$delta_cov
        conditional_law.mean <- -conditional_law$delta_mean
      } else {
        # there is no additional equality constraint if all the selected variables are in the group
        conditional_law.covariance <- con.covariance
        conditional_law.mean <- rep(0,length(initial))
        
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
      
      white_out <- whiten(conditional_law.covariance, linear_part, b, conditional_law.mean)
      forward_map <- white_out$forward_map
      inverse_map <- white_out$inverse_map
      new_A <- white_out$new_A
      new_b <- white_out$new_b
      white_Y <- forward_map(initial)
      white_eta <- forward_map(conditional_law.covariance %*% eta)
      if (max (new_A %*% white_Y - new_b) > 0) stop("Constraints not fulfilled after whitening")
      
      # get a sample of points fulfilling all constraints
      Z <- sample_from_constraints(new_A, new_b, white_Y, white_eta, ndraw = ndraw,
                                   burnin = burnin, how_often = 10, verbose = verbose)
      Z <- t(inverse_map(t(Z)))
      
      etasign <- t(t(eta)*sign(beta[chosen[group.vars]]))
      continue <- FALSE

      i <- 0
      while (!continue) {
        i <- i + 1
        continue <- TRUE
        # project sampled points on to "directions of interest" to calculate p-value
        allnull_statistics <- Z %*% etasign
        allobserved <- initial %*% etasign
        if (i > 5) break()
        sumnull_statistics <- apply(allnull_statistics, 1, sum)
        sumobserved <- sum(allobserved)
        lennull <- length(sumnull_statistics)
        # calculate p-values on two halves of samples to estimate convergence
        pval1 <- (sum(sumnull_statistics[1:floor(lennull / 2)] >= sumobserved) + 1) / (floor(lennull / 2) + 1)
        pval2 <- (sum(sumnull_statistics[(floor(lennull / 2)+1):lennull] >= sumobserved) + 1) / (ceiling(lennull / 2) + 1)
        if (min(pval1, pval2) < sig_Level && ft) {
          # if potentially significant
          if (verbose) print("Checking significance with longer chain")
          ft <<- FALSE
          Z <- sample_from_constraints(new_A, new_b, forward_map(Z[lennull, ]), white_eta,
                                       ndraw = 2 * max(min_size, ndraw), burnin=0,
                                       how_often = 10, verbose = verbose)
          Z <- t(inverse_map(t(Z)))
          continue <- FALSE
        } else if (min(pval1, pval2) < sig_Level && max(pval1, pval2) > 1.5 * sig_Level) {
          # if significance is contradicting on the two parts
          if (verbose) print("Appending the chain")
          if (verbose) {
            print(paste("Adding", 2 * max(min_size, ndraw), "draws"))
          }
          ft <<- FALSE
          Zn <- sample_from_constraints(new_A, new_b, forward_map(Z[lennull, ]), white_eta,
                                        ndraw = 2 * max(min_size, ndraw), burnin = 0,
                                        how_often = 10, verbose = verbose)
          Zn <- t(inverse_map(t(Zn)))
          Z <- rbind(Z, Zn)
          continue <- FALSE
        }
      }
      
      
      # project sampled points on to "directions of interest" to calculate p-value
      allnull_statistics <- Z %*% etasign
      allobserved <- initial %*% etasign
      sumnull_statistics <- apply(allnull_statistics, 1, sum)
      sumobserved <- sum(allobserved)
      add <- 1 # add <- 0 would allow for 0 pvalues
      pvaluessum[j] <- (sum(sumnull_statistics >= sumobserved) + add) / (length(sumnull_statistics) + add)
    } else {
      if (verbose) print(paste("Not checking group", j))
      pvaluessum[j] <- 1
    }
    if (skip) nskipped <- nskipped + 1
    skip <<- FALSE
  }
  warning(paste("Hamiltonian sampler failed for", nskipped, "out of", ngrouptested, "groups"))
  return(list(pv = pvaluessum))
}


constraint_checker<-function(X1, y1, beta, tol.beta, lambda,
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

  
  X_E1 <- X1[, chosen] # active variables on selection data set 
  
  X_Ei1 <- ginv(X_E1)
  
  z_E <- sign(beta[chosen])
  C <- solve(crossprod(X_E1, X_E1))
  A1 <- -diag(x = (z_E), nrow = s) %*% C %*% t(X_E1)
  if (intercept) {
    lamvec <- c(0, rep(lambda, s - 1))
    b1 <- -diag(x = (z_E), nrow = s) %*% C %*% (lamvec * z_E) # linear constraint as derived in Lee et al. 2016
  } else {
    b1 <- -lambda * diag(x = (z_E), nrow = s) %*% C %*% z_E # linear constraint as derived in Lee et al. 2016 
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