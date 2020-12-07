# this file contains the hit-and-run sampler applicable to carving

sample_truncnorm_white <- function( A, b, initial, bias_direction, how_often = 1000,
                                   sigma = 1, burnin = 500, ndraw = 1000, use_A = FALSE) {
  # Sample from a truncated normal with covariance equal to sigma**2 I.
  # Constraint is $Ax \leq b$ where `A` has shape
  # `(q,n)` with `q` the number of constraints and
  # `n` the number of random variables.
  # Parameters
  # A : Linear part of affine constraints.
  # b : Offset part of affine constraints.
  # initial : Initial point for Gibbs draws assumed to satisfy the constraints.
  # bias_direction : Which projections are of most interest?
  # how_often : How often should the sampler make a move along `direction_of_interest`?
  # If negative, defaults to ndraw+burnin  (so it will never be used).
  # sigma : Variance parameter, usually 1 since pre-withened
  # burnin : How many iterations until we start recording samples?
  # ndraw : How many samples should we return?
  # use_A : if true, some of the movement is in the constrained directions 
  
  dims <- dim(A)
  nvar <- dims[2]
  nconstraint <- dims[1]
  trunc_sample <- matrix(NA, nrow = ndraw, ncol = nvar)
  state <- initial 
  tol <- 1.e-7
  
  U <- A %*% state - b
  # max(U) < 0 if KKT actually fulfilled exactly
  if (max(U) > 0) {
    fc <- sum(U > 0) / nconstraint
    stop(paste(fc, "unfulfilled constraints at entry"))
  }
  neta <- dim(bias_direction)[2] # 1 for single variable testing
  # choose random number a priori
  usample <- runif(burnin + ndraw)
  
  # directions not parallel to coordinate axes
  if (use_A) {
    # walk directions are partly from A, partly random
    if (nvar >= 5) {
      directions  <- rbind(A, matrix(rnorm(n = as.integer(nvar / 5) * nvar),
                                     nrow = as.integer(nvar / 5), ncol = nvar))
    } else {
      directions <- A
    }
  } else {
    # walk directions are only random
    directions <- matrix(rnorm(nvar ^ 2), nrow = nvar, ncol = nvar)
  }
  
  directions <- rbind(directions, t(bias_direction)) # last directions are directions of interest
  directions <- t(directions)
  # make columns have unit length since only direction is of interest
  directions <- scale(directions, center = FALSE, scale = sqrt(colSums(directions ^ 2))) 
  directions <-t(directions)
  ndir <- dim(directions)[1]
  
  alphas_dir <- A %*% t(directions)
  alphas_coord <- A
  alphas_max_dir <- apply(abs(alphas_dir), 2 ,max) * tol
  alphas_max_coord <- apply(abs(alphas_coord), 2, max) * tol
  
  
  # choose the order of sampling (randomly)
  random_idx_dir <- sample(1:ndir, burnin + ndraw, replace = TRUE)
  random_idx_coord <- sample(1:nvar, burnin + ndraw, replace = TRUE)

  # for switching between coordinate updates and
  # other directions
  invperiod <- 13 # when to do an update not in main coordinate direction
  docoord <- 0
  iperiod <- 0
  ibias <- 0
  dobias <- 0
  iter_count <- 0
  rand_count <- 0

  while (iter_count  < (ndraw + burnin)) {
    rand_count <- (rand_count) %% (ndraw + burnin) + 1
    docoord <- 1
    iperiod <- iperiod + 1
    ibias <- ibias + 1
  
    if (iperiod == invperiod) {
      # do direction update
      docoord <- 0
      iperiod <- 0
      dobias <- 0
    }
  
    if (ibias == how_often) {
      # do direction of interest update
      docoord <- 0
      ibias <- 0
      dobias <- 1
    }
  
    if (docoord == 1) {
      idx <- random_idx_coord[rand_count]
      V <- state[idx]
      # V is given coordinate of state, i.e. projection on to coordinate axis
    } else if (!dobias) {
      idx <- random_idx_dir[rand_count]
      V <- sum(directions[idx, ] * state)
      # V is projection of state in to direction
    } else {
      idx <- dim(directions)[1] - sample(0:(neta - 1), 1)
      V <- sum(directions[idx, ] * state)
      # last rows of directions are bias_directions
    }
      
    lower_bound <- -1e12
    upper_bound <- 1e12
    # bounds to walk in the given direction, as derived in Lee et al. 2016
    for (irow in 1:nconstraint) {
      if (docoord == 1) {
        alpha <- alphas_coord[irow, idx]
        val <- -U[irow] / alpha + V
        if (alpha > alphas_max_coord[idx] && (val < upper_bound)) {
          upper_bound <- val
        } else if (alpha < (-alphas_max_coord[idx]) && (val > lower_bound)) {
          lower_bound <- val
        }
      } else {
        alpha <- alphas_dir[irow, idx]
        val <- -U[irow] / alpha + V
        if (alpha > alphas_max_dir[idx] && (val < upper_bound)) {
          upper_bound <- val
        } else if (alpha < (-alphas_max_dir[idx]) && (val > lower_bound)) {
          lower_bound <- val
        }
      }
    } 
  
    if (lower_bound > V) {
      lower_bound <- V - tol * sigma
    } else if (upper_bound < V) {
      upper_bound <- V + tol * sigma
    }
  
  
    lower_bound <- lower_bound / sigma
    upper_bound <- upper_bound / sigma
    # 
  
    if (lower_bound > upper_bound) {
      if (max(A %*% initial - b) > 0) {
        stop("boundary error with initial error")
      } else {
        stop("boundary error without initial error") 
      }
      # this should only happen if KKT not exactly fulfilled
    }
  
    if (upper_bound < (-20)) {
      # avoid numerical instability in extreme cases
      unif <- usample[rand_count] * (1 - exp(-abs(lower_bound - upper_bound)
                                             * abs(upper_bound)))
      tnorm <- (upper_bound + log(1 - unif) / abs(upper_bound)) * sigma
    } else if (lower_bound > 20) {
      # avoid numerical instability in extreme cases
      unif <- usample[rand_count] * (1 - exp(-abs(upper_bound - lower_bound)
                                             * lower_bound))
      tnorm <- (lower_bound - log(1 - unif) / lower_bound) * sigma
    } else {
      # standard update
      if (upper_bound < 0){
        # both < 0
        cdfL <- pnorm(lower_bound)
        cdfU <- pnorm(upper_bound)
        unif <- usample[rand_count] * (cdfU - cdfL) + cdfL
        tnorm <- qnorm(unif) * sigma
      } else if(lower_bound > 0) {
        # both > 0
        cdfL <- pnorm(-lower_bound)
        cdfU <- pnorm(-upper_bound)
        unif <- usample[rand_count] * (cdfU - cdfL) + cdfL
        tnorm <- -qnorm(unif) * sigma
      } else {
        # upper_bound >0, lower_bound <0
        cdfU <- 1 - pnorm(-upper_bound)
        cdfL <- pnorm(lower_bound)
        unif <- usample[rand_count] * (cdfU - cdfL) + cdfL
        if (unif < 0.5){
          tnorm <- qnorm(unif) * sigma
        } else {
          tnorm <- -qnorm(1 - unif) * sigma
        }
        
      }
    }
    if (tnorm > upper_bound || tnorm < lower_bound) stop("new value exceeded bounds")

    iter_count <- iter_count + 1
    if (docoord == 1) {
      state[idx] <- tnorm # update that coordinate
      tnorm <- tnorm - V
      U <- U + tnorm * A[, idx] # update U = A * state - b
    } else {
      tnorm <- tnorm - V
      state <- state + tnorm * directions[idx, ] # update state in that direction
      U <- U + tnorm*A%*%directions[idx,] # update U = A * state - b
    }
  
    if (iter_count >= burnin) {
      # add sample after burning period
      trunc_sample[iter_count - burnin, ] <- state
    }
  }
  return (trunc_sample)
}









