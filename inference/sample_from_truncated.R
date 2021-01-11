# this file contains the hit-and-run sampler applicable to carving

sample.truncnorm.white <- function( A, b, initial, bias.direction, how.often = 1000,
                                   sigma = 1, burnin = 500, ndraw = 1000, use.A = FALSE) {
  # Sample from a truncated normal with covariance equal to sigma**2 I.
  # Constraint is $Ax \leq b$ where `A` has shape
  # `(q,n)` with `q` the number of constraints and
  # `n` the number of random variables.
  # Parameters
  # A : Linear part of affine constraints.
  # b : Offset part of affine constraints.
  # initial : Initial point for Gibbs draws assumed to satisfy the constraints.
  # bias.direction : Which projections are of most interest?
  # how.often : How often should the sampler make a move along `bias.direction`?
  # If negative, defaults to ndraw+burnin  (so it will never be used).
  # sigma : Variance parameter, usually 1 since pre-withened
  # burnin : How many iterations until we start recording samples?
  # ndraw : How many samples should we return?
  # use.A : if true, some of the movement is in the constrained directions 
  dims <- dim(A)
  nvar <- dims[2]
  nconstraint <- dims[1]
  trunc.sample <- matrix(NA, nrow = ndraw, ncol = nvar)
  state <- initial 
  tol <- 1.e-7
  
  U <- A %*% state - b
  # max(U) < 0 if KKT actually fulfilled exactly
  if (max(U) > 0) {
    fc <- sum(U > 0)
    stop(paste(fc, "out of", nconstraint, "unfulfilled constraints at entry"))
  }
  neta <- dim(bias.direction)[2] # 1 for single variable testing
  # choose random number a priori
  usample <- runif(burnin + ndraw)
  
  # directions not parallel to coordinate axes
  if (use.A) {
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
  
  directions <- rbind(directions, t(bias.direction)) # last directions are directions of interest
  directions <- t(directions)
  # make columns have unit length since only direction is of interest
  directions <- scale(directions, center = FALSE, scale = sqrt(colSums(directions ^ 2))) 
  directions <-t(directions)
  ndir <- dim(directions)[1]
  
  alphas.dir <- A %*% t(directions)
  alphas.coord <- A
  alphas.max.dir <- apply(abs(alphas.dir), 2 ,max) * tol
  alphas.max.coord <- apply(abs(alphas.coord), 2, max) * tol
  
  
  # choose the order of sampling (randomly)
  random.idx.dir <- sample(1:ndir, burnin + ndraw, replace = TRUE)
  random.idx.coord <- sample(1:nvar, burnin + ndraw, replace = TRUE)

  # for switching between coordinate updates and
  # other directions
  invperiod <- 13 # when to do an update not in main coordinate direction
  docoord <- 0
  iperiod <- 0
  ibias <- 0
  dobias <- 0
  iter.count <- 0
  rand.count <- 0

  while (iter.count  < (ndraw + burnin)) {
    rand.count <- (rand.count) %% (ndraw + burnin) + 1
    docoord <- 1
    iperiod <- iperiod + 1
    ibias <- ibias + 1
  
    if (iperiod == invperiod) {
      # do direction update
      docoord <- 0
      iperiod <- 0
      dobias <- 0
    }
  
    if (ibias == how.often) {
      # do direction of interest update
      docoord <- 0
      ibias <- 0
      dobias <- 1
    }
  
    if (docoord == 1) {
      idx <- random.idx.coord[rand.count]
      V <- state[idx]
      # V is given coordinate of state, i.e. projection on to coordinate axis
    } else if (!dobias) {
      idx <- random.idx.dir[rand.count]
      V <- sum(directions[idx, ] * state)
      # V is projection of state in to direction
    } else {
      idx <- dim(directions)[1] - sample(0:(neta - 1), 1)
      V <- sum(directions[idx, ] * state)
      # last rows of directions are bias.directions
    }
      
    lower.bound <- -1e12
    upper.bound <- 1e12
    # bounds to walk in the given direction, as derived in Lee et al. 2016
    for (irow in 1:nconstraint) {
      if (docoord == 1) {
        alpha <- alphas.coord[irow, idx]
        val <- -U[irow] / alpha + V
        if (alpha > alphas.max.coord[idx] && (val < upper.bound)) {
          upper.bound <- val
        } else if (alpha < (-alphas.max.coord[idx]) && (val > lower.bound)) {
          lower.bound <- val
        }
      } else {
        alpha <- alphas.dir[irow, idx]
        val <- -U[irow] / alpha + V
        if (alpha > alphas.max.dir[idx] && (val < upper.bound)) {
          upper.bound <- val
        } else if (alpha < (-alphas.max.dir[idx]) && (val > lower.bound)) {
          lower.bound <- val
        }
      }
    } 
  
    if (lower.bound > V) {
      lower.bound <- V - tol * sigma
    } else if (upper.bound < V) {
      upper.bound <- V + tol * sigma
    }
  
  
    lower.bound <- lower.bound / sigma
    upper.bound <- upper.bound / sigma
    # 
  
    if (lower.bound > upper.bound) {
      if (max(A %*% initial - b) > 0) {
        stop("boundary error with initial error")
      } else {
        stop("boundary error without initial error") 
      }
      # this should only happen if KKT not exactly fulfilled
    }
  
    if (upper.bound < (-20)) {
      # avoid numerical instability in extreme cases
      unif <- usample[rand.count] * (1 - exp(-abs(lower.bound - upper.bound)
                                             * abs(upper.bound)))
      tnorm <- (upper.bound + log(1 - unif) / abs(upper.bound)) * sigma
    } else if (lower.bound > 20) {
      # avoid numerical instability in extreme cases
      unif <- usample[rand.count] * (1 - exp(-abs(upper.bound - lower.bound)
                                             * lower.bound))
      tnorm <- (lower.bound - log(1 - unif) / lower.bound) * sigma
    } else {
      # standard update
      if (upper.bound < 0){
        # both < 0
        cdfL <- pnorm(lower.bound)
        cdfU <- pnorm(upper.bound)
        unif <- usample[rand.count] * (cdfU - cdfL) + cdfL
        tnorm <- qnorm(unif) * sigma
      } else if(lower.bound > 0) {
        # both > 0
        cdfL <- pnorm(-lower.bound)
        cdfU <- pnorm(-upper.bound)
        unif <- usample[rand.count] * (cdfU - cdfL) + cdfL
        tnorm <- -qnorm(unif) * sigma
      } else {
        # upper.bound >0, lower.bound <0
        cdfU <- 1 - pnorm(-upper.bound)
        cdfL <- pnorm(lower.bound)
        unif <- usample[rand.count] * (cdfU - cdfL) + cdfL
        if (unif < 0.5){
          tnorm <- qnorm(unif) * sigma
        } else {
          tnorm <- -qnorm(1 - unif) * sigma
        }
        
      }
    }
    if (tnorm > upper.bound || tnorm < lower.bound) stop("new value exceeds bounds")

    iter.count <- iter.count + 1
    if (docoord == 1) {
      state[idx] <- tnorm # update that coordinate
      tnorm <- tnorm - V
      U <- U + tnorm * A[, idx] # update U = A * state - b
    } else {
      tnorm <- tnorm - V
      state <- state + tnorm * directions[idx, ] # update state in that direction
      U <- U + tnorm*A%*%directions[idx,] # update U = A * state - b
    }
  
    if (iter.count >= burnin) {
      # add sample after burning period
      trunc.sample[iter.count - burnin, ] <- state
    }
  }
  return (trunc.sample)
}









