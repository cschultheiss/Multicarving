oneSplit_select <- function(b) {
  sel.models <- logical(p)
  try.again <- TRUE
  fLI_error <- 0L
  split_count <- 0
  thresh_count <- 0L
  threshn <- 1e-7
  fLI_problem <- TRUE
  continue <- TRUE
  split_again <- TRUE
  while (split_again) {
    split_again <- FALSE
    split <- sample.int(n, size = n.left)
    x.left <- x[split, ]
    y.left <- y[split]
    x.right <- x[-split, ]
    y.right <- y[-split]
    
    output <- do.call(model.selector, args = c(list(x = x.left, 
                                                    y = y.left), args.model.selector))
    sel.model <- output$sel.model
    beta <- output$beta
    lambda <- output$lambda
    
    fit_again <- TRUE
    thresh_count <- 0
    p.sel <- length(sel.model)
    # for empty model, active constraints are trivially fulfilled
    if (p.sel == 0) fit_again <- FALSE
    
    while (fit_again) {
      fit_again <- FALSE
      checktry <- tryCatch_W_E(constraint_checker(x.left, y.left, beta, 0, lambda,
                                                  family, intercept = args.model.selector$intercept), TRUE)
      if (is.null(checktry$error)) {
        check <- checktry$value
      } else {
        # if run into numerical instability when checking constraints
        check <- TRUE
        split_again <- TRUE
        warning(paste(checktry$error, p.sel, " variables selected with ",
                      n.left, "data points, splitting again"))
      }
      if (!check) {
        if (verbose) 
          cat("......fit again...\n")
        fit_again <- TRUE
        thresh_count <- thresh_count + 1
        if (thresh_count > 2) {
          warning("Giving up reducing threshhold")
          break()
        }
        threshn <- 1e-7 / (100) ^ thresh_count
        fit <- glmnet(x = x.left, y = y.left, standardize = args.model.selector$standardize,
                      intercept = args.model.selector$intercept, thresh = threshn,family = family)
        if (verbose) cat(threshn,"\n")
        coefs <- coef(fit,x = x.left,y = y.left,s = lambda/n.left,exact = TRUE,
                      standardize = args.model.selector$standardize, 
                      intercept = args.model.selector$intercept, thresh = threshn, family = family)
        beta <- coefs[-1]
        sel.model <- which(abs(beta) > 0)
        
        if (family == "binomial") beta <- coefs
        p.sel <- length(sel.model)
        if (p.sel == 0) fit_again <- FALSE
        warning(paste("reducing threshold", thresh_count, "to", threshn, sep = " "))
      }
    }
    
    p.sel <- length(sel.model)
    # use new split in case of singularity. This is mostly an issue for discrete x.
    if (args.model.selector$intercept) {
      if ((p.sel > 0 && (rankMatrix(cbind(rep(1, n.left), x.left[, sel.model]))[[1]]< (p.sel + 1) ||
                         (p.sel < n.right - 1 && rankMatrix(cbind(rep(1, n.right), x.right[, sel.model]))[[1]] < (p.sel + 1)))) ||
          fit_again) split_again <- TRUE
    } else {
      if ((p.sel > 1 && (rankMatrix(x.left[, sel.model])[[1]] < (p.sel) ||
                         (p.sel < n.right  && rankMatrix(x.right[, sel.model])[[1]]< (p.sel)))) ||
          fit_again) split_again <- TRUE
    }
    if (split_again) {
      reason <- numeric(0)
      if (args.model.selector$intercept){
        if (rankMatrix(cbind(rep(1, n.left), x.left[,sel.model]))[[1]] < (p.sel + 1)) reason <- c(reason, "1")
        if (p.sel < n.right - 1 && rankMatrix(cbind(rep(1, n.right), x.right[,sel.model]))[[1]] < (p.sel + 1)) reason <- c(reason, "2")
      } else {
        if (rankMatrix( x.left[,sel.model])[[1]] < (p.sel)) reason <- c(reason, "1")
        if (p.sel < n.right && rankMatrix(x.right[,sel.model])[[1]] < (p.sel)) reason <- c(reason, "2")
      }
      
      if (fit_again) reason <- c(reason, "3")
      split_count <- split_count + 1
      if (split_count > 4) {
        stop(paste("More than 5 splits needed, final reason:", reason))
      }
      if (verbose) 
        cat("......Splitting again...\n")
      warning(paste("Splitting again ", split_count, "reason", reason))
    }
  }
  sel.models[sel.model] <- TRUE
  return (list(sel.models = sel.models, split = split, beta = beta, lambda = lambda))
}


oneSplit_infer<- function(sel) {
  if (split_pval) {
    pvals.v <- matrix(1, nrow = 2, ncol = p)
  } else {
    pvals.v <- rep(1, p)
  }
  sel.models <- sel$sel.models
  sel.model <- which(sel.models)
  p.sel <- length(sel.model)
  beta <- sel$beta
  split <- sel$split
  lambda <- sel$lambda
  
  if (se.estimator == "modwise" && family == "gaussian") {
    if (length(beta) == p + 1) beta <- beta[-1]
    if (args.model.selector$intercept){
      RSS <- sum((scale(y, T, F) - scale(x, T, F) %*% beta) ^ 2)
      if (args.se.estimator$df_corr) {
        den <- n - p.sel - 1
      } else {
        den <- n
      }
      sigma_model <- sqrt(RSS / den)
    } else {
      RSS <- sum((y- x %*% beta) ^ 2)
      if (args.se.estimator$df_corr) {
        den <- n - p.sel
      } else {
        den <- n
      }
      sigma_model <- sqrt(RSS / den)
    }
    estSigma <- sigma_model
    args.lasso.inference$sigma <- sigma_model
    args.classical.fit$Sigma <- sigma_model
  }
  if (p.sel > 0) {
    fLItry <- tryCatch_W_E(do.call(OptimalFixedLasso, args = c(list(X = x, y = y, ind = split, beta = beta, tol.beta = 0,
                                                                    lambda = lambda, intercept = args.model.selector$intercept),
                                                               args.lasso.inference)), 0)
    if (!is.null(fLItry$error)) {
      warning(paste("Failed to infer a split, due to:", fLItry$error, sep=" "))
      pvals.v[] = NA
      list(pvals = pvals.v, sel.models = sel.models, 
           split = split)
    } else if (!is.null(fLItry$warning)) {
      for (war in unique(fLItry$warning)) {
        warning(paste(war, sep = " "))
      }
    }
    fLI <- fLItry$value
    sel.pval1 <- fLI$pv
    if (any(is.na(sel.pval1))) {
      stop("The carve procedure returned a p-value NA")
    } 
    if (length(sel.pval1) != p.sel) { 
      stop(paste("The carve procedure didn't return the correct number of p-values for the provided submodel.",
                 p.sel, length(sel.pval1)))
    }
    if (!all(sel.pval1 >= 0 & sel.pval1 <= 1)) {
      stop("The carve procedure returned values below 0 or above 1 as p-values")
    }
    if (FWER) {
      sel.pval1 <- pmin(sel.pval1 * p.sel, 1) # for FWER
    } else {
      sel.pval1 <- pmin(sel.pval1, 1) # for FCR
    }
    if (split_pval) {
      x.right <- x[-split, ]
      y.right <- y[-split]
      if (args.model.selector$intercept) {
        bound <- n.right-1
      } else {
        bound <- n.right
      }
      if (p.sel < bound) {
        sel.pval2try <- tryCatch_W_E(do.call(classical.fit, 
                                             args = c(list(x = x.right[, sel.model], y = y.right), args.classical.fit)),
                                     rep(NA, p.sel))
        sel.pval2 <- sel.pval2try$value
        if (!is.null(sel.pval2try$error)) {
          warning(paste(sel.pval2try$error, "while caluclatng split p-values", sep=" "))
        }
        NAs <- FALSE
        if (any(is.na(sel.pval2))) NAs <- TRUE
        # do not stop if splitting leads to NA
        if (length(sel.pval2) != p.sel) 
          stop("The classical.fit function didn't return the correct number of p-values for the provided submodel.")
        if (!all(sel.pval2 >= 0 & sel.pval2 <= 1) && !NAs) 
          stop("The classical.fit function returned values below 0 or above 1 as p-values")
        if (FWER) {
          sel.pval2 <- pmin(sel.pval2 * p.sel, 1) # for FWER
        } else {
          sel.pval2 <- pmin(sel.pval2, 1) # for FCR
        }
        pvals.v[1, sel.model] <- sel.pval1
        pvals.v[2, sel.model] <- sel.pval2
      } else {
        # if split p-values can not be determined, leave them at 1
        pvals.v[1,sel.model] <- sel.pval1
      }
    } else {
      pvals.v[sel.model] <- sel.pval1
    }
    
  }
  
  if (p.sel == 0) {
    # leave all p-values to be 1
    if (verbose) 
      cat("......Empty model selected. That's ok...\n")
  }
  
  
  list(pvals = pvals.v, sel.models = sel.models, 
       split = split)
}

multi.carve.re <- function (x, y, B = 50, fraction = 0.9,
                         model.selector = lasso.cvcoef, classical.fit = lm.pval.flex,
                         parallel = FALSE, ncores = getOption("mc.cores", 2L), 
                         gamma = ((1:B)/B)[((1:B)/B) >= 0.05],
                         family = "gaussian",
                         args.model.selector = list(intercept = TRUE, standardize = FALSE),
                         se.estimator = "1se", args.se.estimator = list(df_corr = FALSE, intercept = TRUE, standardize = FALSE),
                         args.classical.fit = list(ttest = FALSE), return.nonaggr = FALSE, return.selmodels = FALSE,
                         verbose = FALSE, FWER = TRUE, split_pval= TRUE,
                         use_sigma_modwise = FALSE,
                         args.lasso.inference = list(sigma = NA)) {
  # routine to split the data, select a model and calculate carving p-values B times
  # x: matrix of predictors
  # y: response vector
  # B: number of splits
  # fraction: fraction used for selection
  # model.selector: how the model is chosen
  # classical.fit: function to calculate splitting p-values
  # parallel: whether to parallelize the splits
  # ncores: number of cores for parallelization
  # gamma: quantiles to consider, if several, additional penalty is applied
  # family: gaussian or binomial
  # args.model.selector: additional arguments for selection process
  # args.classical.fit: additional arguments for calculating splitting p-values
  # return.nonaggr: shall raw p-values be returned
  # return sel.models: shall the information, which model was selected be returned
  # verbose: whether to print key steps
  # FWER: shall a FWER correction be applied
  # split_pval: shall p-values for splitting be determined as well
  # use_sigma_modwise: shall sigma be calculated on a per model basis
  # args.lasso.inference: additional arguments for inference after Lasso
  
  
  args.model.selector$family <- family
  if (family == "gaussian"){
    if (se.estimator == "None" && is.na(args.lasso.inference$sigma)) stop("Neither SE estimator type nor sigma provided for Gaussian family. This is not ok")
    if (is.na(args.lasso.inference$sigma)) {
      if (se.estimator %in% c("1se", "modwise")){
        use_lambda.min = FALSE
      } else {
        use_lambda.min = TRUE
      }
      estSigma <- do.call(estimateSigma.flex,
                          args = c(list(x = x, y = y, use_lambda.min = use_lambda.min), args.se.estimator))
      globalSigma <- estSigma$sigmahat
      args.lasso.inference$sigma <- globalSigma
      args.classical.fit$Sigma <- globalSigma
    } else {
      # provided sigma has priority over se estimator
      se.estimator <- "None"
      globalSigma <- args.lasso.inference$sigma
      args.classical.fit$Sigma <- globalSigma
    }
  }
  
  n <- nrow(x)
  p <- ncol(x)
  n.left <- floor(n * fraction)
  n.right <- n - n.left
  stopifnot(n.left >= 1, n.right >= 0)
  oneSplit_select <- function(b) {
    sel.models <- logical(p)
    try.again <- TRUE
    fLI_error <- 0L
    split_count <- 0
    thresh_count <- 0L
    threshn <- 1e-7
    fLI_problem <- TRUE
    continue <- TRUE
    split_again <- TRUE
    while (split_again) {
      split_again <- FALSE
      split <- sample.int(n, size = n.left)
      x.left <- x[split, ]
      y.left <- y[split]
      x.right <- x[-split, ]
      y.right <- y[-split]
      
      output <- do.call(model.selector, args = c(list(x = x.left, 
                                                      y = y.left), args.model.selector))
      sel.model <- output$sel.model
      beta <- output$beta
      lambda <- output$lambda
      
      fit_again <- TRUE
      thresh_count <- 0
      p.sel <- length(sel.model)
      # for empty model, active constraints are trivially fulfilled
      if (p.sel == 0) fit_again <- FALSE
      
      while (fit_again) {
        fit_again <- FALSE
        checktry <- tryCatch_W_E(constraint_checker(x.left, y.left, beta, 0, lambda,
                                                    family, intercept = args.model.selector$intercept), TRUE)
        if (is.null(checktry$error)) {
          check <- checktry$value
        } else {
          # if run into numerical instability when checking constraints
          check <- TRUE
          split_again <- TRUE
          warning(paste(checktry$error, p.sel, " variables selected with ",
                        n.left, "data points, splitting again"))
        }
        if (!check) {
          if (verbose) 
            cat("......fit again...\n")
          fit_again <- TRUE
          thresh_count <- thresh_count + 1
          if (thresh_count > 2) {
            warning("Giving up reducing threshhold")
            break()
          }
          threshn <- 1e-7 / (100) ^ thresh_count
          fit <- glmnet(x = x.left, y = y.left, standardize = args.model.selector$standardize,
                        intercept = args.model.selector$intercept, thresh = threshn,family = family)
          if (verbose) cat(threshn,"\n")
          coefs <- coef(fit,x = x.left,y = y.left,s = lambda/n.left,exact = TRUE,
                        standardize = args.model.selector$standardize, 
                        intercept = args.model.selector$intercept, thresh = threshn, family = family)
          beta <- coefs[-1]
          sel.model <- which(abs(beta) > 0)
          
          if (family == "binomial") beta <- coefs
          p.sel <- length(sel.model)
          if (p.sel == 0) fit_again <- FALSE
          warning(paste("reducing threshold", thresh_count, "to", threshn, sep = " "))
        }
      }
      
      p.sel <- length(sel.model)
      # use new split in case of singularity. This is mostly an issue for discrete x.
      if (args.model.selector$intercept) {
        if ((p.sel > 0 && (rankMatrix(cbind(rep(1, n.left), x.left[, sel.model]))[[1]]< (p.sel + 1) ||
                           (p.sel < n.right - 1 && rankMatrix(cbind(rep(1, n.right), x.right[, sel.model]))[[1]] < (p.sel + 1)))) ||
            fit_again) split_again <- TRUE
      } else {
        if ((p.sel > 1 && (rankMatrix(x.left[, sel.model])[[1]] < (p.sel) ||
                           (p.sel < n.right  && rankMatrix(x.right[, sel.model])[[1]]< (p.sel)))) ||
            fit_again) split_again <- TRUE
      }
      if (split_again) {
        reason <- numeric(0)
        if (args.model.selector$intercept){
          if (rankMatrix(cbind(rep(1, n.left), x.left[,sel.model]))[[1]] < (p.sel + 1)) reason <- c(reason, "1")
          if (p.sel < n.right - 1 && rankMatrix(cbind(rep(1, n.right), x.right[,sel.model]))[[1]] < (p.sel + 1)) reason <- c(reason, "2")
        } else {
          if (rankMatrix( x.left[,sel.model])[[1]] < (p.sel)) reason <- c(reason, "1")
          if (p.sel < n.right && rankMatrix(x.right[,sel.model])[[1]] < (p.sel)) reason <- c(reason, "2")
        }
        
        if (fit_again) reason <- c(reason, "3")
        split_count <- split_count + 1
        if (split_count > 4) {
          stop(paste("More than 5 splits needed, final reason:", reason))
        }
        if (verbose) 
          cat("......Splitting again...\n")
        warning(paste("Splitting again ", split_count, "reason", reason))
      }
    }
    sel.models[sel.model] <- TRUE
    return (list(sel.models = sel.models, split = split, beta = beta, lambda = lambda))
  }
  oneSplit_infer<- function(sel) {
    if (split_pval) {
      pvals.v <- matrix(1, nrow = 2, ncol = p)
    } else {
      pvals.v <- rep(1, p)
    }
    sel.models <- sel$sel.models
    sel.model <- which(sel.models)
    p.sel <- length(sel.model)
    beta <- sel$beta
    split <- sel$split
    lambda <- sel$lambda
    
    if (se.estimator == "modwise" && family == "gaussian") {
      if (length(beta) == p + 1) beta <- beta[-1]
      if (args.model.selector$intercept){
        RSS <- sum((scale(y, T, F) - scale(x, T, F) %*% beta) ^ 2)
        if (args.se.estimator$df_corr) {
          den <- n - p.sel - 1
        } else {
          den <- n
        }
        sigma_model <- sqrt(RSS / den)
      } else {
        RSS <- sum((y- x %*% beta) ^ 2)
        if (args.se.estimator$df_corr) {
          den <- n - p.sel
        } else {
          den <- n
        }
        sigma_model <- sqrt(RSS / den)
      }
      estSigma <- sigma_model
      args.lasso.inference$sigma <- sigma_model
      args.classical.fit$Sigma <- sigma_model
    }
    if (p.sel > 0) {
      fLItry <- tryCatch_W_E(do.call(OptimalFixedLasso, args = c(list(X = x, y = y, ind = split, beta = beta, tol.beta = 0,
                                                                      lambda = lambda, intercept = args.model.selector$intercept),
                                                                 args.lasso.inference)), 0)
      if (!is.null(fLItry$error)) {
        warning(paste("Failed to infer a split, due to:", fLItry$error, sep=" "))
        pvals.v[] = NA
        list(pvals = pvals.v, sel.models = sel.models, 
             split = split)
      } else if (!is.null(fLItry$warning)) {
        for (war in unique(fLItry$warning)) {
          warning(paste(war, sep = " "))
        }
      }
      fLI <- fLItry$value
      sel.pval1 <- fLI$pv
      if (any(is.na(sel.pval1))) {
        stop("The carve procedure returned a p-value NA")
      } 
      if (length(sel.pval1) != p.sel) { 
        stop(paste("The carve procedure didn't return the correct number of p-values for the provided submodel.",
                   p.sel, length(sel.pval1)))
      }
      if (!all(sel.pval1 >= 0 & sel.pval1 <= 1)) {
        stop("The carve procedure returned values below 0 or above 1 as p-values")
      }
      if (FWER) {
        sel.pval1 <- pmin(sel.pval1 * p.sel, 1) # for FWER
      } else {
        sel.pval1 <- pmin(sel.pval1, 1) # for FCR
      }
      if (split_pval) {
        x.right <- x[-split, ]
        y.right <- y[-split]
        if (args.model.selector$intercept) {
          bound <- n.right-1
        } else {
          bound <- n.right
        }
        if (p.sel < bound) {
          sel.pval2try <- tryCatch_W_E(do.call(classical.fit, 
                                               args = c(list(x = x.right[, sel.model], y = y.right), args.classical.fit)),
                                       rep(NA, p.sel))
          sel.pval2 <- sel.pval2try$value
          if (!is.null(sel.pval2try$error)) {
            warning(paste(sel.pval2try$error, "while caluclatng split p-values", sep=" "))
          }
          NAs <- FALSE
          if (any(is.na(sel.pval2))) NAs <- TRUE
          # do not stop if splitting leads to NA
          if (length(sel.pval2) != p.sel) 
            stop("The classical.fit function didn't return the correct number of p-values for the provided submodel.")
          if (!all(sel.pval2 >= 0 & sel.pval2 <= 1) && !NAs) 
            stop("The classical.fit function returned values below 0 or above 1 as p-values")
          if (FWER) {
            sel.pval2 <- pmin(sel.pval2 * p.sel, 1) # for FWER
          } else {
            sel.pval2 <- pmin(sel.pval2, 1) # for FCR
          }
          pvals.v[1, sel.model] <- sel.pval1
          pvals.v[2, sel.model] <- sel.pval2
        } else {
          # if split p-values can not be determined, leave them at 1
          pvals.v[1,sel.model] <- sel.pval1
        }
      } else {
        pvals.v[sel.model] <- sel.pval1
      }
      
    }
    
    if (p.sel == 0) {
      # leave all p-values to be 1
      if (verbose) 
        cat("......Empty model selected. That's ok...\n")
    }
    
    
    list(pvals = pvals.v, sel.models = sel.models, 
         split = split)
  }
  split.out <- if (parallel) {
    stopifnot(isTRUE(is.finite(ncores)), ncores >= 1L)
    if (verbose) 
      cat("...starting parallelization of sample-splits\n")
    mclapply(1:B, oneSplit, mc.cores = ncores)
  } else {
    lapply(1:B, function(b) {oneSplit_infer(oneSplit_select())})

  }
  myExtract <- function(name) {
    matrix(unlist(lapply(split.out, "[[", name)), nrow = B, 
           byrow = TRUE)
  }
  if (split_pval) { 
    ls <- list()
    pvalsall <- array(unlist(lapply(split.out, "[[", "pvals")), dim = c(2, p, B))
    for (icf in  1:2) {
      pvals <- t(pvalsall[icf, , ])
      colnames(pvals) <- colnames(x)
      if (return.selmodels) {
        sel.models <- myExtract("sel.models")
        colnames(sel.models) <- colnames(x)
      } else {
        sel.models <- NA
      }
      pvals.current <- which.gamma <- numeric(p)
      for (j in 1:p) {
        if (any(!is.na(pvals[, j]))) {
          quant.gamma <- quantile(pvals[, j], gamma, na.rm = TRUE, type = 1) / gamma
          penalty <- if (length(gamma) > 1) 
            (1 - log(min(gamma)))
          else 1
          pvals.pre <- min(quant.gamma) * penalty
          pvals.current[j] <- pmin(pvals.pre, 1)
          which.gamma[j] <- which.min(quant.gamma)
        } else {
          pvals.current[j] <- NA
          which.gamma[j] <- NA
        }
      }
      names(pvals.current) <- names(which.gamma) <- colnames(x)
      
      if (!return.nonaggr) 
        pvals <- NA
      if (return.selmodels) {
        if (icf==2) {
          keep <- c("return.selmodels", "x", "y", "gamma", "split.out", 
                    "pvals", "pvals.current", "which.gamma", "sel.models","ls","icf")
          rm(list = setdiff(names(environment()), keep))
        }
      }
      ls[[icf]] <- structure(list(pval = NA, pval.corr = pvals.current, pvals.nonaggr = pvals, 
                                  gamma.min = gamma[which.gamma], sel.models = sel.models,
                                  method = "multi.split", call = match.call()), class = "hdi")
      
    }
    return(ls)
  } else {
    pvals <- myExtract("pvals")
    colnames(pvals) <- colnames(x)
    if (return.selmodels) {
      sel.models <- myExtract("sel.models")
      colnames(sel.models) <- colnames(x)
    } else {
      sel.models <- NA
    }
    pvals.current <- which.gamma <- numeric(p)
    for (j in 1:p) {
      quant.gamma <- quantile(pvals[, j], gamma) / gamma
      penalty <- if (length(gamma) > 1) 
        (1 - log(min(gamma)))
      else 1
      pvals.pre <- min(quant.gamma) * penalty
      pvals.current[j] <- pmin(pvals.pre, 1)
      which.gamma[j] <- which.min(quant.gamma)
    }
    names(pvals.current) <- names(which.gamma) <- colnames(x)
    if (!return.nonaggr) 
      pvals <- NA
    if (return.selmodels) {
      keep <- c("return.selmodels", "x", "y", "gamma", "split.out", 
                "pvals", "pvals.current", "which.gamma", "sel.models")
      rm(list = setdiff(names(environment()), keep))
    }
    structure(list(pval = NA, pval.corr = pvals.current, pvals.nonaggr = pvals, 
                   gamma.min = gamma[which.gamma], sel.models = sel.models,
                   method = "multi.split", call = match.call()), class = "hdi")
  }
}