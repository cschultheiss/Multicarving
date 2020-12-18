# this file contains functions adapted from the functions in hdi for multisplitting
# making them applicable to multicarving

multi.carve <- function (x, y, B = 50, fraction = 0.9,
                            model.selector = lasso.cvcoef, classical.fit = lm.pval.flex,
                            parallel = FALSE, ncores = getOption("mc.cores", 2L), 
                            gamma = ((1:B)/B)[((1:B)/B) >= 0.05],
                            family = "gaussian",
                            args.model.selector = list(intercept = TRUE, standardize = FALSE),
                            se.estimator = "1se", args.se.estimator = list(df.corr = FALSE, intercept = TRUE, standardize = FALSE),
                            args.classical.fit = list(ttest = FALSE), return.nonaggr = FALSE, return.selmodels = FALSE, skip.variables = TRUE,
                            verbose = FALSE, FWER = TRUE, split.pval= TRUE,
                            args.lasso.inference = list(sigma = NA, sig_Level = 0.05, FWER = FWER, aggregation = min(gamma))) {
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
  # skip.variables: shall carving p-values for variables selected less than min(gamma)*B times be omitted
  # verbose: whether to print key steps
  # FWER: shall a FWER correction be applied
  # split.pval: shall p-values for splitting be determined as well
  # args.lasso.inference: additional arguments for inference after Lasso
  
  
  args.model.selector$family <- family
  args.lasso.inference$family <- family
  
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
  sel.out <- if (parallel) {
    stopifnot(isTRUE(is.finite(ncores)), ncores >= 1L)
    if (verbose) 
      cat("...starting parallelization of sample-splits\n")
    mclapply(1:B, oneSplit_select, mc.cores = ncores)
  } else {
    lapply(1:B, oneSplit_select)
  }
  if (skip.variables){
    myExtract.sel <- function(name) {
      matrix(unlist(lapply(sel.out, "[[", name)), nrow = B, 
             byrow = TRUE)
    }
    sel.models <- myExtract.sel("sel.models")
    times.selected <- apply(sel.models, 2, sum)
    which.check <- which(times.selected >= min(gamma) *B)
    warning(paste("Reducing number of tests from", sum(sel.models), "to", sum(sel.models[, which.check])))
  } else {
    which.check <- NULL
  }
  
  oneSplit_infer<- function(sel) {
    if (split.pval) {
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
        if (args.se.estimator$df.corr) {
          den <- n - p.sel - 1
        } else {
          den <- n
        }
        sigma_model <- sqrt(RSS / den)
      } else {
        RSS <- sum((y- x %*% beta) ^ 2)
        if (args.se.estimator$df.corr) {
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
                                                                      lambda = lambda, intercept = args.model.selector$intercept, which.check = which.check),
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
      if (split.pval) {
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
  inf.out <- if (parallel) {
    stopifnot(isTRUE(is.finite(ncores)), ncores >= 1L)
    if (verbose) 
      cat("...starting parallelization of sample-splits\n")
    mclapply(sel.out, oneSplit_infer, mc.cores = ncores)
  } else {
    lapply(sel.out, oneSplit_infer)
  }
  myExtract <- function(name) {
    matrix(unlist(lapply(inf.out, "[[", name)), nrow = B, 
           byrow = TRUE)
  }
  if (split.pval) { 
    ls <- list()
    pvalsall <- array(unlist(lapply(inf.out, "[[", "pvals")), dim = c(2, p, B))
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
          quant.gamma <- quantile(pvals[, j], gamma, na.rm = TRUE, type = 3) / gamma
          quant.gamma1 <- quantile(pvals[, j], gamma, na.rm = TRUE, type = 1) / gamma
          if (min(quant.gamma) < min (quant.gamma1)) warning(paste("Could reduce p-value for variable", j,  "with other quantile function"))
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
          keep <- c("return.selmodels", "x", "y", "gamma", "inf.out", 
                    "pvals", "pvals.current", "which.gamma", "sel.models", "FWER","ls","icf")
          rm(list = setdiff(names(environment()), keep))
        }
      }
      method = "multi.carve"
      if (icf == 2) method = "multi.split"
      ls[[icf]] <- structure(list(pval.corr = pvals.current, pvals.nonaggr = pvals, 
                                  gamma.min = gamma[which.gamma], sel.models = sel.models, FWER = FWER,
                                  method = method, call = match.call()), class = "carve")
      
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
      quant.gamma <- quantile(pvals[, j], gamma, na.rm = TRUE, type = 3) / gamma
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
      keep <- c("return.selmodels", "x", "y", "gamma", "inf.out", 
                "pvals", "pvals.current", "which.gamma", "sel.models", "FWER")
      rm(list = setdiff(names(environment()), keep))
    }
    structure(list(pval.corr = pvals.current, pvals.nonaggr = pvals, 
                   gamma.min = gamma[which.gamma], sel.models = sel.models, FWER = FWER,
                   method = "multi.carve", call = match.call()), class = "carve")
  }
}

carve100 <- function (x, y, model.selector = lasso.cvcoef, family = "gaussian", args.model.selector = list(intercept = TRUE, standardize = FALSE, tol.beta = 1e-5),
                      return.selmodels = FALSE, verbose = FALSE, FWER = TRUE, estimate.sigma = TRUE, df.corr = FALSE,
                      args.lasso.inference = list(sigma = NA)) {

  args.model.selector$family <- family
  args.lasso.inference$family <- family
  
  if (family == "gaussian" && !estimate.sigma && is.na(args.lasso.inference$sigma)) stop("Sigma not provided and estimation not enabled for Gaussian family. This is not ok")
  if (!is.na(args.lasso.inference$sigma)) estimate.sigma <- FALSE

  n <- nrow(x)
  p <- ncol(x)
  pvals.v <- rep(1,p)
  sel.models <- logical(p)

  thresh_count <- 0L
  threshn <- 1e-7
  split_again <- FALSE

  output <- do.call(model.selector, 
                    args = c(list(x = x, y = y), args.model.selector))
  sel.model <- output$sel.model
  beta <- output$beta
  lambda <- output$lambda

  abort <- FALSE
  fit_again <- TRUE
  thresh_count <- 0
  p.sel <- length(sel.model)
  if (p.sel == 0) fit_again <- FALSE
  
  while (fit_again) {
    fit_again <- FALSE
    checktry <- tryCatch_W_E(constraint_checker(x, y, beta, 0, lambda, family,
                                                intercept = args.model.selector$intercept), TRUE)
    if (is.null(checktry$error)) {
      check <- checktry$value
    } else {
      check <- TRUE
      abort <- TRUE
      warning(paste(checktry$error, p.sel, " variables selected with ",
                    n, "data points, bad selection"))
      break()
    }
    if (!check) {
      fit_again <- TRUE
      thresh_count <- thresh_count + 1
      if (thresh_count > 4) {
        warning("Giving up reducing threshhold")
        abort <- TRUE
        break()
      }
      threshn <- 1e-7 / (100) ^ thresh_count
      fit <- glmnet(x = x, y = y, standardize = args.model.selector$standardize,
                    intercept = args.model.selector$intercept, thresh = threshn,family = family)
      cat(threshn, "\n")
      coefs <- coef(fit, x = x,y = y, s = lambda / n, exact = TRUE, standardize = args.model.selector$standardize,
                    intercept = args.model.selector$intercept, thresh = threshn, family = family)
      beta <- coefs[-1]
      
      if (args.model.selector$intercept) {
        sel.model <- which(abs(beta) > args.model.selector$tol.beta * sqrt(nrow(x) / colSums(scale(x, T, F) ^ 2))) # model indices
      } else {
        sel.model <- which(abs(beta) > args.model.selector$tol.beta * sqrt(nrow(x) / colSums(x ^ 2))) # model indices
      }
      
      if (family == "binomial") beta <- coefs
      p.sel <- length(sel.model)
      if (p.sel == 0) fit_again <- FALSE
      warning(paste("reducing threshold", thresh_count, "from", threshn))
    }
  }
  if (abort) {
    stop("Can not fit that data")
  }
  
  p.sel <- length(sel.model)
  if (args.model.selector$intercept) {
    if ((p.sel > 1 && (rankMatrix(cbind(rep(1, n), x[, sel.model]))[[1]] < (p.sel + 1))) ||
        fit_again) abort <- TRUE
  } else {
    if ((p.sel > 1 && (rankMatrix(x[, sel.model])[[1]] < (p.sel))) ||
        fit_again) abort <- TRUE
  }
  if (abort) {
    stop("Can not fit that data")
  }
  
  if (estimate.sigma && family == "gaussian") {
    if (length(beta) == p + 1) beta <- beta[-1]
    if (args.model.selector$intercept){
      RSS <- sum((scale(y, T, F) - scale(x, T, F) %*% beta) ^ 2)
      if (df.corr) {
        den <- n - p.sel - 1
      } else {
        den <- n
      }
      sigma_model <- sqrt(RSS / den)
    } else {
      RSS <- sum((y- x %*% beta) ^ 2)
      if (df.corr) {
        den <- n - p.sel
      } else {
        den <- n
      }
      sigma_model <- sqrt(RSS / den)
    }
    args.lasso.inference$sigma <- sigma_model
  }
  if (p.sel > 0) {
    if (length(beta) == p+1 && family == "gaussian") beta <- beta[-1]
    fLItry <- tryCatch_W_E(do.call(fixedLassoInf, args = c(list(x = x, y = y, beta = beta, 
                                                                tol.beta = args.model.selector$tol.beta,
                                                                lambda = lambda, intercept = args.model.selector$intercept),
                                                           args.lasso.inference)), 0)
    if (!is.null(fLItry$error)) {
      stop(paste(fLItry$error, "stopping carve100", sep=" "))
    } else if (!is.null(fLItry$warning)) {
      war <- unique(fLItry$warning)
      stop(paste(war, "stopping carve100", sep = " "))
    }
    fLI <- fLItry$value
    if (p.sel > 0 && p.sel < nrow(x) - 1) {
      sel.pval1 <- fLI$pv
      if (any(is.na(sel.pval1))) 
        stop("The carve 100 procedure returned a p-value NA")
      if (length(sel.pval1) != p.sel) { 
        stop(paste("The carve 100 procedure didn't return the correct number of p-values for the provided submodel.",
                   p.sel, length(sel.pval1)))
      }
      if (!all(sel.pval1 >= 0 & sel.pval1 <= 1)) 
        stop("The carve 100 procedure returned values below 0 or above 1 as p-values")
      if (FWER) {
        sel.pval1 <- pmin(sel.pval1 * p.sel, 1) #for FWER
      } else {
        sel.pval1 <- pmin(sel.pval1, 1) #for FCR
      }
      pvals.v[sel.model] <- sel.pval1
      
      if (return.selmodels) 
        sel.models[sel.model] <- TRUE
      }
  }
  
  if (p.sel == 0) {
    if (verbose) 
      cat("......Empty model selected. That's ok...\n")
  }
  pvals <- pvals.v
  colnames(pvals) <- colnames(x)
  if (return.selmodels) {
    colnames(sel.models) <- colnames(x)
  } else {
    sel.models <- NA
  }
  
  if (return.selmodels) {
    keep <- c("return.selmodels", "x", "y", 
              "pvals", "sel.models", "FWER")
    rm(list = setdiff(names(environment()), keep))
  }
  structure(list(pval.corr = pvals, FWER= FWER,
                 sel.models = sel.models, method = "carve100", call = match.call()), class = "carve")
}

multi.carve.group <- function (x, y, groups, B = 50, fraction = 0.9, family = "gaussian",
                               model.selector = lasso.cvcoef, parallel = FALSE, ncores = getOption("mc.cores", 2L),
                               gamma = ((1:B)/B)[((1:B)/B) >= 0.05], args.model.selector = list(intercept = TRUE, standardize = FALSE),
                               return.nonaggr = FALSE, return.selmodels = FALSE, skip.groups = TRUE,
                               se.estimator = "1se", args.se.estimator = list(df.corr = FALSE, intercept = TRUE, standardize = FALSE),
                               verbose = FALSE, FWER = FALSE,
                               args.lasso.inference = list(sigma = NA, sig_Level = 0.05, FWER = FWER, aggregation = min(gamma))) {
  args.lasso.inference$family <- family
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
    } else {
      # provided sigma has priority over se estimator
      se.estimator <- "None"
      globalSigma <- args.lasso.inference$sigma
    }
  }
  
  n <- nrow(x)
  p <- ncol(x)
  ngroup <- length(groups)
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
  sel.out <- if (parallel) {
    stopifnot(isTRUE(is.finite(ncores)), ncores >= 1L)
    if (verbose) 
      cat("...starting parallelization of sample-splits\n")
    mclapply(1:B, oneSplit_select, mc.cores = ncores)
  } else {
    lapply(1:B, oneSplit_select)
  }
  if (skip.groups){
    myExtract.sel <- function(name) {
      matrix(unlist(lapply(sel.out, "[[", name)), nrow = B, 
             byrow = TRUE)
    }
    sel.models <- myExtract.sel("sel.models")
    sel.groups <- matrix(FALSE, ncol = length(groups), nrow = B)
    for (i in 1:B) {
      j <- 0
      for (group in groups){
        j <- j + 1
        sel.groups[i, j] <- sum(sel.models[i, group]) > 0
      }
    }
    times.selected <- apply(sel.groups, 2, sum)
    which.check <- which(times.selected >= min(gamma) * B)
    warning(paste("Reducing number of tests from", sum(sel.groups), "to", sum(sel.groups[, which.check])))
  } else {
    which.check <- NULL
  }
  
  oneSplit_infer<- function(sel) {
    pvals.v <- rep(1, ngroup)
    sel.models <- sel$sel.models
    sel.model <- which(sel.models)
    p.sel <- length(sel.model)
    beta <- sel$beta
    split <- sel$split
    lambda <- sel$lambda
    ngrouptested <- sum(unlist(lapply(lapply(groups, intersect, sel.model), length)) > 0)
    if (se.estimator == "modwise" && family == "gaussian") {
      if (length(beta) == p + 1) beta <- beta[-1]
      if (args.model.selector$intercept){
        RSS <- sum((scale(y, T, F) - scale(x, T, F) %*% beta) ^ 2)
        if (args.se.estimator$df.corr) {
          den <- n - p.sel - 1
        } else {
          den <- n
        }
        sigma_model <- sqrt(RSS / den)
      } else {
        RSS <- sum((y- x %*% beta) ^ 2)
        if (args.se.estimator$df.corr) {
          den <- n - p.sel
        } else {
          den <- n
        }
        sigma_model <- sqrt(RSS / den)
      }
      estSigma <- sigma_model
      args.lasso.inference$sigma <- sigma_model
    }
    if (ngrouptested > 0) {
      fLItry <- tryCatch_W_E(do.call(OptimalFixedLassoGroup, args = c(list(X = x, y = y, ind = split, beta = beta, tol.beta = 0,
                                                                      lambda = lambda, intercept = args.model.selector$intercept, groups = groups, which.check = which.check),
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
      if (length(sel.pval1) != length(groups)) { 
        stop(paste("The carve procedure didn't return the correct number of p-values for the provided groups.",
                   length(groups), length(sel.pval1)))
      }
      if (!all(sel.pval1 >= 0 & sel.pval1 <= 1)) {
        stop("The carve procedure returned values below 0 or above 1 as p-values")
      }
      if (FWER) {
        sel.pval1 <- pmin(sel.pval1 * ngrouptested, 1) #for FWER
      } else {
        sel.pval1 <- pmin(sel.pval1, 1) #for FCR
      }
      pvals.v <- sel.pval1

      
    }
    
    if (ngrouptested == 0) {
      # leave all p-values to be 1
      if (verbose) 
        cat("......No variable from the groups selected. That's ok...\n")
    }
    
    
    list(pvals = pvals.v, sel.models = sel.models, 
         split = split)
  }
  inf.out <- if (parallel) {
    stopifnot(isTRUE(is.finite(ncores)), ncores >= 1L)
    if (verbose) 
      cat("...starting parallelization of sample-splits\n")
    mclapply(sel.out, oneSplit_infer, mc.cores = ncores)
  } else {
    lapply(sel.out, oneSplit_infer)
  }
  myExtract <- function(name) {
    matrix(unlist(lapply(inf.out, "[[", name)), nrow = B, 
           byrow = TRUE)
  }
  pvals <- myExtract("pvals")
  colnames(pvals) <- colnames(x)
  if (return.selmodels) {
    sel.models <- myExtract("sel.models")
    colnames(sel.models) <- colnames(x)
  }
  else {
    sel.models <- NA
  }
  pvals.current <- which.gamma <- numeric(ngroup)
  for (j in 1:ngroup) {
    quant.gamma <- quantile(pvals[, j], gamma, na.rm = TRUE, type = 3) / gamma
    penalty <- if (length(gamma) > 1) 
      (1 - log(min(gamma)))
    else 1
    pvals.pre <- min(quant.gamma) * penalty
    pvals.current[j] <- pmin(pvals.pre, 1)
    which.gamma[j] <- which.min(quant.gamma)
  }
  names(pvals.current) <- names(which.gamma) <- colnames(paste(groups))
  if (!return.nonaggr) 
    pvals <- NA
  if (return.selmodels) {
    keep <- c("return.selmodels", "x", "y", "gamma", "inf.out", 
              "pvals", "pvals.current", "which.gamma", "sel.models", "FWER")
    rm(list = setdiff(names(environment()), keep))
  }
  structure(list(pval.corr = pvals.current, pvals.nonaggr = pvals, 
                 gamma.min = gamma[which.gamma], FWER = FWER,
                 sel.models = sel.models, method = "multi.carve.group", call = match.call()), class = "carve")
}

multi.carve.ci.saturated <- function(x, y, B = 50, fraction = 0.9, ci.level = 0.95, model.selector = lasso.cvcoef,
                                     classical.fit = lm.pval, parallel = FALSE, ncores = getOption("mc.cores", 2L),
                                     gamma = ((1:B)/B)[((1:B)/B) >= 0.05], family = "gaussian",
                                     args.model.selector = list(intercept = TRUE, standardize = FALSE),
                                     se.estimator = "modwise", args.se.estimator = list(df.corr = TRUE, intercept = TRUE, standardize = FALSE),
                                     args.classical.fit = NULL, args.classical.ci = NULL, return.nonaggr = FALSE, 
                                     return.selmodels = FALSE, verbose = FALSE, ci.timeout = 10,
                                     FWER = FALSE, split.pval = TRUE, ttest = TRUE,
                                     args.lasso.inference = list(sigma = NA)) {

  args.model.selector$family <- family
  args.lasso.inference$family <- family
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
    } else {
      # provided sigma has priority over se estimator
      se.estimator <- "None"
      globalSigma <- args.lasso.inference$sigma
    }
  }
  
  n <- nrow(x)
  p <- ncol(x)
  n.left <- floor(n * fraction)
  n.right <- n - n.left
  stopifnot(n.left >= 1, n.right >= 0)
  oneSplit <- function(b) {
    if (split.pval) {
      pvals.v <- matrix(1, nrow = 2, ncol = p)
    } else {
      pvals.v <- rep(1, p)
    }
    sel.models <- logical(p)
    vlo.v <- rep(-Inf, p)
    vup.v <- rep(Inf, p)
    estimates.v <- rep(NA, p)
    sescarve.v <- rep(Inf, p)
    lci.v <- rep(-Inf, p)
    uci.v <- rep(Inf, p)
    centers.v <- rep(NA, p)
    ses.v <- rep(Inf, p)
    df.res <- NA
    try.again <- TRUE
    thresh_count <- 0L
    threshn <- 1e-7
    continue <- TRUE
    split_again <- TRUE
    split_count <- 0
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
      if (p.sel == 0) fit_again <- FALSE
      
      while (fit_again) {
        fit_again <- FALSE
        checktry <- tryCatch_W_E(constraint_checker(x.left, y.left, beta, 0, lambda, family,
                                                    intercept = args.model.selector$intercept), TRUE)
        if (is.null(checktry$error)) {
          check <- checktry$value
        } else {
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
          fit <- glmnet(x = x.left,y = y.left,standardize = args.model.selector$standardize,
                        intercept = args.model.selector$intercept, thresh = threshn, family = family)
          if (verbose) cat(threshn,"\n")
          coefs <- coef(fit, x = x.left, y = y.left, s = lambda / n.left, exact = TRUE,
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

    if (se.estimator == "modwise" && family == "gaussian") {
      if (length(beta) == p + 1) beta <- beta[-1]
      if (args.model.selector$intercept){
        RSS <- sum((scale(y, T, F) - scale(x, T, F) %*% beta) ^ 2)
        if (args.se.estimator$df.corr) {
          den <- n - p.sel - 1
        } else {
          den <- n
        }
        sigma_model <- sqrt(RSS / den)
      } else {
        RSS <- sum((y- x %*% beta) ^ 2)
        if (args.se.estimator$df.corr) {
          den <- n - p.sel
        } else {
          den <- n
        }
        sigma_model <- sqrt(RSS / den)
      }
      estSigma <- sigma_model
      args.lasso.inference$sigma <- sigma_model
    }
      
    if (p.sel > 0) {
      fLItry <- tryCatch_W_E(do.call(OptimalFixedLasso,
                                     args = c(list(X = x, y = y,ind = split, beta = beta, tol.beta = 0,
                                                   lambda = lambda, intercept = args.model.selector$intercept,
                                                   selected = FALSE), args.lasso.inference)), 0)
      if (!is.null(fLItry$error)) {
        warning(paste("Failed to infer a split, due to:", fLItry$error, sep=" "))
        pvals.v[] = NA
        list(pvals = pvals.v, sel.models = sel.models, centers = centers.v, 
             ses = ses.v, df.res = df.res, lci = lci.v, uci = uci.v, sescarve = sescarve.v,
             vlo = vlo.v, vup = vup.v, estimates= estimates.v, split = split)
      } else if (!is.null(fLItry$warning)) {
        for (war in unique(fLItry$warning)) {
          warning(paste(war, sep = " "))
        }
      }
      fLI <- fLItry$value
      sel.pval1 <- fLI$pv
      sel.vlo <- fLI$vlo
      sel.vup <- fLI$vup
      sel.sescarve <- fLI$ses
      sel.estimates <- fLI$estimates
      
      if (any(is.na(sel.pval1))) {
        warning(sel.pval1)
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
      if (split.pval) {
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
          tmp.fit.lm <- lm(y.right ~ x.right[, sel.model], 
                           args.classical.fit)
          a <- (1 - ci.level)/2
          a <- c(a, 1 - a)
          fac <- qt(a, tmp.fit.lm$df.residual)
          sel.ses <- sqrt(diag(vcov(tmp.fit.lm)))[-1]
          sel.centers <- coef(tmp.fit.lm)[-1]
          sel.ci <- sel.centers + sel.ses %o% fac
          centers.v[sel.model] <- sel.centers
          lci.v[sel.model] <- sel.ci[, 1]
          uci.v[sel.model] <- sel.ci[, 2]
          ses.v[sel.model] <- sel.ses
          df.res <- tmp.fit.lm$df.residual
          pvals.v[1, sel.model] <- sel.pval1
          pvals.v[2, sel.model] <- sel.pval2
          vlo.v[sel.model] <- sel.vlo
          vup.v[sel.model] <- sel.vup
          estimates.v[sel.model] <- sel.estimates
          sescarve.v[sel.model] <- sel.sescarve
        } else {
          pvals.v[1, sel.model] <- sel.pval1
          vlo.v[sel.model] <- sel.vlo
          vup.v[sel.model] <- sel.vup
          estimates.v[sel.model] <- sel.estimates
          sescarve.v[sel.model] <- sel.sescarve
        }
      } else {
        pvals.v[sel.model] <- sel.pval1
        vlo.v[sel.model] <- sel.vlo
        vup.v[sel.model] <- sel.vup
        estimates.v[sel.model] <- sel.estimates
        sescarve.v[sel.model] <- sel.sescarve
      }
          
      if (return.selmodels) 
        sel.models[sel.model] <- TRUE
    }
    
    if (p.sel == 0) {
      if (verbose) 
        cat("......Empty model selected. That's ok...\n")
    }
    list(pvals = pvals.v, sel.models = sel.models, centers = centers.v, 
         ses = ses.v, df.res = df.res, lci = lci.v, uci = uci.v, sescarve = sescarve.v,
         vlo = vlo.v, vup = vup.v, estimates= estimates.v, split = split)
  }
  split.out <- if (parallel) {
    stopifnot(isTRUE(is.finite(ncores)), ncores >= 1L)
    if (verbose) 
      cat("...starting parallelization of sample-splits\n")
    mclapply(1:B, oneSplit, mc.cores = ncores)
  } else {
    if (verbose) 
      lapply(1:B, function(b) {
        cat("...Split", b, "\n")
        oneSplit()
      })
    else replicate(B, oneSplit(), simplify = FALSE)
    
  }
  myExtract <- function(name) {
    matrix(unlist(lapply(split.out, "[[", name)), nrow = B, 
           byrow = TRUE)
  }
  if (split.pval) { 
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
      if (icf == 1) {
        vlo <- myExtract("vlo")
        vars <- ncol(vlo)
        vup <- myExtract("vup")
        sescarve <- myExtract("sescarve")
        estimates <- myExtract("estimates")
      } else {
        lci <- myExtract("lci")
        uci <- myExtract("uci")
        centers <- myExtract("centers")
        ses <- myExtract("ses")
      }
      
      df.res <- unlist(lapply(split.out, `[[`, "df.res"))
      pvals.current <- which.gamma <- numeric(p)
      for (j in 1:p) {
        if (any(!is.na(pvals[, j]))) {
          quant.gamma <- quantile(pvals[, j], gamma, na.rm = TRUE, type = 3) / gamma
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
      s0 <- if (any(is.na(sel.models)))
        NA
      else apply(sel.models, 1, sum)
      if (icf == 1) {
        new.ci <- mapply(aggregate.ci.saturated, vlo = split(vlo, rep(1:vars, each = B)),
                         vup = split(vup, rep(1:vars, each = B)), 
                         centers = split(estimates, rep(1:vars, each = B)), 
                         ses = split(sescarve, rep(1:ncol(sescarve), each = B)), 
                         gamma.min = min(gamma), multi.corr = FWER, verbose = FALSE, timeout = ci.timeout,
                         s0 = list(s0 = s0), ci.level = ci.level, var = 1:vars)
      } else {
        new.ci <- mapply(hdi:::aggregate.ci, lci = split(lci, rep(1:vars, each = B)),
                         rci = split(uci, rep(1:vars, each = B)),
                         centers = split(centers, rep(1:vars, each = B)),
                         ses = split(ses, rep(1:ncol(ses), each = B)), df.res = list(df.res = df.res),
                         gamma.min = min(gamma), multi.corr = FALSE, verbose = FALSE,
                         s0 = list(s0 = s0), ci.level = ci.level, var = 1:vars)
      }
      lci.current <- t(new.ci)[, 1]
      uci.current <- t(new.ci)[, 2]
      names(lci.current) <- names(uci.current) <- names(pvals.current)
      
      if (!return.nonaggr)
        pvals <- NA
      if (return.selmodels) {
        if (icf == 2) {
          keep <- c("return.selmodels", "x", "y", "gamma", "split.out",
                    "pvals", "pvals.current", "which.gamma", "sel.models",
                    "ls", "icf", "ci.level", "estimates", "FWER",
                    "lci.current", "uci.current", "vlo", "vup", "sescarve", "ses")
          rm(list = setdiff(names(environment()), keep))
        }
      }
      if (icf == 1) {
        ls[[icf]] <- structure(list(pval.corr = pvals.current, pvals.nonaggr = pvals, 
                                    ci.level = ci.level, lci = lci.current, vlo = vlo, vup = vup, ses = sescarve, centers = estimates,
                                    uci = uci.current, gamma.min = gamma[which.gamma], FWER = FWER,
                                    sel.models = sel.models, method = "multi.carve", call = match.call()), class = "carve")
      } else {
        ls[[icf]] <- structure(list(pval.corr = pvals.current, pvals.nonaggr = pvals,
                                    ci.level = ci.level, lci =  lci.current, ses=ses,
                                    uci = uci.current, gamma.min = gamma[which.gamma], FWER = FWER,
                                    sel.models = sel.models, method = "multi.split", call = match.call()), class = "carve")
      }
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
    vlo <- myExtract("vlo")
    vup <- myExtract("vup")
    estimates <- myExtract("estimates")
    sescarve <- myExtract("sescarve")
    pvals.current <- which.gamma <- numeric(p)
    for (j in 1:p) {
      quant.gamma <- quantile(pvals[, j], gamma, na.rm = TRUE,type = 3) / gamma
      penalty <- if (length(gamma) > 1) 
        (1 - log(min(gamma)))
      else 1
      pvals.pre <- min(quant.gamma) * penalty
      pvals.current[j] <- pmin(pvals.pre, 1)
      which.gamma[j] <- which.min(quant.gamma)
    }
    names(pvals.current) <- names(which.gamma) <- colnames(x)
    vars <- ncol(vlo)
    s0 <- if (any(is.na(sel.models))) 
      NA
    else apply(sel.models, 1, sum)
    new.ci <- mapply(aggregate.ci.saturated, vlo = split(vlo, rep(1:vars, each = B)),
                     vup = split(vup, rep(1:vars, each = B)), 
                     centers = split(estimates, rep(1:vars, each = B)), 
                     ses = split(sescarve, rep(1:ncol(sescarve), each = B)), 
                     gamma.min = min(gamma), multi.corr = FALSE, verbose = FALSE, timeout = ci.timeout,
                     s0 = list(s0 = s0), ci.level = ci.level, var = 1:vars)
    lci.current <- t(new.ci)[, 1]
    uci.current <- t(new.ci)[, 2]
    if (!return.nonaggr) 
      pvals <- NA
    names(lci.current) <- names(uci.current) <- names(pvals.current)
    if (return.selmodels) {
      keep <- c("return.selmodels", "x", "y", "gamma", "split.out", "FWER",
                "pvals", "pvals.current", "which.gamma", "sel.models", "ci.level", 
                "lci.current", "uci.current", "vlo", "vup", "sescarve", "estimates")
      rm(list = setdiff(names(environment()), keep))
    }
    structure(list(pval.corr = pvals.current, pvals.nonaggr = pvals, 
                   ci.level = ci.level, lci = lci.current, vlo = vlo, vup = vup, ses = sescarve, centers = estimates,
                   uci = uci.current, gamma.min = gamma[which.gamma], FWER = FWER,
                   sel.models = sel.models, method = "multi.carve", call = match.call()), class = "carve")
  }
}



pval.aggregator <- function(pval.list, gamma, cutoff = TRUE) {
  aggregated.list <- list()
  i <- 0
  for (pvals in pval.list) {
    i <- i+1
    p <- dim(pvals)[2]
    pvals.current <- numeric(p)
    for (j in 1:p) {
      quant.gamma <- quantile(pvals[, j], gamma, na.rm = TRUE, type = 3) / gamma
      penalty <- if (length(gamma) > 1) 
        (1 - log(min(gamma)))
      else 1
      pvals.pre <- min(quant.gamma) * penalty
      if (cutoff) pvals.pre <- pmin(pvals.pre, 1)
      pvals.current[j] <- pvals.pre
    }
    aggregated.list[[i]] <- pvals.current
  }
  return(aggregated.list)
}

# different Lasso selector

lasso.firstqcoef <- function (x, y, q, tol.beta = 0, return_intercept = NULL, ...) {
  fit <- glmnet(x, y, dfmax = q,...)
  m <- predict(fit, type = "nonzero")
  delta <- q - unlist(lapply(m, length))
  delta[delta < 0] <- Inf
  take <- which.min(delta)
  sel.model <- m[[take]]
  lambda <- fit$lambda[take]
  coefs <- coef(fit, s = lambda)
  beta <- coefs[2 : (dim(x)[2] + 1)]
  if (is.null(return_intercept)) return_intercept <- (abs(coefs[1]) > 0)
  if (return_intercept) {
    return_beta <- coefs[1 : (dim(x)[2] + 1)]
  } else {
    return_beta <- beta
  }
  chosen <- which(abs(beta) > tol.beta * sqrt(nrow(x) / colSums(x ^ 2))) # model indices
  return(list(sel.model = chosen, beta = return_beta, lambda = lambda * dim(x)[1]))
}

lasso.cvcoef<-function (x, y, nfolds = 10, grouped = nrow(x) > 3 * nfolds,
                        tol.beta = 0, use_lambda.min = FALSE, return_intercept = NULL, ...) {
  fit.cv <- cv.glmnet(x, y, nfolds = nfolds, grouped = grouped, thresh = 1e-7, ...)
  if (use_lambda.min) {
    lambda <- fit.cv$lambda.min
  } else {
    lambda <- fit.cv$lambda.1se
  }
  coefs <- coef(fit.cv,x=x,y=y,s=lambda,exact=TRUE)
  beta <- coefs[2:(dim(x)[2]+1)]
  if (is.null(return_intercept)) return_intercept <- (abs(coefs[1])>0)
  if (return_intercept) {
    return_beta <- coefs[1 : (dim(x)[2] + 1)]
  } else {
    return_beta <- beta
  }
  if (coefs[1] != 0) {
    chosen <- which(abs(beta) > tol.beta * sqrt(nrow(x) / colSums(scale(x, T, F) ^ 2 ))) # model indices
  } else {
    chosen <- which(abs(beta) > tol.beta * sqrt(nrow(x) / colSums(x ^ 2))) # model indices
  }
  return(list(sel.model = chosen,beta = return_beta,lambda = lambda * dim(x)[1]))
}

fixedLasso.modelselector<-function(x, y, lambda, tol.beta, thresh = 1e-7, exact = FALSE, ...) {
  fit<-glmnet(x, y, alpha = 1, thresh = thresh,...)
  if (exact) {
    coefs <-coef(fit, s = lambda / (dim(x)[1]),x = x,
                 y = y, exact = exact, thresh = thresh) # lambda/n1 due to different definition of the LASSO loss function 
  } else {
    coefs <- coef(fit, s = lambda / (dim(x)[1])) # lambda/n1 due to different definition of the LASSO loss function
  }
  
  beta <- coefs[2 : (dim(x)[2] + 1)]
  if (coefs[1] != 0) {
    chosen <- which(abs(beta) > tol.beta * sqrt(nrow(x) / colSums(scale(x, T, F) ^ 2 ))) # model indices
  } else {
    chosen <- which(abs(beta) > tol.beta * sqrt(nrow(x) / colSums(x ^ 2))) # model indices
  }
  return(list(sel.model = chosen, beta = beta, lambda = lambda))
}

estimateSigma.flex <- function (x, y, intercept = TRUE, standardize = FALSE, use_lambda.min = FALSE, df.corr = FALSE) {
  selectiveInference:::checkargs.xy(x, rep(0, nrow(x)))
  if (nrow(x) < 10) 
    stop("Number of observations must be at least 10 to run estimateSigma")
  cvfit <- cv.glmnet(x, y, intercept = intercept, standardize = standardize)
  lamhat <- cvfit$lambda.1se
  if (use_lambda.min) lamhat <- cvfit$lambda.min
  fit <- glmnet(x, y, intercept = intercept, standardize = standardize)
  yhat <- predict(fit, x, s = lamhat)
  nz <- sum(predict(fit, s = lamhat, type = "coef") != 
             0)
  den <- length(y)
  if (df.corr) {
    den <- den - nz
  }

  sigma = sqrt(sum((y - yhat) ^ 2) / den)
  return(list(sigmahat = sigma, df = nz))
}

lm.pval.flex <- function (x, y, exact = TRUE, intercept = TRUE, Sigma = NA, ttest = TRUE, ...) {
  if (intercept) {
    fit.lm <- lm(y ~ x, ...)
    fit.summary <- summary(fit.lm)
    tstat <- coef(fit.summary)[-1, "t value"]
  } else {
    fit.lm <- lm(y ~ -1 + x,...) 
    fit.summary <- summary(fit.lm)
    tstat <- coef(fit.summary)[, "t value"]
  }
  
  if (is.na(Sigma) || ttest) {
    setNames(2 * (if (exact) 
      pt(abs(tstat), df = fit.lm$df.residual, lower.tail = FALSE)
      else pnorm(abs(tstat), lower.tail = FALSE)), colnames(x)) 
  } else {
    sigma_hat<-sqrt(sum((fit.lm$residuals) ^ 2) / fit.lm$df.residual)
    setNames(2 * pnorm(abs(tstat * sigma_hat / Sigma), lower.tail = FALSE), colnames(x))
  }
}

glm.pval.pseudo<-function(x, y, maxit = 100, delta = 0.01, epsilon = 1e-06) {
  increase <- TRUE
  incs <- 0
  while (increase) {
    increase <- FALSE
    pi.hat <- max(delta, min(1 - delta, mean(y)))
    delta.0 <- (pi.hat * delta)/(1 + delta)
    delta.1 <- (1 + pi.hat * delta)/(1 + delta)
    y.tilde <- delta.0 * (1 - y) + delta.1 * y
    pseudo.y <- cbind(y.tilde, 1 - y.tilde)
    s <- dim(x)[2]
    str<-paste("pseudo.y~", paste("x[,", 1:s, "]", collapse = "+"), collapse = "")
    tryfit <- tryCatch_W_E(glm(formula = as.formula(str), family = "binomial", maxit = maxit, epsilon = epsilon), 0)
    if ("glm.fit: fitted probabilities numerically 0 or 1 occurred" %in% tryfit$warning) {
      delta <- 1.1*delta
      incs <- incs+1
      increase <- TRUE
      if (delta > 0.1) break("Increased delta too much without succes")
    } else {
      fit <- tryfit$value
    }
  }
  if (incs>0) warning(paste("Increased delta to ", delta))
  return(glm.pval(x, pseudo.y, maxit = maxit, epsilon = epsilon))
}



aggregate.ci.saturated <- function(vlo, vup, centers, ses, gamma.min, multi.corr = FALSE,
                                   verbose = FALSE, var, s0, ci.level, timeout = 10) {
  inf.ci <- is.infinite(vlo) & is.infinite(vup)
  no.inf.ci <- sum(inf.ci)
  B <- length(vlo)
  if (verbose) {
    cat("number of Inf ci:", no.inf.ci, "\n")
  }
  if ((no.inf.ci == B) || (no.inf.ci > (1 - gamma.min) * 
                                     B)) {
    # variable not select frequently enough
    return(c(-Inf, Inf))
  }
  vlo <- vlo[!inf.ci]
  vup <- vup[!inf.ci]
  centers <- centers[!inf.ci]
  ses <- ses[!inf.ci]
  s0 <- s0[!inf.ci]
  ci.info <- list(vlo = vlo, vup = vup, centers = centers, 
                  no.inf.ci = no.inf.ci, ses = ses, 
                  s0 = s0,  gamma.min = gamma.min, multi.corr = multi.corr, 
                  ci.level = ci.level)
  # find a starting point within the interval
  inner <- find.inside.point.gammamin.saturated(low = min(centers), high = max(centers), 
                                      ci.info = ci.info, verbose = verbose)

  outer <- min(vlo[is.finite(vlo)])
  newboundstry <- tryCatch_W_E(find.bisection.bounds.gammamin.saturated(shouldcover = inner, shouldnotcover = outer,
                                                                        ci.info = ci.info, verbose = verbose, timeout = timeout), NA)
  if (!is.null(newboundstry$error)) {
    warning(paste(newboundstry$error, " setting lower interval limit to -Inf for variable ", var))
    l.bound  <- -Inf
  } else {
    new.bounds <- newboundstry$value
    inner <- new.bounds$shouldcover
    outer <- new.bounds$shouldnotcover
    
    l.bound <- bisection.gammamin.coverage.saturated(outer = outer, inner = inner, 
                                                     ci.info = ci.info, verbose = verbose)
  }

  if (verbose) {
    cat("lower bound ci aggregated is", l.bound, "\n")
  }
  outer <- max(vup[is.finite(vup)])
  newboundstry <- tryCatch_W_E(find.bisection.bounds.gammamin.saturated(shouldcover = inner, shouldnotcover = outer, 
                                                                        ci.info = ci.info, verbose = verbose, timeout=timeout), NA)
  if (!is.null(newboundstry$error)) {
    warning(paste(newboundstry$error, " setting upper interval limit to Inf for variable ", var))
    u.bound  <- Inf
  } else {
    new.bounds <- newboundstry$value
    inner <- new.bounds$shouldcover
    outer <- new.bounds$shouldnotcover
    u.bound <- bisection.gammamin.coverage.saturated(inner = inner, outer = outer, 
                                                     ci.info = ci.info, verbose = verbose)
  }

  if (verbose) {
    cat("upper bound ci aggregated is", u.bound, "\n")
  }
  return(c(l.bound, u.bound))
}

find.inside.point.gammamin.saturated <- function (low, high, ci.info, verbose) {
  range.length <- 10
  test.range <- seq(low, high, length.out = range.length)
  cover <- mapply(does.it.cover.gammamin.saturated, beta.j = test.range, 
                  ci.info = list(ci.info = ci.info))
  while (!any(cover)) {
    range.length <- 10 * range.length
    test.range <- seq(low, high, length.out = range.length)
    cover <- mapply(does.it.cover.gammamin.saturated, beta.j = test.range, 
                    ci.info = list(ci.info = ci.info))
    if (range.length > 10^3) {
      message("FOUND NO INSIDE POINT")
      message("number of splits")
      message(length(ci.info$centers))
      message("centers")
      message(ci.info$centers)
      message("ses")
      message(ci.info$ses)
      stop("couldn't find an inside point between low and high. The confidence interval doesn't exist!")
    }
  }
  if (verbose) {
    cat("Found an inside point at granularity of", 
        range.length, "\n")
  }
  min(test.range[cover])
}

does.it.cover.gammamin.saturated <- function (beta.j, ci.info) {
  # is beta.j within the interval
  # i.e. is th two-sided pvalue of beta \geq 1 - ci.level
  if (missing(ci.info)) 
    stop("ci.info is missing to the function does.it.cover.gammamin")
  centers <- ci.info$centers
  npv <- length(centers)
  no.inf.ci <- ci.info$no.inf.ci
  ses <- ci.info$ses
  vlo <- ci.info$vlo
  vup <- ci.info$vup
  gamma.min <- ci.info$gamma.min
  multi.corr <- ci.info$multi.corr
  s0 <- ci.info$s0
  alpha <- 1 - ci.info$ci.level
  pvals <- numeric(npv)
  for (i in 1:npv) {
    # calculate the per-split p-values, this can be numerically tricky!
    pvals[i] <- selectiveInference:::tnorm.surv(centers[i], beta.j, ses[i], vlo[i], vup[i])
    if (pvals[i] == 0 || pvals[i] == 1 || is.na(pvals[i])) {
      pvals[i] <- selectiveInference:::tnorm.surv(centers[i], beta.j, ses[i], vlo[i], vup[i], bits = 2)
      if (is.na(pvals[i])) pvals[i] <- selectiveInference:::tnorm.surv(centers[i], beta.j, ses[i], vlo[i], vup[i], bits = 100)
    } 
  }
  if (any(is.na(pvals))) stop("At least one p-value is NA")
  pvals <- pmin(pvals, 1-pvals)
  if (multi.corr) pvals <- pvals * s0
  pval.rank <- rank(pvals, ties.method = "first")
  nsplit <- length(pval.rank) + no.inf.ci
  gamma.b <- pval.rank / nsplit

  alpha.b <- (alpha * gamma.b / (1 - log(gamma.min)))
  if (all(gamma.b < gamma.min)) {
    # we should not get to this point, since the number of infinite CI's is checked
    return(TRUE)
  } else {
    coveredpre <- all(alpha.b[gamma.b >= gamma.min] < 2 * pvals[gamma.b >= gamma.min])
    return(coveredpre)
  }
}

find.bisection.bounds.gammamin.saturated <- function (shouldcover, shouldnotcover, ci.info, verbose, timeout = 10) {
  # find a point that is definitely outside the interval
  reset.shouldnotcover <- FALSE
  if (does.it.cover.gammamin.saturated(beta.j = shouldnotcover, ci.info = ci.info)) {
    reset.shouldnotcover <- TRUE
    if (verbose) 
      cat("finding a new shouldnotcover bound\n")
    start_user_time <- proc.time()[["elapsed"]]
    i <- 0
    while (does.it.cover.gammamin.saturated(beta.j = shouldnotcover, 
                                            ci.info = ci.info)) {
      i <- i + 1
      old <- shouldnotcover
      shouldnotcover <- shouldnotcover + i * (shouldnotcover - shouldcover)
      shouldcover <- old
      if ((proc.time()[["elapsed"]] - start_user_time) > timeout) 
        stop(paste("Searched outside point for more than ", timeout, " seconds, reached ", shouldnotcover))
    }
    if (verbose) {
      cat("new\n")
      cat("shouldnotcover", shouldnotcover, "\n")
      cat("shouldcover", shouldcover, "\n")
    }
  }
  if (!does.it.cover.gammamin.saturated(beta.j = shouldcover, ci.info = ci.info)) {
    if (reset.shouldnotcover) 
      stop("Problem: we first reset shouldnotcover and are now resetting shouldcover, this is not supposed to happen")
    if (verbose) 
      cat("finding a new shouldcover bound\n")
    while (!does.it.cover.gammamin.saturated(beta.j = shouldcover, 
                                   ci.info = ci.info)) {
      old <- shouldcover
      shouldcover <- shouldcover + (shouldcover - shouldnotcover)
      shouldnotcover <- old
    }
    if (verbose) {
      cat("new\n")
      cat("shouldnotcover", shouldnotcover, "\n")
      cat("shouldcover", shouldcover, "\n")
    }
  }
  return(list(shouldcover = shouldcover, shouldnotcover = shouldnotcover))
}

bisection.gammamin.coverage.saturated <- function (outer, inner, ci.info, verbose, eps.bound = 10 ^ (-7)) {
  # find parameter with a p-value of exactly 1-ci.level with a precision of eps.bound (precision of the parameter)
  check.bisection.bounds.gammamin.saturated(shouldcover = inner, shouldnotcover = outer, 
                                  ci.info = ci.info, verbose = verbose)
  eps <- 1
  while (eps > eps.bound) {
    middle <- (outer + inner)/2
    if (does.it.cover.gammamin.saturated(beta.j = middle, ci.info = ci.info)) {
      inner <- middle
    } else {
      outer <- middle
    }
    eps <- abs(inner - outer)
  }
  solution <- (inner + outer)/2
  if (verbose) {
    cat("finished bisection...eps is", eps, "\n")
  }
  return(solution)
}

check.bisection.bounds.gammamin.saturated <- function (shouldcover, shouldnotcover, ci.info, verbose) {
  # this function should not fail
  if (does.it.cover.gammamin.saturated(beta.j = shouldnotcover, ci.info = ci.info)) {
    stop("shouldnotcover bound is covered! we need to decrease it even more! (PLZ implement)")
  } else {
    if (verbose) 
      cat("shouldnotcover bound is not covered, this is good")
  }
  if (does.it.cover.gammamin.saturated(beta.j = shouldcover, ci.info = ci.info)) {
    if (verbose) 
      cat("shouldcover is covered!, It is a good covered bound")
  } else {
    stop("shouldcover is a bad covered bound, it is not covered!")
  }
}


pval.creator <- function(beta, gamma, vlo, vup, centers, ses, s0 = NA, multi.corr = FALSE) {
  B <- length(centers)
  test <- is.finite(vlo)
  npv <- sum(test)
  vlo <- vlo[test]
  vup <- vup[test]
  ses <- ses[test]
  centers <- centers[test]
  no.inf.ci <- B - npv
  gamma.min <- min(gamma)
  pvals <- numeric(npv)
  for (i in 1:npv) {
    pvals[i] <- selectiveInference:::tnorm.surv(centers[i], beta, ses[i], vlo[i], vup[i])
  }
  pvals <- 2 * pmin(pvals, 1 - pvals)
  pvals <- c(pvals, rep(1, no.inf.ci))
  if (multi.corr) {
    if (any(is.na(s0))) 
      stop("need s0 information to be able to create multiple testing corrected pvalues")
    s0 <- s0[test]
    pvals[1:npv] <- pvals[1:npv] * s0
  }
  pval <- pval.aggregator(list(as.matrix(pvals, ncol = 1)), gamma)
  return(pval)
}

print.carve <- function(x){
  if (is.list(x)) {
    if (isTRUE(x$method == "multi.carve" || x$method == "multi.carve.group")) {
      cat("Result from multicarving \n")
    } else if (isTRUE(x$method == "multi.split")){
      cat("Result from multisplitting \n")
    } else if (isTRUE(x$method == "carve100")){
      cat("Result from carve 100 \n")
    } else {
      return (print.default(x))
    }
    if (x$method == "multi.carve.group"){
      cat("alpha = 0.01:")
      cat(" Selected predictor groups:", which(x$pval.corr <= 0.01), 
          "\n")
      cat("alpha = 0.05:")
      cat(" Selected predictor groups:", which(x$pval.corr <= 0.05), 
          "\n")
    } else {
      cat("alpha = 0.01:")
      cat(" Selected predictors:", which(x$pval.corr <= 0.01), 
          "\n")
      cat("alpha = 0.05:")
      cat(" Selected predictors:", which(x$pval.corr <= 0.05), 
          "\n")
    }

    cat("------\n")
    if (isTRUE(x$FWER)){
      cat("Familywise error rate controlled at level alpha.\n")
    } else if (isFALSE(x$FWER)){
      if (x$method == "multi.carve.group"){
        cat("Single group error rate controlled at level alpha.\n")
      } else {
        cat("Single variable error rate controlled at level alpha.\n")
      }
      cat("No multiplicity correction applied.\n")
    }
   
  } else {
    print.default(x)
  }
}