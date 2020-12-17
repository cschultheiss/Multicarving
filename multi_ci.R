rm(list = ls(all = TRUE))
local <- FALSE # on local machine or on D-MATH server
save <- FALSE
if (!save && !local) {
  save <- TRUE
  warning("Do not use cluster without saving, save was set to TRUE")
}
# create save location, adjust depending on folder structure
if (save) {
  newdir <- format(Sys.time(), "%d-%b-%Y %H.%M")
  if (local) {
    dir.create(paste("C:/Users/Christoph/Documents/ETH/MA/R_trials/multi_carve/", newdir, sep = ""))
  } else {
    dir.create(paste("multi_carve/",newdir,sep="")) 
  }
}

require(MASS)
require(glmnet)
require(Matrix)
require(tictoc)
require(hdi)
require(selectiveInference)
require(doSNOW)
require(parallel)
require(doRNG)
require(tmg)
require(truncnorm)


if (local) {
  # adjust depending on folder structure
  source("C:/Users/Christoph/Documents/ETH/MA/R_trials/hdi_adjustments.R")
  source("C:/Users/Christoph/Documents/ETH/MA/R_trials/optimal_inference.R")
  source("C:/Users/Christoph/Documents/ETH/MA/R_trials/tryCatch-W-E.R")
} else {
  source("hdi_adjustments.R")
  source("optimal_inference.R")
  source("tryCatch-W-E.R")
}

if (local) {
  library(tcltk)
}


# toeplitz
n <- 100
p <- 200
rho <- 0.6
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
sel.index <- c(1, 5, 10, 15, 20)
ind <- sel.index
beta <- rep(0, p)
beta[sel.index] <- 1
sparsity <- 5
set.seed(42) # to make different methods comparable, fix the x-matrix
x <- mvrnorm(n, rep(0, p), Cov) 
# should create the right x on D-MATH server, x[1 ,1] = 0.958

# if x-data is presaved
# if (local) {
# load("~/ETH/MA/R_trials/x.RData")
# } else {
#   load("x.RData")
# }
y_true <- x %*% beta
sigma <- 2
report_sigma <- FALSE

frac_vec <- c(0.5, 0.75, 0.9, 0.95, 0.99) # selection fraction
nsim <- 200
ntasks <- nsim
progress <- function(n, tag) {
  mod <- if (local) 4
  else
    16
  if (local) {
    setTkProgressBar(pb, n)
  }
  if (n %% mod == 0 ) {
    cat(sprintf('tasks completed: %d; tag: %d\n', n, tag))
  }
  if (n %% mod == 0 ) {
    toc()
    tic()
  }
}
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
seed_v <- sample(1:10000,length(frac_vec))
print(seed_v) # 3588 3052 2252 5257 8307
seed_n <- 0
usei <- TRUE # should sigma_hat be calculated modelwise
B <- 50
for (frac in frac_vec) {
  seed_n <- seed_n+1
  set.seed(seed_v[seed_n])
  # check set-up
  print(usei)
  print(frac)
  print(B)
  print(report_sigma)
  if (local) {
    pb <- tkProgressBar(max = ntasks, title = paste("B=", B, " split=", frac))
  }
  opts <- list(progress = progress)
  
  # parallelization
  if (local) {
    cl<-makeSOCKcluster(4)
  } else {
    cl<-makeSOCKcluster(16) 
  }
  rseed <- seed_v[seed_n]
  clusterSetRNGStream(cl, iseed = rseed) #make things reproducible
  registerDoSNOW(cl)
  tic()
  res<-foreach(gu=1:nsim, .combine = rbind,
               .packages = c("MASS","selectiveInference","glmnet","Matrix",
                             "hdi","tmg","truncnorm"),.options.snow=opts) %dorng%{
   # alternative if sequential computation is preferred
   # res<-foreach(gu=1:nsim,.combine = rbind) %do%{
    y_true <- x %*% beta
    y <- y_true + sigma * rnorm(n)
    reported_sigma <- NA
    gammavec <- round(seq(ceiling(0.05 * B) / B, 1 - 1 / B, by = 1 / B), 2)
    gammavec[1] <- 0.05999 # due to inconsistency for multisplitting CI
    mcrtry <- tryCatch_W_E(multi.carve.ci.saturated(x, y, B = B, fraction = frac, ci.level = 0.95, 
                                     model.selector = lasso.cvcoef, classical.fit = lm.pval.flex,
                                     parallel = FALSE, ncores = getOption("mc.cores", 2L), gamma = gammavec,
                                     args.model.selector = list(standardize = FALSE, intercept = TRUE,
                                                                tol.beta = 0, use_lambda.min = FALSE),
                                     args.classical.fit = list(Sigma = reported_sigma),
                                     verbose = FALSE, ci.timeout = 10, FWER = FALSE, split_pval = TRUE,
                                     return.selmodels = TRUE, return.nonaggr = TRUE,
                                     args.lasso.inference = list(sigma = reported_sigma, verbose = FALSE)),0)
    
    if(!is.null(mcrtry$error)) {
      err <- paste("mcr:", mcrtry$error)
      war <- if (is.null(mcrtry$warning)) NA
      else c(mcrtry$warning)
      out_list <- list()
      out_list$splitlci <- rep(NA,p)
      out_list$splituci <- rep(NA,p)
      out_list$carvelci <- rep(NA,p)
      out_list$carveuci <- rep(NA,p)
      out_list$betamin <- rep(NA,p)
      out_list$betamax <- rep(NA,p)
      out_list$exception <- list(err, paste(1:length(war), ":", war, collapse = ", "))
      out_list
    } else {
      mcr <- mcrtry$value
      
      betamin <- rep(-Inf,p)
      betamax <- rep(Inf,p)
      for (j in 1:p) {
        betas <- rep(NA,50)
        for (i in 1:B) {
          if (mcr[[1]]$sel.models[i, j]) {
            mod <- which(mcr[[1]]$sel.models[i, ])
            coefn <- which(mod == j)
            xsub <- cbind(rep(1, n), x[, mod])
            betas[i] <- (ginv(xsub) %*% y_true)[coefn + 1]
          }
        }
        if (any(!is.na(betas))) {
          betamin[j] <- min(betas, na.rm = TRUE)
          betamax[j] <- max(betas, na.rm = TRUE)
        }
      }
      slci <- mcr[[2]]$lci
      suci <- mcr[[2]]$uci
      clci <- mcr[[1]]$lci
      cuci <- mcr[[1]]$uci

      out_list <- list()
      out_list$splitlci <- slci
      out_list$splituci <- suci
      out_list$carvelci <- clci
      out_list$carveuci <- cuci
      out_list$betamin <- betamin
      out_list$betamax <- betamax
      err <- if (is.null(mcrtry$error)) NA
      else c(mcrtry$error) # should not happen due to earlier check
      war <- if (is.null(mcrtry$warning)) NA
      else c(mcrtry$warning)
      out_list$exception <- list(err, paste(1:length(war), ":", war, collapse = ", "))
      out_list
      # end of simulation run
    }
    # end of simulation for given fraction
  }  
  toc()
  stopCluster(cl)
  if (local) {
    close(pb)
  }

  expmatr <- matrix(unlist(res[, "exception"]),nrow = dim(res)[1],
                    ncol = 2,byrow = TRUE)
  print(sum(is.na(expmatr[, 1])))
  succ <- which(is.na(expmatr[, 1]))
  print("succesful runs")
  
  splitlci <- matrix(unlist(res[, "splitlci"]), nrow = dim(res)[1],
                     ncol = p, byrow = TRUE)
  splituci <- matrix(unlist(res[, "splituci"]), nrow = dim(res)[1],
                     ncol = p, byrow = TRUE)
  carvelci <- matrix(unlist(res[, "carvelci"]), nrow = dim(res)[1],
                     ncol = p, byrow = TRUE)
  carveuci <- matrix(unlist(res[, "carveuci"]), nrow = dim(res)[1],
                     ncol = p, byrow = TRUE)
  betamin <- matrix(unlist(res[, "betamin"]), nrow = dim(res)[1],
                    ncol = p, byrow = TRUE)
  betamax <- matrix(unlist(res[, "betamax"]), nrow = dim(res)[1],
                    ncol = p, byrow = TRUE)
  
  colnames(splitlci) <- paste(1:p)
  colnames(splituci) <- paste(1:p)
  colnames(carvelci) <- paste(1:p)
  colnames(carveuci) <- paste(1:p)
  colnames(betamin) <- paste(1:p)
  colnames(betamax) <- paste(1:p)
  spliterr <- (splitlci > betamin | splituci < betamax)
  print(mean(apply(spliterr, 1, sum)))
  carveerr <- (carvelci > betamin | carveuci < betamax)
  print(mean(apply(carveerr, 1, sum)))
  splitdis <- splituci - splitlci
  splitdis[spliterr] <- NA
  carvedis <- carveuci - carvelci
  carvedis[carveerr] <- NA
  print(quantile(splitdis, seq(0.01, 0.1, 0.01), na.rm = TRUE))
  print(quantile(carvedis, seq(0.01, 0.1, 0.01), na.rm = TRUE))
  results <- list("splitlci" = splitlci, "splituci" = splituci, "carvelci" = carvelci,
                  "carveuci" = carveuci, "betamin" = betamin, "betamax" = betamax)
  
  simulation <- list("results" = results, "exceptions" = expmatr, "B" = B,
                     "split" = frac, "nsim" = nsim, "seed" = rseed)
  resname <- paste0("results ", format(Sys.time(), "%d-%b-%Y %H.%M"), " split=", frac, " B=", B, " seed=", rseed)
  # adjust depending on folder structure
  if (save) {
    if (local) {
      save(simulation, file = paste("C:/Users/Christoph/Documents/ETH/MA/R_trials/multi_carve/",
                                    newdir, "/", resname, ".RData", sep = ""))
    } else {
      save(simulation,file = paste("multi_carve/", newdir, "/", resname, ".RData", sep = ""))
    }
  }
}

      
      