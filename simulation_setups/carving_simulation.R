rm(list = ls(all = TRUE))
save <- TRUE

# create save location, adjust depending on folder structure
if (save) {
  newdir <- format(Sys.time(), "%d-%b-%Y %H.%M")
  dir.create(paste("simulation_setups/multi_carve/", newdir, sep="")) 
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
require(git2r)

commit <- revparse_single(revision = "HEAD")
print(paste("Run on commit", commit$sha, 'i.e.:', commit$summary))


# adjust depending on folder structure
source("inference/hdi_adjustments.R")
source("inference/carving.R")
source("inference/sample_from_truncated.R")
source("inference/tryCatch-W-E.R")



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
print (x[1,1])
# should create the right x on D-MATH server, x[1 ,1] = 0.958 for Toeplitz 0.6
y.true <- x %*% beta
SNR <- 1.713766 # value created for Toeplitz 0.6
sigma <- sqrt(drop(var(y.true)) / SNR)
if (rho == 0.6) sigma <- 2
report.sigma <- FALSE

# # riboflavin
# # adjust depending on folder structure
# riboflavin <- read.csv("riboflavin.csv",
#                          stringsAsFactors = FALSE)
# riboflavin.tmp <- t(riboflavin[, -1])
# colnames(riboflavin.tmp) <- riboflavin[, 1]
# x <- riboflavin.tmp[, -1]
# rm(riboflavin)
# rm(riboflavin.tmp)
# n <- dim(x)[1]
# p <- dim(x)[2]
# report.sigma <- FALSE
# SNR <- 16
# sparsity <- 2 # 4 in other set-up

B.vec <- c(1, 50) # c(1, (1:5) * 10) # number of splits
frac.vec <- c(0.5, 0.75, 0.9, 0.95, 0.99) # selection fraction
nsim <- 200
ntasks <- nsim
progress <- function(n, tag) {
  mod <- 16
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
seed.vec <- sample(1:10000, length(frac.vec))
print(seed.vec) # 3588 3052 2252 5257 8307
seed.n <- 0
B <- max(B.vec)
for (frac in frac.vec) {
  seed.n <- seed.n + 1
  set.seed(seed.vec[seed.n])
  # check set-up
  print(frac)
  print(B)
  print(report.sigma)
  print(sigma)
  opts <- list(progress = progress)
  
  # parallelization
  # choose different number of cores if wished
  cl<-makeSOCKcluster(16) 
  
  rseed <- seed.vec[seed.n]
  clusterSetRNGStream(cl, iseed = rseed) #make things reproducible
  registerDoSNOW(cl)
  tic()
  res<-foreach(gu = 1:nsim, .combine = rbind,
               .packages = c("MASS", "selectiveInference", "glmnet", "Matrix",
                             "hdi", "tmg", "truncnorm", "tictoc"), .options.snow = opts) %dorng%{
  # alternative if sequential computation is preferred
  # res <- foreach(gu = 1:nsim, .combine = rbind) %do%{

    # # Riboflavin
    # beta <- rep(0, p)
    # ind <- sample(1:p, sparsity)
    # beta[ind] <- 1
    # y.true <- x %*% beta
    # sigma <- sqrt(drop(var(y.true)) / SNR)
    
    y <- y.true + sigma * rnorm(n)

    reported.sigma <- NA
    if (report.sigma) {
      reported.sigma <- sigma
    } else {
      # not acutally necessary if sigma is estimated within the routines
      estSigma <- estimateSigma.flex(scale(x, T, F), scale(y, T, F),
                                    intercept = FALSE, standardize = FALSE)
      reported.sigma <- estSigma$sigmahat
    }

    mcrtry <- tryCatch_W_E(multi.carve(x, y, B = B, fraction = frac, model.selector = lasso.cvcoef,
                                       classical.fit = lm.pval.flex, parallel = FALSE,
                                       ncores = getOption("mc.cores", 2L),
                                       args.model.selector = list(standardize = FALSE, intercept = TRUE, tol.beta = 0, use.lambda.min = FALSE),
                                       args.classical.fit = list(Sigma = reported.sigma, t.test = FALSE), verbose = FALSE,
                                       FWER = FALSE, split.pval = TRUE, return.selmodels = TRUE, return.nonaggr = TRUE,
                                       args.lasso.inference = list(sigma = reported.sigma,
                                                                                             verbose = TRUE, selected = TRUE)), 0)
    c100try <- tryCatch_W_E(carve100(x, y, model.selector = lasso.cvcoef,
                                     args.model.selector = list(standardize = FALSE, intercept = TRUE, tol.beta = 1e-5, use.lambda.min = FALSE),
                                     verbose = FALSE, FWER = FALSE, return.selmodels = TRUE,
                                     estimate.sigma = FALSE, args.lasso.inference = list(sigma = reported.sigma)), 0)
    
    out.list <- list()
    out.list$y <- y
    if (!is.null(mcrtry$error) || !is.null(c100try$error)) {
      # error handling
      err <- paste("mcr:", mcrtry$error, "carve100:", c100try$error)
      war <- if (is.null(mcrtry$warning) || is.null(c100try$warning)) NA
      else c(mcrtry$warning, c100try$warning)
      for (b in B.vec) {
        if (b > 1) {
          out.list[[as.character(b)]] <- rep(NA, (length(ind) + 1) * 16)
        } else {
          out.list[[as.character(b)]] <- rep(NA, (length(ind) + 1) * 6 + 9)
        }
      }
      out.list$exception <- list(err, paste(1:length(war), ":", war, collapse = ", "))
      out.list
    } else {
      mcr <- mcrtry$value
      pcarve.nofwer <- mcr[[1]]$pvals.nonaggr
      psplit.nofwer <- mcr[[2]]$pvals.nonaggr
      model.size <- apply(mcr[[1]]$sel.models, 1, sum)
      model.size[model.size == 0]  <- 1 # if no variable is selected, p-values are 1
      # ommit the clipping to calculate adjusted power
      # pcarve.fwer <- pmin(pcarve.nofwer * model.size, 1)
      # psplit.fwer <- pmin(psplit.nofwer * model.size, 1)
      pcarve.fwer <- pcarve.nofwer * model.size
      psplit.fwer <- psplit.nofwer * model.size
      c100 <- c100try$value
      pc100.nofwer <- c100$pval.corr
      model.size100 <- sum(c100$sel.models)
      model.size100[model.size100 == 0]  <- 1
      # ommit the clipping to calculate adjusted power
      # pc100.fwer <- pmin(pc100.nofwer * model.size100, 1)
      pc100.fwer <- pc100.nofwer * model.size100
      
      for (B in B.vec) {
        if (B > 1) {
          use <- 1:B
          pvals.aggregated <- pval.aggregator(list(pcarve.nofwer[use, ], pcarve.fwer[use, ], psplit.nofwer[use, ], psplit.fwer[use, ]),
                                              round(seq(ceiling(0.05 * B)/B, 1, by = 1/B), 2), cutoff = FALSE)
          pvals.aggregated2 <- pval.aggregator(list(pcarve.nofwer[use, ], pcarve.fwer[use, ], psplit.nofwer[use, ], psplit.fwer[use, ]),
                                               round(seq(ceiling(0.3 * B)/B, 1, by = 1/B), 2), cutoff = FALSE)
          pvals.aggregated3 <- pval.aggregator(list(pcarve.nofwer[use, ], pcarve.fwer[use, ], psplit.nofwer[use, ], psplit.fwer[use, ]),
                                               round(ceiling(0.05 * B)/B, 2), cutoff = FALSE)
          pvals.aggregated4 <- pval.aggregator(list(pcarve.nofwer[use, ], pcarve.fwer[use, ], psplit.nofwer[use, ], psplit.fwer[use, ]),
                                               round(ceiling(0.3 * B)/B, 2), cutoff = FALSE)
        } else {
          pvals.aggregated <- list(pcarve.nofwer[1, ], pcarve.fwer[1, ], psplit.nofwer[1, ], psplit.fwer[1, ])
        }
        
        
        run.res <- vector(length = 0) # store important quantities
        np <- length(pvals.aggregated)
        for (i in 1:np) {
          true.pv <- pvals.aggregated[[i]][ind] # p-values of active variables
          bad.pv <- min(pvals.aggregated[[i]][-ind]) # lowest p-value of inactive variables to check FWER
          run.res <- c(run.res, true.pv, bad.pv)
        }
        if (B > 1) {
          # for multicarving, test different aggregation methods
          for (i in 1:np) {
            true.pv <- pvals.aggregated2[[i]][ind]
            bad.pv <- min(pvals.aggregated2[[i]][-ind])
            run.res <- c(run.res, true.pv, bad.pv)
          }
          for (i in 1:np) {
            true.pv <- pvals.aggregated3[[i]][ind]
            bad.pv <- min(pvals.aggregated3[[i]][-ind])
            run.res <- c(run.res, true.pv, bad.pv)
          }
          for (i in 1:np) {
            true.pv <- pvals.aggregated4[[i]][ind]
            bad.pv <- min(pvals.aggregated4[[i]][-ind])
            run.res <- c(run.res, true.pv, bad.pv)
          }
        }
        if (B == 1) {
          # analyse first split specially for B = 1 and analyse carve100
          R <- length(which(mcr[[1]]$sel.models[1, ])) # number of variables selected in first split
          TS <- sum(ind %in% which(mcr[[1]]$sel.models[1, ])) # number of active variables selected
          V <- R - TS # number of inactive variables selected
          carve.err <- sum(pvals.aggregated[[1]][-ind] < 0.05) # number of false rejection from single-carving
          split.err <- sum(pvals.aggregated[[3]][-ind] < 0.05) # number of false rejection from single-splitting
          carve100.err <- sum(pc100.nofwer[-ind] < 0.05) # number of false rejection from carve100
          run.res <- c(run.res, R, V, TS) 
          true.pv <- pc100.nofwer[ind] # p-values of active variables
          bad.pv <- min(pc100.nofwer[-ind]) # lowest p-value of inactive variables to check FWER
          run.res <- c(run.res, true.pv, bad.pv)
          true.pv <- pc100.fwer[ind]
          bad.pv <- min(pc100.fwer[-ind])
          R100 <- length(which(c100$sel.models)) # number of variables selected using all data
          TS100 <- sum(ind %in% which(c100$sel.models)) # number of active variables selected
          V100 <- R100 - TS100 # number of inactive variables selected
          run.res <- c(run.res, true.pv, bad.pv, R100, V100, TS100,
                       carve.err, split.err, carve100.err)
        }
        out.list[[as.character(B)]] <- run.res
      }
      err <- if (is.null(mcrtry$error) && is.null(c100try$error)) NA
      else c(mcrtry$error, c100try$error) # should not happen due to earlier check
      war <- if (is.null(mcrtry$warning) && is.null(c100try$warning)) NA
      else c(mcrtry$warning, c100try$warning)
      out.list$exception <- list(err, paste(1:length(war), ":", war, collapse = ", "))
      out.list
      # end of simulation run
    }
    # end of simulation for given fraction
  }
  toc()
  stopCluster(cl)
  
  # analyse results for given fraction
  # get matrix of errors and warnings
  expmatr <- matrix(unlist(res[, "exception"]), nrow = dim(res)[1],
                    ncol = 2, byrow = TRUE)
  print(sum(is.na(expmatr[, 1])))
  succ = which(is.na(expmatr[, 1]))
  print("succesful runs")
  
  all.y <- matrix(unlist(res[,"y"]), nrow = dim(res), byrow = TRUE)
  sd <- attr(res, "rng")
  
  for (B in B.vec) {
    if (B == 1) {
      names1 <- c("carve","carvefw", "split", "splitfw")
      names2 <- c("carve100","carve100fw")
      subres <- matrix(unlist(res[,as.character(B)]), nrow = dim(res)[1],
                       ncol = 6 * (sparsity + 1) + 9, byrow = TRUE)
      if (any(!is.na(subres[-succ, ]))) print("not as it should be") # sanity check
      subres <- subres[succ,]
      colnames(subres) <- c(rep(names1, each = (sparsity + 1)), "R", "V", "R-V",
                            rep(names2, each = (sparsity + 1)), "R100", "V100",
                            "R-V100", "carve.err", "split.err", "carve100.err")
      names <- c(names1, names2)
    } else {
      names <- c("carve5", "carvefw5", "split5", "splitfw5", "carve30",
               "carvefw30", "split30", "splitfw30", "carvefix5",
               "carvefwfix5", "splitfix5", "splitfwfix5", "carvefix30",
               "carvefwfix30", "splitfix30", "splitfwfix30")
      subres <- matrix(unlist(res[,as.character(B)]), nrow = dim(res)[1],
                       ncol = 16 * (sparsity + 1), byrow = TRUE)
      if (any(!is.na(subres[-succ, ]))) print("not as it should be")
      subres <- subres[succ,]
      colnames(subres) <- c(rep(names, each = (sparsity + 1)))
    }
    subres <- as.data.frame(subres)

    simulation <- list("results" = subres, "exceptions" = expmatr, "y" = all.y, "B" = B, "split" = frac,
                       "nsim" = nsim, "seed" = rseed, "All used B" = B.vec, "sd" = sd, "commit" = commit)
    print(paste("results using fraction ", frac, " and B=", B, sep = ""))
    if (B == 1) {
      print(mean(subres$`R-V` == sparsity)) # probability of screening
      good <- which(subres$`R-V` == sparsity)
      print(apply(subres[, c("R", "V", "R-V")], 2, mean)) # totally, active, inactive selected
      # probability of screening using all data for selection
      print(mean(subres$`R-V100` == sparsity)) 
      good100 <- which(subres$`R-V100` == sparsity)
      # totally, active and inacitve selected using all data for selection
      print(apply(subres[, c("R100", "V100", "R-V100")], 2, mean))
      print(c(c(sum(subres$carve.err), sum(subres$split.err)) / sum(subres$V),
              sum(subres$carve100.err) / sum(subres$V100))) # Rejection amongst falsely selected
      # Rejection amongst falsely selected conditioned on screening, should be below 0.05
      print(c(c(sum(subres$carve.err[good]), sum(subres$split.err[good])) / sum(subres$V[good]),
              sum(subres$carve100.err[good100]) / sum(subres$V100[good100]))) 
    } 
    allrej <- matrix(NA, nrow = 1,ncol = length(names))
    colnames(allrej) <- names
    for (name in names) {
      nameind <- which(colnames(subres) == name)
      mat <- subres[, nameind]
      rej <- quantile(mat[, sparsity + 1], 0.05, na.rm = TRUE)
      rejmat <- mat[, 1:sparsity] < rej
      allrej[, as.character(name)] <- mean(rejmat, na.rm = TRUE)
      
    }
    print("Adjusted")
    print(allrej) # adjusted power
    fwer <- numeric(length(names))
    names(fwer) <- names
    for (name in names) {
      nameind <- which(colnames(subres) == name)
      mat <- subres[, nameind]
      fwer[as.character(name)] <- mean(mat[, sparsity + 1] < 0.05, na.rm = TRUE)
      rejmat <- mat[, 1:sparsity] < 0.05
      allrej[, as.character(name)] <- mean(rejmat, na.rm = TRUE)
    }
    print("Unadjusted")
    print(fwer) # FWER
    print(allrej) # power
    resname <- paste0("results ", format(Sys.time(), "%d-%b-%Y %H.%M"),
                      " split=", frac, " B=", B, " seed=", rseed)
    # adjust depending on folder structure
    if (save) save(simulation, file = paste("simulation_setups/multi_carve/", newdir, "/", resname, ".RData", sep = ""))
    # end of analysis for given B
  }
  # end of analysis for given fraction
}

print("Finale")



