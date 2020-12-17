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
require(git2r)

commit <- revparse_single(revision = "HEAD")
print(paste("Run on commit", commit$sha, 'i.e.:', commit$summary))


if (local) {
  # adjust depending on folder structure
  source("C:/Users/Christoph/Documents/ETH/MA/R_trials/hdi_adjustments.R")
  source("C:/Users/Christoph/Documents/ETH/MA/R_trials/optimal_inference.R")
  source("C:/Users/Christoph/Documents/ETH/MA/R_trials/tryCatch-W-E.R")
} else {
  source("hdi_adjustments.R")
  source("optimal_inference.R")
  source("tryCatch-W-E.R")
  source("mutli_restructured.R")
  source("sample_from_truncated.R")
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
print (x[1,1])
# should create the right x on D-MATH server, x[1 ,1] = 0.958 for toeplitz 0.6
# # if x-data is presaved
# if (local) {
# load("~/ETH/MA/R_trials/x.RData")
# } else {
#   load("x.RData")
# }
y_true <- x %*% beta
SNR <- 1.713766 # value created for Toeplitz 0.6
sigma <- sqrt(drop(var(y_true)) / SNR)
if (rho==0.6) sigma <- 2
report_sigma <- FALSE

# # riboflavin
# # adjust depending on folder structure
# if (local) {
#   riboflavin <- read.csv("C:/Users/Christoph/Documents/ETH/MA/data/riboflavin.csv",
#                          stringsAsFactors = FALSE)
# } else {
#   riboflavin <- read.csv("riboflavin.csv",
#                          stringsAsFactors = FALSE)
# }
# riboflavin.tmp <- t(riboflavin[, -1])
# colnames(riboflavin.tmp) <- riboflavin[, 1]
# x <- riboflavin.tmp[, -1]
# rm(riboflavin)
# rm(riboflavin.tmp)
# n <- dim(x)[1]
# p <- dim(x)[2]
# report_sigma <- FALSE
# SNR <- 16
# sparsity <- 2 # 4 in other set-up

# # geno data, not in Paper
# # adjust depending on folder structure
# if (local) {
#   ind1 <- readRDS("C:/Users/Christoph/Documents/ETH/MA/Data/ind1.rds")
#   geno <- readRDS("C:/Users/Christoph/Documents/ETH/MA/Data/geno-n600-p1000.rds")
# } else {
#   ind1 <- readRDS("ind1.rds")
#   geno <- readRDS("geno-n600-p1000.rds")
# }
# x <- geno[ind1,]
# rm(geno)
# rm(ind1)
# n <- dim(x)[1]
# p <- dim(x)[2]
# report_sigma <- FALSE
# SNR <-2
# sparsity <- 10

B_vec <- c(1, 50) # c(1, (1:5) * 10) # number of splits
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
usei <- FALSE # should sigma_hat be calculated modelwise, change for Riboflavin
B <- max(B_vec)
for (frac in frac_vec) {
  seed_n <- seed_n+1
  set.seed(seed_v[seed_n])
  # check set-up
  print(usei)
  print(frac)
  print(B)
  print(report_sigma)
  print(sigma)
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

    # # Riboflavin or Geno
    # beta <- rep(0, p)
    # ind <- sample(1:p, sparsity)
    # beta[ind] <- 1
    # y_true <- x %*% beta
    # sigma <- sqrt(drop(var(y_true)) / SNR)
    
    y <- y_true + sigma * rnorm(n)

    reported_sigma <- NA
    if (report_sigma) {
      reported_sigma <- sigma
    } else {
      # not acutally necessary if usei = TRUE
      estSigma <- estimateSigma.flex(scale(x, T, F), scale(y, T, F),
                                    intercept = FALSE, standardize = FALSE)
      reported_sigma <- estSigma$sigmahat
    }

    mcrtry <- tryCatch_W_E(multi.carve(x, y, B = B, fraction = frac, model.selector = lasso.cvcoef,
                                       classical.fit = lm.pval.flex, parallel = FALSE,
                                       ncores = getOption("mc.cores", 2L), # aggregate outside to test different methods
                                       args.model.selector = list(standardize = FALSE, intercept = TRUE, tol.beta = 0, use_lambda.min = FALSE),
                                       args.classical.fit = list(Sigma = reported_sigma, ttest = FALSE), verbose = FALSE,
                                       FWER = FALSE, split_pval = TRUE, return.selmodels = TRUE, return.nonaggr = TRUE,
                                       args.lasso.inference = list(sigma = reported_sigma,
                                                                                             verbose = TRUE, selected = TRUE)), 0)
    c100try <- tryCatch_W_E(carve100(x, y, model.selector = lasso.cvcoef,
                                     args.model.selector = list(standardize = FALSE, intercept = TRUE, tol.beta = 1e-5, use_lambda.min = FALSE),
                                     verbose = FALSE, FWER = FALSE, return.selmodels = TRUE,
                                     estimate.sigma = FALSE, args.lasso.inference = list(sigma = reported_sigma)), 0)
    
    out_list <- list()
    out_list$y <- y
    if (!is.null(mcrtry$error) || !is.null(c100try$error)) {
      # error handling
      err <- paste("mcr:", mcrtry$error, "carve100:", c100try$error)
      war <- if (is.null(mcrtry$warning) || is.null(c100try$warning)) NA
      else c(mcrtry$warning, c100try$warning)
      for (b in B_vec) {
        if (b > 1) {
          out_list[[as.character(b)]] <- rep(NA, (length(ind) + 1) * 16)
        } else {
          out_list[[as.character(b)]] <- rep(NA, (length(ind) + 1) * 6 + 9)
        }
      }
      out_list$exception <- list(err, paste(1:length(war), ":", war, collapse = ", "))
      out_list
    } else {
      mcr <- mcrtry$value
      pcarve_nofwer <- mcr[[1]]$pvals.nonaggr
      psplit_nofwer <- mcr[[2]]$pvals.nonaggr
      model_size <- apply(mcr[[1]]$sel.models, 1, sum)
      model_size[model_size == 0]  <- 1 # if no variable is selected, p-values are 1
      # ommit the clipping to calculate adjusted power
      # pcarve_fwer <- pmin(pcarve_nofwer * model_size, 1)
      # psplit_fwer <- pmin(psplit_nofwer * model_size, 1)
      pcarve_fwer <- pcarve_nofwer*model_size
      psplit_fwer <- psplit_nofwer*model_size
      c100 <- c100try$value
      pc100_nofwer <- c100$pval.corr
      model_size100 <- sum(c100$sel.models)
      model_size100[model_size100 == 0]  <- 1
      # pc100_fwer <- pmin(pc100_nofwer * model_size100, 1)
      pc100_fwer <- pc100_nofwer * model_size100
      
      for (B in B_vec) {
        if (B > 1) {
          use <- 1:B
          pvals.aggregated <- pval.aggregator(list(pcarve_nofwer[use, ], pcarve_fwer[use, ], psplit_nofwer[use, ], psplit_fwer[use, ]),
                                              round(seq(ceiling(0.05 * B)/B, 1, by = 1/B), 2), cutoff = FALSE)
          pvals.aggregated2 <- pval.aggregator(list(pcarve_nofwer[use, ], pcarve_fwer[use, ], psplit_nofwer[use, ], psplit_fwer[use, ]),
                                               round(seq(ceiling(0.3 * B)/B, 1, by = 1/B), 2), cutoff = FALSE)
          pvals.aggregated3 <- pval.aggregator(list(pcarve_nofwer[use, ], pcarve_fwer[use, ], psplit_nofwer[use, ], psplit_fwer[use, ]),
                                               round(ceiling(0.05 * B)/B, 2), cutoff = FALSE)
          pvals.aggregated4 <- pval.aggregator(list(pcarve_nofwer[use, ], pcarve_fwer[use, ], psplit_nofwer[use, ], psplit_fwer[use, ]),
                                               round(ceiling(0.3 * B)/B, 2), cutoff = FALSE)
        } else {
          pvals.aggregated <- list(pcarve_nofwer[1, ], pcarve_fwer[1, ], psplit_nofwer[1, ], psplit_fwer[1, ])
        }
        
        
        run_res <- vector(length = 0) # store important quantities
        np <- length(pvals.aggregated)
        for (i in 1:np) {
          true_pv <- pvals.aggregated[[i]][ind] # p-values of active variables
          bad_pv <- min(pvals.aggregated[[i]][-ind]) # lowest p-value of inactive variables to check FWER
          run_res <- c(run_res, true_pv, bad_pv)
        }
        if (B > 1) {
          for (i in 1:np) {
            true_pv <- pvals.aggregated2[[i]][ind]
            bad_pv <- min(pvals.aggregated2[[i]][-ind])
            run_res <- c(run_res, true_pv, bad_pv)
          }
          for (i in 1:np) {
            true_pv <- pvals.aggregated3[[i]][ind]
            bad_pv <- min(pvals.aggregated3[[i]][-ind])
            run_res <- c(run_res, true_pv, bad_pv)
          }
          for (i in 1:np) {
            true_pv <- pvals.aggregated4[[i]][ind]
            bad_pv <- min(pvals.aggregated4[[i]][-ind])
            run_res <- c(run_res, true_pv, bad_pv)
          }
        }
        if (B == 1) {
          # analyse first split specially for B = 1 and analyse carve100
          R <- length(which(mcr[[1]]$sel.models[1, ])) # number of variables selected in first split
          TS <- sum(ind %in% which(mcr[[1]]$sel.models[1, ])) # number of active variables selected
          V <- R - TS # number of inactive variables selected
          carve_err <- sum(pvals.aggregated[[1]][-ind] < 0.05)
          split_err <- sum(pvals.aggregated[[3]][-ind] < 0.05)
          carve100_err <- sum(pc100_nofwer[-ind] < 0.05)
          run_res <- c(run_res, R, V, TS) 
          true_pv <- pc100_nofwer[ind]
          bad_pv <- min(pc100_nofwer[-ind])
          R100 <- length(which(c100$sel.models))
          TS100 <- sum(ind %in% which(c100$sel.models))
          V100 <- R100 - TS100
          run_res <- c(run_res, true_pv, bad_pv)
          true_pv <- pc100_fwer[ind]
          bad_pv <- min(pc100_fwer[-ind])
          run_res <- c(run_res, true_pv, bad_pv, R100, V100, TS100,
                       carve_err, split_err, carve100_err)
        }
        out_list[[as.character(B)]] <- run_res
      }
      err <- if (is.null(mcrtry$error) && is.null(c100try$error)) NA
      else c(mcrtry$error, c100try$error) # should not happen due to earlier check
      war <- if (is.null(mcrtry$warning) && is.null(c100try$warning)) NA
      else c(mcrtry$warning, c100try$warning)
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
  # analyse results for given fraction
  expmatr <- matrix(unlist(res[, "exception"]), nrow = dim(res)[1],
                    ncol = 2, byrow = TRUE)
  print(sum(is.na(expmatr[, 1])))
  succ = which(is.na(expmatr[, 1]))
  print("succesful runs")
  all_y <- matrix(unlist(res[,"y"]), nrow = dim(res), byrow = TRUE)
  for (B in B_vec) {
    if (B == 1) {
      names1 <- c("carve","carvefw", "split", "splitfw")
      names2 <- c("carve100","carve100fw")
      subres <- matrix(unlist(res[,as.character(B)]), nrow = dim(res)[1],
                       ncol = 6 * (sparsity + 1) + 9, byrow = TRUE)
      if (any(!is.na(subres[-succ, ]))) print("not as it should be") # sanity check
      subres <- subres[succ,]
      colnames(subres) <- c(rep(names1, each = (sparsity + 1)), "R", "V", "R-V",
                            rep(names2, each = (sparsity + 1)), "R100", "V100",
                            "R-V100", "carve_err", "split_err", "carve100_err")
      names <- c(names1, names2)
    } else {
      names<-c("carve5", "carvefw5", "split5", "splitfw5", "carve30",
               "carvefw30", "split30", "splitfw30", "carvefix5",
               "carvefwfix5", "splitfix5", "splitfwfix5", "carvefix30",
               "carvefwfix30", "splitfix30", "splitfwfix30")
      subres <- matrix(unlist(res[,as.character(B)]), nrow = dim(res)[1],
                       ncol = 16 * (sparsity+1), byrow = TRUE)
      if (any(!is.na(subres[-succ, ]))) print("not as it should be")
      subres <- subres[succ,]
      colnames(subres) <- c(rep(names, each = (sparsity + 1)))
    }
    subres <- as.data.frame(subres)

    simulation <- list("results" = subres, "exceptions" = expmatr, "y" = all_y, "B" = B, "split" = frac,
                       "nsim" = nsim, "seed" = rseed, "All used B" = B_vec, "estimated sigma" = usei, "commit" = commit)
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
      print(c(c(sum(subres$carve_err), sum(subres$split_err)) / sum(subres$V),
              sum(subres$carve100_err) / sum(subres$V100))) # Rejection amongst falsely selected
      # Rejection amongst falsely selected conditioned on screening, should be below 0.05
      print(c(c(sum(subres$carve_err[good]), sum(subres$split_err[good])) / sum(subres$V[good]),
              sum(subres$carve100_err[good100]) / sum(subres$V100[good100]))) 
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
    print("Unadjusted")
    fwer <- numeric(length(names))
    names(fwer) <- names
    for (name in names) {
      nameind <- which(colnames(subres) == name)
      mat <- subres[, nameind]
      fwer[as.character(name)] <- mean(mat[,sparsity + 1] < 0.05, na.rm = TRUE)
      rejmat <- mat[, 1:sparsity] < 0.05
      allrej[, as.character(name)] <- mean(rejmat, na.rm = TRUE)
    }
    print(fwer) # FWER
    print(allrej) # power
    resname <- paste0("results ", format(Sys.time(), "%d-%b-%Y %H.%M"),
                      " split=", frac, " B=", B, " seed=", rseed)
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
}

print("Finale")



