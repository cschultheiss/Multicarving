rm(list = ls(all = TRUE))
local <- FALSE
save <- FALSE
if (!save && !local) {
  save <- TRUE
  warning("Do not use cluster without saving, save was set to TRUE")
}
# create save location, adjust depending on folder structure
if (save) {
  newdir <- format(Sys.time(), "%d-%b-%Y %H.%M")
  if (local) {
    dir.create(paste("C:/Users/Christoph/Documents/ETH/MA/R_trials/group_example/", newdir, sep = ""))
  } else {
    dir.create(paste("group_example/",newdir,sep="")) 
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
  source("sample_from_truncated.R")
  source("tryCatch-W-E.R")
}

if (local) {
  library(tcltk)
}

p <- 500
rho <- 0.6
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
# sparse alternative
subcov <- matrix(0.8, 5, 5)
diag(subcov) <- 1
Cov[1:5, 1:5] <- subcov

truebeta <- rep(0, p)
# delta_vec <- seq(0, 0.06, 0.02)
# n_vec <- c(250, 350, 500, 800)
# sparse alternative
delta_vec <- seq(0, 0.5, 0.1)
n_vec <- c(250, 350, 500)
frac_vec <- c(0.5, 0.75, 0.9, 0.95, 0.99, 1)
B_vec <- c(1, (1:5) * 10)
b_vec <- B_vec
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
seed_v <- sample(1:10000,length(frac_vec)*length(delta_vec)*length(n_vec))
print(seed_v) # 3588 3052 2252 5257 8307 ...
seed_n <- 0

sigma <- 1
# ind <- c(25:50)
# groupstotest <- list(30:200)
# sparse alternative
ind <- c(1, 3)
groupstotest <- list(c(1:5))

for (n in n_vec) {
  for (delta in delta_vec) {
    truebeta[ind] <- delta
    for (frac in frac_vec) {
      B_vec <- b_vec
      if (frac == 1) B_vec <- 1 # can not do multicarving with f = 1
      seed_n <- seed_n+1
      set.seed(seed_v[seed_n])
      B <- max(B_vec)
      print(delta)
      print(n)
      print(frac)
      print(B)
      if (local) {
        pb <- tkProgressBar(max = ntasks, title = paste("split=", frac))
      }
      opts <- list(progress = progress)
      
      if (local) {
        cl<-makeSOCKcluster(4)
      } else {
        cl<-makeSOCKcluster(16) 
      }
      rseed <- seed_v[seed_n]
      clusterSetRNGStream(cl, iseed = rseed) # make things reproducible
      registerDoSNOW(cl)
      tic()
      res<-foreach(gu=1:nsim, .combine = rbind,
                   .packages = c("MASS","selectiveInference","glmnet","Matrix",
                                 "hdi", "tmg","truncnorm"),.options.snow=opts) %dorng%{
      # alternative if sequential computation is preferred
      # res<-foreach(gu=1:nsim,.combine = rbind) %do%{
                                   
        x <- mvrnorm(n, rep(0, p), Cov)
        u <- rnorm(n) * sigma
        y <- x %*% truebeta + u                          
        sigmahat <- estimateSigma(scale(x, T, F), scale(y, T, F),
                                  intercept = FALSE, standardize = FALSE)$sigmahat
        mcgtry <- tryCatch_W_E(multi.carve_group(x, y, B, frac, model.selector = lasso.cvcoef, gamma = 1,
                                                 return.nonaggr = TRUE, return.selmodels = TRUE, skip.groups = FALSE,
                                                 args.model.selector = list(standardize = FALSE, intercept = TRUE,
                                                                            tol.beta = 0, use_lambda.min = TRUE),
                                                 args.lasso.inference = list(verbose = TRUE, sigma = sigmahat),
                                                 groups = groupstotest), 0)
        out_list <- list()
        out_list$y <- y
        if (!is.null(mcgtry$error)) {
          # error handling
          err <- paste("mcr:", mcgtry$error)
          war <- if (is.null(mcgtry$warning)) NA
          else c(mcgtry$warning)
          for (b in B_vec) {
            if (b > 1) {
              out_list[[as.character(b)]] <- rep(NA, (length(groupstotest)) * 4)
            } else {
              out_list[[as.character(b)]] <- rep(NA, (length(groupstotest)) * 2)
            }
          }
          out_list$exception <- list(err, paste(1:length(war), ":", war, collapse = ", "))
          out_list
        } else {
          mcg <- mcgtry$value
          pcarve <- mcg$pvals.nonaggr
          for (B in B_vec) {
            if (B > 1) {
              use <- 1:B
              pvals.aggregated <- pval.aggregator(list(matrix(pcarve[use,],nrow = B)),
                                                  round(seq(ceiling(0.05 * B)/B, 1 - 1/B, by = 1/B), 2))
              pvals.aggregated2 <- pval.aggregator(list(matrix(pcarve[use,],nrow = B)),
                                                 round(seq(ceiling(0.3 * B)/B, 1 - 1/B, by = 1/B), 2))
              pvals.aggregated3 <- pval.aggregator(list(matrix(pcarve[use,],nrow = B)),
                                                 round(ceiling(0.05 * B)/B, 2))
              pvals.aggregated4 <- pval.aggregator(list(matrix(pcarve[use,],nrow = B)),
                                                 round(ceiling(0.3 * B)/B ,2))
            } else {
              pvals.aggregated <- list(pcarve[1, ])
            }
            run_res <- vector(length = 0)
            np <- length(pvals.aggregated)
            for (i in 1:np) {
              pv <- pvals.aggregated[[i]]
              run_res <- c(run_res, pv)
            }
            if (B > 1) {
              for (i in 1:np) {
                pv <- pvals.aggregated2[[i]]
                run_res <- c(run_res, pv)
              }
              for (i in 1:np) {
                pv <- pvals.aggregated3[[i]]
                run_res <- c(run_res, pv)
              }
              for (i in 1:np) {
                pv <- pvals.aggregated4[[i]]
                run_res <- c(run_res, pv)
              }
            }
            if (B == 1) {
              for (i in length(groupstotest)) {
                tested <- sum((mcg$sel.models[1, groupstotest[[i]]]))
                run_res <- c(run_res, tested)
              }
            }
            out_list[[as.character(B)]] <- run_res
          }
          err <- if (is.null(mcgtry$error)) NA
          else c(mcgtry$error) # should not happen due to earlier check
          war <- if (is.null(mcgtry$warning)) NA
          else c(mcgtry$warning)
          out_list$exception <- list(err, paste(1:length(war), ":", war, collapse = ", "))
          out_list
        }
        # end of simulation run
      }
      toc()
      stopCluster(cl)
      if (local) {
        close(pb)
      }
      expmatr <- matrix(unlist(res[, "exception"]),nrow = dim(res)[1],
                        ncol = 2, byrow = TRUE)
      print(sum(is.na(expmatr[, 1])))
      succ = which(is.na(expmatr[, 1]))
      print("succesful runs")
      all_y <- matrix(unlist(res[,"y"]), nrow = dim(res), byrow = TRUE)
      sd <- attr(res, "rng")
      for (B in B_vec) {
        if (B == 1) {
          names <- c("carve")
          subres <- matrix(unlist(res[, as.character(B)]), nrow = dim(res)[1],
                           ncol = length(groupstotest) * 2, byrow = TRUE)
          if (any(!is.na(subres[-succ, ]))) print("not as it should be")
          subres <- subres[succ, ]
          colnames(subres) <- c(rep(names, each = (length(groupstotest))), "|G|")
        } else {
          names <- c("carve5", "carve30", "carvefix5", "carvefix30")
          subres <- matrix(unlist(res[, as.character(B)]), nrow = dim(res)[1],
                           ncol = 4 * (length(groupstotest)), byrow = TRUE)
          if (any(!is.na(subres[-succ, ]))) print("not as it should be")
          subres = subres[succ, ]
          colnames(subres) <- c(rep(names, each = (length(groupstotest))))
        }
        subres <- as.data.frame(subres)
        simulation <- list("results" = subres, "B" = B, "n" = n ,
                           "exceptions" = expmatr, "y" = all_y, "split" = frac,
                           "delta" = delta, "nsim" = nsim, "seed" = rseed, "sd" = sd, "commit" = commit )
        print(paste("results using fraction ", frac, " B=", B, " delta=",
                    delta, " and n=", n, sep=""))
        options(digits = 3)
        if (B == 1) {
          # probability that group is tested for
          print(apply(matrix(subres[ ,(length(groupstotest) + 1):( 2 *length(groupstotest))],
                             ncol = length(groupstotest)) > 0, 2, mean))
          # average number of variables in group selected
          print(apply(matrix(subres[, (length(groupstotest) + 1):( 2 *length(groupstotest))],
                             ncol = length(groupstotest)), 2, mean))
        }
        rej <- numeric(length(names) * length(groupstotest))
        names(rej) <- rep(names, length(groupstotest))
        for (name in names) {
          nameind <- which(colnames(subres) == name)
          mat <- matrix(subres[, nameind], ncol = length(groupstotest))
          rej[nameind] <- apply(mat < 0.05, 2, mean)
        }
        print(rej)
        resname <- paste0("results ", format(Sys.time(), "%d-%b-%Y %H.%M")," n=",n,
                          " split=", frac, " B=", B, " delta=", delta, " seed=", rseed)
        # adjust depending on folder structure
        if (save) {
          if (local) {
            save(simulation, file = paste("C:/Users/Christoph/Documents/ETH/MA/R_trials/group_example/",
                                          newdir, "/", resname, ".RData", sep = ""))
          } else {
            save(simulation, file = paste("group_example/", newdir, "/", resname, ".RData", sep = ""))
          }
        }
      }
      # end of simulation for given fraction
    }
    # end of simulation for given delta
  }
  # end of simulation for given n
}



