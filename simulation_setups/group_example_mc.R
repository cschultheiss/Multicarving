rm(list = ls(all = TRUE))
save <- TRUE

# create save location, adjust depending on folder structure
if (save) {
  newdir <- format(Sys.time(), "%d-%b-%Y %H.%M")
  dir.create(paste("simulation_setups/group_example/", newdir, sep="")) 
  
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

source("hdi_adjustments.R")
source("optimal_inference.R")
source("sample_from_truncated.R")
source("tryCatch-W-E.R")

p <- 500
rho <- 0.6
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
# sparse alternative
subcov <- matrix(0.8, 5, 5)
diag(subcov) <- 1
Cov[1:5, 1:5] <- subcov

truebeta <- rep(0, p)
# delta.vec <- seq(0, 0.06, 0.02)
# n.vec <- c(250, 350, 500, 800)
# sparse alternative
delta.vec <- seq(0, 0.5, 0.1)
n.vec <- c(250, 350, 500)
frac.vec <- c(0.5, 0.75, 0.9, 0.95, 0.99, 1)
B.vec <- c(1, (1:5) * 10)
b.vec <- B.vec
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
seed.vec <- sample(1:10000, length(frac.vec) * length(delta.vec) * length(n.vec))
print(seed.vec) # 3588 3052 2252 5257 8307 ...
seed.n <- 0

sigma <- 1
# ind <- c(25:50)
# groups.to.test <- list(30:200)
# sparse alternative
ind <- c(1, 3)
groups.to.test <- list(c(1:5))

for (n in n.vec) {
  for (delta in delta.vec) {
    truebeta[ind] <- delta
    for (frac in frac.vec) {
      B.vec <- b.vec
      if (frac == 1) B.vec <- 1 # can not do multicarving with f = 1
      seed.n <- seed.n+1
      set.seed(seed.vec[seed.n])
      B <- max(B.vec)
      print(delta)
      print(n)
      print(frac)
      print(B)
      opts <- list(progress = progress)
      

      # parallelization
      # choose different number of cores if wished
      cl<-makeSOCKcluster(16) 

      rseed <- seed.vec[seed.n]
      clusterSetRNGStream(cl, iseed = rseed) # make things reproducible
      registerDoSNOW(cl)
      tic()
      res<-foreach(gu = 1:nsim, .combine = rbind,
                   .packages = c("MASS", "selectiveInference", "glmnet", "Matrix",
                                 "hdi", "tmg", "truncnorm", "tictoc"), .options.snow=opts) %dorng%{
      # alternative if sequential computation is preferred
      # res <- foreach(gu = 1:nsim, .combine = rbind) %do%{
                                   
        x <- mvrnorm(n, rep(0, p), Cov)
        u <- rnorm(n) * sigma
        y <- x %*% truebeta + u                          
        sigmahat <- estimateSigma.flex(x, y,
                                  intercept = TRUE, standardize = FALSE, use.lambda.min = TRUE, df.corr = TRUE)$sigmahat
        mcgtry <- tryCatch_W_E(multi.carve.group(x, y, B, frac, model.selector = lasso.cvcoef,
                                                 return.nonaggr = TRUE, return.selmodels = TRUE, skip.groups = FALSE,
                                                 args.model.selector = list(standardize = FALSE, intercept = TRUE,
                                                                            tol.beta = 0, use.lambda.min = TRUE),
                                                 args.lasso.inference = list(verbose = TRUE, sigma = sigmahat),
                                                 groups = groups.to.test), 0)
        out.list <- list()
        out.list$y <- y
        if (!is.null(mcgtry$error)) {
          # error handling
          err <- paste("mcr:", mcgtry$error)
          war <- if (is.null(mcgtry$warning)) NA
          else c(mcgtry$warning)
          for (b in B.vec) {
            if (b > 1) {
              out.list[[as.character(b)]] <- rep(NA, (length(groups.to.test)) * 4)
            } else {
              out.list[[as.character(b)]] <- rep(NA, (length(groups.to.test)) * 2)
            }
          }
          out.list$exception <- list(err, paste(1:length(war), ":", war, collapse = ", "))
          out.list
        } else {
          mcg <- mcgtry$value
          pcarve <- mcg$pvals.nonaggr
          for (B in B.vec) {
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
            run.res <- vector(length = 0)
            np <- length(pvals.aggregated)
            for (i in 1:np) {
              pv <- pvals.aggregated[[i]]
              run.res <- c(run.res, pv)
            }
            if (B > 1) {
              # for multicarving, test different aggregation methods
              for (i in 1:np) {
                pv <- pvals.aggregated2[[i]]
                run.res <- c(run.res, pv)
              }
              for (i in 1:np) {
                pv <- pvals.aggregated3[[i]]
                run.res <- c(run.res, pv)
              }
              for (i in 1:np) {
                pv <- pvals.aggregated4[[i]]
                run.res <- c(run.res, pv)
              }
            }
            if (B == 1) {
              for (i in length(groups.to.test)) {
                # analyse first split specially for B = 1
                tested <- sum((mcg$sel.models[1, groups.to.test[[i]]])) # number of variables tested from group
                run.res <- c(run.res, tested)
              }
            }
            out.list[[as.character(B)]] <- run.res
          }
          err <- if (is.null(mcgtry$error)) NA
          else c(mcgtry$error) # should not happen due to earlier check
          war <- if (is.null(mcgtry$warning)) NA
          else c(mcgtry$warning)
          out.list$exception <- list(err, paste(1:length(war), ":", war, collapse = ", "))
          out.list
        }
        # end of simulation run
      }
      toc()
      stopCluster(cl)

      # get matrix of errors and warnings
      expmatr <- matrix(unlist(res[, "exception"]),nrow = dim(res)[1],
                        ncol = 2, byrow = TRUE)
      print(sum(is.na(expmatr[, 1])))
      succ = which(is.na(expmatr[, 1]))
      print("succesful runs")
      all.y <- matrix(unlist(res[,"y"]), nrow = dim(res), byrow = TRUE)
      sd <- attr(res, "rng")
      for (B in B.vec) {
        if (B == 1) {
          names <- c("carve")
          subres <- matrix(unlist(res[, as.character(B)]), nrow = dim(res)[1],
                           ncol = length(groups.to.test) * 2, byrow = TRUE)
          if (any(!is.na(subres[-succ, ]))) print("not as it should be")
          subres <- subres[succ, ]
          colnames(subres) <- c(rep(names, each = (length(groups.to.test))), "|G|")
        } else {
          names <- c("carve5", "carve30", "carvefix5", "carvefix30")
          subres <- matrix(unlist(res[, as.character(B)]), nrow = dim(res)[1],
                           ncol = 4 * (length(groups.to.test)), byrow = TRUE)
          if (any(!is.na(subres[-succ, ]))) print("not as it should be")
          subres = subres[succ, ]
          colnames(subres) <- c(rep(names, each = (length(groups.to.test))))
        }
        subres <- as.data.frame(subres)
        simulation <- list("results" = subres, "B" = B, "n" = n ,
                           "exceptions" = expmatr, "y" = all.y, "split" = frac,
                           "delta" = delta, "nsim" = nsim, "seed" = rseed, "sd" = sd, "commit" = commit )
        print(paste("results using fraction ", frac, " B=", B, " delta=",
                    delta, " and n=", n, sep=""))
        options(digits = 3)
        if (B == 1) {
          # probability that group is tested for
          print(apply(matrix(subres[ ,(length(groups.to.test) + 1):( 2 *length(groups.to.test))],
                             ncol = length(groups.to.test)) > 0, 2, mean))
          # average number of variables in group selected
          print(apply(matrix(subres[, (length(groups.to.test) + 1):( 2 *length(groups.to.test))],
                             ncol = length(groups.to.test)), 2, mean))
        }
        rej <- numeric(length(names) * length(groups.to.test))
        names(rej) <- rep(names, length(groups.to.test))
        for (name in names) {
          nameind <- which(colnames(subres) == name)
          mat <- matrix(subres[, nameind], ncol = length(groups.to.test))
          rej[nameind] <- apply(mat < 0.05, 2, mean)
        }
        print(rej)
        resname <- paste0("results ", format(Sys.time(), "%d-%b-%Y %H.%M")," n=",n,
                          " split=", frac, " B=", B, " delta=", delta, " seed=", rseed)
        # adjust depending on folder structure
        if (save) save(simulation, file = paste("simulation_setups/group_example/", newdir, "/", resname, ".RData", sep = ""))

      }
      # end of simulation for given fraction
    }
    # end of simulation for given delta
  }
  # end of simulation for given n
}



