rm(list = ls(all = TRUE))
save <- TRUE

# create save location, adjust depending on folder structure

if (save) {
  newdir <- format(Sys.time(), "%d-%b-%Y %H.%M")
  dir.create(paste("multi_carve/",newdir,sep="")) 

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
source("hdi_adjustments.R")
source("optimal_inference.R")
source("tryCatch-W-E.R")

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
# should create the right x on D-MATH server, x[1 ,1] = 0.958


y.true <- x %*% beta
sigma <- 2

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
B <- 50
for (frac in frac.vec) {
  seed.n <- seed.n+1
  set.seed(seed.vec[seed.n])
  # check set-up
  print(frac)
  print(B)
  opts <- list(progress = progress)
  
  # parallelization
  # choose different number of cores if wished
  cl <- makeSOCKcluster(16) 

  rseed <- seed.vec[seed.n]
  clusterSetRNGStream(cl, iseed = rseed) # make things reproducible
  registerDoSNOW(cl)
  tic()
  res<-foreach(gu = 1:nsim, .combine = rbind,
               .packages = c("MASS", "selectiveInference", "glmnet", "Matrix",
                             "hdi", "tmg", "truncnorm"), .options.snow=opts) %dorng%{
   # alternative if sequential computation is preferred
   # res <- foreach(gu = 1:nsim, .combine = rbind) %do%{
    y.true <- x %*% beta
    y <- y.true + sigma * rnorm(n)
    gammavec <- round(seq(ceiling(0.05 * B) / B, 1 - 1 / B, by = 1 / B), 2)
    gammavec[1] <- 0.05999 # due to inconsistency for multisplitting CI
    mcrtry <- tryCatch_W_E(multi.carve.ci.saturated(x, y, B = B, fraction = frac, ci.level = 0.95, 
                                     model.selector = lasso.cvcoef, classical.fit = lm.pval,
                                     parallel = FALSE, ncores = getOption("mc.cores", 2L), gamma = gammavec,
                                     args.model.selector = list(standardize = FALSE, intercept = TRUE,
                                                                tol.beta = 0, use.lambda.min = FALSE),
                                     verbose = FALSE, ci.timeout = 10, FWER = FALSE, split.pval = TRUE,
                                     return.selmodels = TRUE, return.nonaggr = TRUE), 0)
    
    out.list <- list()
    out.list$y <- y
    if(!is.null(mcrtry$error)) {
      err <- paste("mcr:", mcrtry$error)
      war <- if (is.null(mcrtry$warning)) NA
      else c(mcrtry$warning)
      out.list$split.lci <- rep(NA,p)
      out.list$split.uci <- rep(NA,p)
      out.list$carve.lci <- rep(NA,p)
      out.list$carve.uci <- rep(NA,p)
      out.list$beta.min <- rep(NA,p)
      out.list$beta.max <- rep(NA,p)
      out.list$exception <- list(err, paste(1:length(war), ":", war, collapse = ", "))
      out.list
    } else {
      mcr <- mcrtry$value
      
      beta.min <- rep(-Inf, p)
      beta.max <- rep(Inf, p)
      for (j in 1:p) {
        betas <- rep(NA, 50)
        for (i in 1:B) {
          if (mcr[[1]]$sel.models[i, j]) {
            mod <- which(mcr[[1]]$sel.models[i, ])
            coefn <- which(mod == j)
            x.sub <- cbind(rep(1, n), x[, mod])
            betas[i] <- (ginv(x.sub) %*% y.true)[coefn + 1]
          }
        }
        if (any(!is.na(betas))) {
          beta.min[j] <- min(betas, na.rm = TRUE)
          beta.max[j] <- max(betas, na.rm = TRUE)
        }
      }

      out.list$split.lci <- mcr[[2]]$lci # lower end of splitting interval
      out.list$split.uci <- mcr[[2]]$uci # upper end of splitting interval
      out.list$carve.lci <- mcr[[1]]$lci # lower end of carving interval
      out.list$carve.uci <- mcr[[1]]$uci # upper end of carving interval
      out.list$beta.min <- beta.min
      out.list$beta.max <- beta.max
      err <- if (is.null(mcrtry$error)) NA
      else c(mcrtry$error) # should not happen due to earlier check
      war <- if (is.null(mcrtry$warning)) NA
      else c(mcrtry$warning)
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
  expmatr <- matrix(unlist(res[, "exception"]),nrow = dim(res)[1],
                    ncol = 2,byrow = TRUE)
  print(sum(is.na(expmatr[, 1])))
  succ <- which(is.na(expmatr[, 1]))
  print("succesful runs")
  
  all.y <- matrix(unlist(res[,"y"]), nrow = dim(res), byrow = TRUE)
  sd <- attr(res, "rng")
  
  split.lci <- matrix(unlist(res[, "split.lci"]), nrow = dim(res)[1],
                     ncol = p, byrow = TRUE)
  split.uci <- matrix(unlist(res[, "split.uci"]), nrow = dim(res)[1],
                     ncol = p, byrow = TRUE)
  carve.lci <- matrix(unlist(res[, "carve.lci"]), nrow = dim(res)[1],
                     ncol = p, byrow = TRUE)
  carve.uci <- matrix(unlist(res[, "carve.uci"]), nrow = dim(res)[1],
                     ncol = p, byrow = TRUE)
  beta.min <- matrix(unlist(res[, "beta.min"]), nrow = dim(res)[1],
                    ncol = p, byrow = TRUE)
  beta.max <- matrix(unlist(res[, "beta.max"]), nrow = dim(res)[1],
                    ncol = p, byrow = TRUE)
  
  colnames(split.lci) <- paste(1:p)
  colnames(split.uci) <- paste(1:p)
  colnames(carve.lci) <- paste(1:p)
  colnames(carve.uci) <- paste(1:p)
  colnames(beta.min) <- paste(1:p)
  colnames(beta.max) <- paste(1:p)
  split.err <- (split.lci > beta.min | split.uci < beta.max)
  print(mean(apply(split.err, 1, sum))) # average number of false coverage per run for splitting
  carve.err <- (carve.lci > beta.min | carve.uci < beta.max)
  print(mean(apply(carve.err, 1, sum))) # average number of false coverage per run for carving
  split.dis <- split.uci - split.lci
  split.dis[split.err] <- NA
  carve.dis <- carve.uci - carve.lci
  carve.dis[carve.err] <- NA
  print(quantile(split.dis, seq(0.01, 0.1, 0.01), na.rm = TRUE)) # quantile of interval lengths for covering intervals
  print(quantile(carve.dis, seq(0.01, 0.1, 0.01), na.rm = TRUE))
  results <- list("split.lci" = split.lci, "split.uci" = split.uci, "carve.lci" = carve.lci,
                  "carve.uci" = carve.uci, "beta.min" = beta.min, "beta.max" = beta.max)
  
  simulation <- list("results" = results, "exceptions" = expmatr, "y" = all.y, "B" = B,
                     "split" = frac, "nsim" = nsim, "seed" = rseed,
                     "sd" = sd, "commit" = commit)
  resname <- paste0("results ", format(Sys.time(), "%d-%b-%Y %H.%M"), " split=", frac, " B=", B, " seed=", rseed)
  # adjust depending on folder structure
  if (save) save(simulation, file = paste("multi_carve/", newdir, "/", resname, ".RData", sep = ""))
}

      
      