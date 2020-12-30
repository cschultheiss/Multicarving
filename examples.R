require(MASS)
require(glmnet)
require(Matrix)
require(tictoc)
require(hdi)
require(selectiveInference)
require(tmg)

# adjust depending on folder structure
source("inference/hdi_adjustments.R")
source("inference/optimal_inference.R")
source("inference/sample_from_truncated.R")
source("inference/tryCatch-W-E.R")

n <- 100
p <- 200
rho <- 0.6
Cov <- toeplitz(rho ^ (seq(0, p - 1)))
sel.index <- c(1, 5, 10, 15, 20)
ind <- sel.index
beta <- rep(0, p)
beta[sel.index] <- 1
sparsity <- 5
RNGkind("Mersenne-Twister")
set.seed(42)
x <- mvrnorm(n, rep(0, p), Cov)
print (x[1,1])
# should create the right x, x[1 ,1] = 0.958 for Toeplitz 0.6
y.true <- x %*% beta
sigma <- 2
y <- y.true + sigma * rnorm(n)

set.seed(12)
mc <- multi.carve(x, y)
mc

c100 <- carve100(x, y)
c100

set.seed(12)
mc.ci <- multi.carve.ci.saturated(x, y)
mc.ci

set.seed(12)
groups <- list(c(1:5), c(21:30), c(50:200))
mc.g <- multi.carve.group(x, y, groups)
mc.g
