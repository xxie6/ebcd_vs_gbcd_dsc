# pop_sizes
# branch_sds
# indiv_sd
# n_genes
# constrain_F
# seed


sim_4pops <- function(args) {
  set.seed(args$seed)
  
  n <- sum(args$pop_sizes)
  p <- args$n_genes
  
  FF <- matrix(rnorm(7 * p, sd = rep(args$branch_sds, each = p)), ncol = 7)
  # if (args$constrain_F) {
  #   FF_svd <- svd(FF)
  #   FF <- FF_svd$u
  #   FF <- t(t(FF) * branch_sds * sqrt(p))
  # }
  
  LL <- matrix(0, nrow = n, ncol = 7)
  LL[, 1] <- 1
  LL[, 2] <- rep(c(1, 1, 0, 0), times = args$pop_sizes)
  LL[, 3] <- rep(c(0, 0, 1, 1), times = args$pop_sizes)
  LL[, 4] <- rep(c(1, 0, 0, 0), times = args$pop_sizes)
  LL[, 5] <- rep(c(0, 1, 0, 0), times = args$pop_sizes)
  LL[, 6] <- rep(c(0, 0, 1, 0), times = args$pop_sizes)
  LL[, 7] <- rep(c(0, 0, 0, 1), times = args$pop_sizes)
  
  E <- matrix(rnorm(n * p, sd = args$indiv_sd), nrow = n)
  Y <- LL %*% t(FF) + E
  YYt <- (1/p)*tcrossprod(Y)
  return(list(Y = Y, YYt = YYt, LL = LL, FF = FF, K = ncol(LL)))
}

data <- sim_4pops(args)