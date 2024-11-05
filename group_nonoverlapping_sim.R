# n
# pop_sizes
# branch_sds
# indiv_sd
# n_genes
# seed

sim_star_data <- function(args) {
  set.seed(args$seed)
  
  n <- sum(args$pop_sizes)
  p <- args$n_genes
  K <- length(args$pop_sizes)
  
  FF <- matrix(rnorm(K * p, sd = rep(args$branch_sds, each = p)), ncol = K)
  
  LL <- matrix(0, nrow = n, ncol = K)
  for (k in 1:K) {
    vec <- rep(0, K)
    vec[k] <- 1
    LL[, k] <- rep(vec, times = args$pop_sizes)
  }
  
  E <- matrix(rnorm(n * p, sd = args$indiv_sd), nrow = n)
  Y <- LL %*% t(FF) + E
  YYt <- (1/p)*tcrossprod(Y)
  
  return(list(Y = Y, YYt = YYt, LL = LL, FF = FF, K = ncol(LL)))
}

data <- sim_star_data(args)