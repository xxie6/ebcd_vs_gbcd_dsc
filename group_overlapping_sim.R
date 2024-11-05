# n
# p
# indiv_sd
# seed

sim_binary_loadings_data <- function(args) {
  set.seed(args$seed)
  
  FF <- matrix(rnorm(args$k * args$p, sd = args$group_sd), ncol = args$k)
  LL <- matrix(rbinom(args$n*args$k, 1, args$pi0), nrow = args$n, ncol = args$k)
  E <- matrix(rnorm(args$n * args$p, sd = args$indiv_sd), nrow = args$n)
  
  Y <- LL %*% t(FF) + E
  YYt <- (1/args$p)*tcrossprod(Y)
  
  return(list(Y = Y, YYt = YYt, LL = LL, FF = FF, K = ncol(LL)))
}

data <- sim_binary_loadings_data(args)