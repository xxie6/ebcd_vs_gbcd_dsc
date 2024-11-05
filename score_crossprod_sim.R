# the metric in essence is the average cosine similarity
library(gtools)

compute_crossprod_single <- function(perm, est, truth, norms_est, norms_truth){
  crossprod_vals <- abs(diag(crossprod(est[,perm], truth)))
  crossprod_vals <- crossprod_vals/(norms_est[perm]*norms_truth)
  return(crossprod_vals)
}

compute_crossprod_similarity <- function(est, truth){
  K_est <- ncol(est)
  K_truth <- ncol(truth)
  n <- nrow(est)
  
  #if estimates don't have same number of columns, try padding the estimate with zeros and make cosine similarity zero
  if (K_est < K_truth){
    est <- cbind(est, matrix(rep(0, n*(K_truth-K_est)), nrow = n))
  }
  
  norms_est <- sqrt(colSums(est^2))
  norms_truth <- sqrt(colSums(truth^2))
  norms_est[norms_est == 0] <- Inf
  
  if (K_est <= K_truth){
    all_perms <- permutations(K_truth, K_truth, repeats.allowed = FALSE)
  }
  if (K_est > K_truth){
    all_perms <- permutations(K_est, K_truth, repeats.allowed = FALSE)
  }
  
  perm_crossprods <- apply(all_perms, 1, compute_crossprod_single, est = est, truth = truth, norms_est = norms_est, norms_truth = norms_truth)
  which.best.perm <- which.max(colSums(perm_crossprods))
  
  return(mean(perm_crossprods[,which.best.perm]))
}

result <- compute_crossprod_similarity(est, truth)