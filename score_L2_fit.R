compute_L2_fit <- function(est, dat){
  score <- sum((dat - est)^2) - sum((diag(dat) - diag(est))^2)
  return(score)
}
result <- compute_L2_fit(est,dat)