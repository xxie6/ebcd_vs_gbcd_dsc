library(ebcd)
ebcd.wrapper <- function(input, args=NULL){
  set.seed(args$seed)
  Kmax <- input$K
  ebcd_fit <- ebcd(X = t(input$Y), Kmax = Kmax, maxiter_backfit = 10000, ebnm_fn = args$ebnm_fn)
  return(list(ebcd_fit = ebcd_fit, est_LLt = tcrossprod(ebcd_fit$EL)))
}
ebcd_data <- ebcd.wrapper(input, args)