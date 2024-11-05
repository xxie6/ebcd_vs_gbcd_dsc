source("sympres.R")
sympres.wrapper <- function(input, args = NULL){
  Kmax <- input$K
  sympres_fit <- fit_sympres(input$YYt, 
                               Kmax, 
                               args$init_method, 
                               conv_tol = 1e-6,
                               seed = 6, 
                               args$additive_term, 
                               args$off_diagonal)
  return(list(sympres_fit = sympres_fit, 
              est_LLt = (sympres_fit$P %*% diag(sympres_fit$W) %*% t(sympres_fit$Q))))
}
sympres_data <- sympres.wrapper(input, args)