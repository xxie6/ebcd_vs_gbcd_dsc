library(flashier)
library(dplyr)
flashcov.wrapper <- function(input, args){
  Kmax <- input$K
  flash_cov_fit <- flash_init(data = input$YYt, var_type = 0) %>%
    flash_greedy(ebnm_fn = args$ebnm_fn, Kmax = Kmax) %>%
    flash_backfit()
  flash_cov_ldf <- ldf(flash_cov_fit)
  flash_cov_L <- flash_cov_ldf$L %*% diag(sqrt(flash_cov_ldf$D))
  return(list(flash_cov_fit = flash_cov_fit, est_LLt = tcrossprod(flash_cov_L)))
}
flash_cov_data <- flashcov.wrapper(input, args)