library(ebnm)
library(flashier)
library(ashr)
library(Matrix)
library(stats)
source("gbcd_functions_dsc.R")
gbcd.wrapper <- function(input, args = NULL){
  Kmax <- input$K
  gbcd_fit <- fit_gbcd(Y = input$Y, Kmax = ceiling(Kmax/2))
  return(list(fit_obj = gbcd_fit$res, est_LLt = tcrossprod(gbcd_fit$scaled_L)))
}
gbcd_data <- gbcd.wrapper(input)