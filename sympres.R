# SYMPRES

sympres_init <- function(S, 
                          r, 
                          init_method, 
                          seed = 6, 
                          additive_term = TRUE, 
                          off_diagonal = TRUE){
  # initialize P
  n <- nrow(S)
  
  if (init_method == 'default'){
    X <- svd(S, nu = r)
    P <- matrix(as.numeric(X$u > mean(X$u)), nrow = n)
    W <- rep(1, r) # W is a diagonal matrix, just store the diagonal entries
  }
  else if (init_method == 'random'){
    set.seed(seed)
    P <- matrix(round(runif(n*r)), nrow = n)
    W <- runif(r)
  }
  Q <- P
  
  if (additive_term == TRUE){
    if (off_diagonal == TRUE){
      c <- sum((S - P %*% diag(W) %*% t(Q))*(1 - diag(n)))/sum(1 - diag(n))
    }
    else{
      c <- mean(S - P %*% diag(W) %*% t(Q))
    }
  }
  else{
    c <- 0 # C is a matrix where all the values are just one value
  }
  
  obj_func_val <- sum((S - (P%*% diag(W) %*% t(Q)))^2)
  obj_func_list <- c(obj_func_val)
  return(list(P = P, Q = Q, W = W, c = c, 
              obj_func_val = obj_func_val, obj_func_list = obj_func_list, 
              init_method = init_method))
}

sympres_update <- function(sind_obj, S, r, off_diagonal, additive_term, conv_tol){
  n <- nrow(S)
  fold <- sind_obj$obj_func_val + 2*conv_tol*sind_obj$obj_func_val
  num_iter <- 0
  W <- sind_obj$W
  P <- sind_obj$P
  Q <- sind_obj$Q
  c <- sind_obj$c
  
  while (((fold - sind_obj$obj_func_val) > conv_tol*sind_obj$obj_func_val) & (num_iter < 3)){
    fold <- sind_obj$obj_func_val
    num_iter <- num_iter + 1
    for (i in 1:r){
      # define S-bar
      if (r == 1){
        SS <- S - sind_obj$c
      } else{
        SS <- S - P[,-i] %*% diag(W[-c(i)]) %*% t(Q[,-i]) - c
      }
      
      # Update W(i)
      T1 <- c(SS)
      
      # Sub step: Correct P[,i] or Q[,i] in case column is zero; replace it by 1
      while (sum(P[,i]) == 0){
        P[,i] <- round(runif(n))
      }
      Q[,i] <- P[,i]
      
      # Update W(i)
      if (off_diagonal == TRUE){
        #W[i] <- T1/(t(kronecker(Q[,i], P[,i])) * c(1 - diag(n))) #NEED TO CHANGE TO SOLVE LINEAR SYSTEM; This is matlab notation
        g <- (t(kronecker(Q[,i], P[,i])) * c(1 - diag(n)))
        W[i] <- sum((T1)*(g))/sum((g)^2)
      }
      else{
        #W[i] <- T1/(t(kronecker(Q[,i], P[,i])))
        g <- (t(kronecker(Q[,i], P[,i])))
        W[i] <- sum((T1)*(g))/sum((g)^2)
      }
      # Set negative elements to zero
      if (W[i] < 0){
        W[i] <- 0
      }
      
      # Update P[,i]
      if (W[i] == 0){
        P[,i] <- round(runif(n))
      }
      else{
        if (off_diagonal == TRUE){
          G <- (W[i]^2) * (matrix(rep(1, n^2), nrow = n) - diag(n))
          G <- G - 2 * W[i] * (SS * (matrix(rep(1, n^2), nrow = n) - diag(n)))
        }
        else{
          G <- (W[i]^2)*matrix(rep(1, n^2), nrow = n)
          G <- G - 2*W[i] * SS
        }
        for (j in 1:n){
          P[j,i] <- 0.5
          P[j,i] <- as.numeric((G[j,] %*% P[,i])  < 0)
        }
      }
      Q[,i] <- P[,i]
    }
    
    # Update c
    if (additive_term == TRUE){
      if (off_diagonal == TRUE){
        c <- sum((S - P %*% diag(W) %*% t(Q))*(1-diag(n)))/sum(1 - diag(n))
      }
      else{
        c <- mean(S - P %*% diag(W) %*% t(Q))
      }
    }
    
    if (off_diagonal == TRUE){
      sind_obj$obj_func_val <- sum(((S - (P%*% diag(W) %*% t(Q)) - c) * (1 - diag(n)))^2)
    }
    else{
      sind_obj$obj_func_val <- sum((S - (P%*% diag(W) %*% t(Q)) - c)^2)
    }
    sind_obj$obj_func_list <- c(sind_obj$obj_func_list, sind_obj$obj_func_val)
  }
  sind_obj$W <- W
  sind_obj$P <- P
  sind_obj$Q <- Q
  sind_obj$c <- c
  sind_obj$off_diagonal <- off_diagonal
  sind_obj$additive_term <- additive_term
  return(sind_obj)
}

fit_sympres <- function(S, 
                         r, 
                         init_method, 
                         conv_tol = 1e-6,
                         seed = 6, 
                         additive_term = TRUE, 
                         off_diagonal = TRUE){
  sym_obj_init <- sympres_init(S, r, init_method, seed,
                                 additive_term, off_diagonal)
  sym_obj <- sympres_update(sym_obj_init, S, r, off_diagonal, additive_term, 
                              conv_tol)
  return(sym_obj)
}



