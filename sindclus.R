# SINDCLUS

sindclus_init <- function(S, 
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
              obj_func_val = obj_func_val, obj_func_list = obj_func_list, init_method = init_method))
}

sindclus_update <- function(sind_obj, S, r, off_diagonal, additive_term, conv_tol){
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
      T1 <- t(c(SS))
      
      # Sub step: Correct P[,i] or Q[,i] in case column is zero; replace it by 1
      while (sum(P[,i]) == 0){
        P[,i] <- round(runif(n))
      }
      
      while (sum(Q[,i]) == 0){
        Q[,i] <- round(runif(n))
      }
      
      # Update W(i)
      if (off_diagonal == TRUE){
        #W[i] <- T1/(t(kronecker(Q[,i], P[,i])) * c(1 - diag(n))) #This is matlab notation
        g <- (t(kronecker(Q[,i], P[,i])) * c(1 - diag(n)))
        W[i] <- sum((T1)*(g))/sum((g)^2)
      }
      else{
        #W[i] <- T1/(t(kronecker(Q[,i], P[,i]))) #This is matlab notation
        g <- (t(kronecker(Q[,i], P[,i])))
        W[i] <- sum((T1)*(g))/sum((g)^2)
      }
      # Set negative elements to zero
      if (W[i] < 0){
        W[i] <- 0
      }
      
      if (off_diagonal == TRUE){
        # Update P (original SINDCLUS algo)
        V1 <- rowSums((SS - t(t(rep(1,n)))%*%kronecker(W[i], Q[,i]) )^2 * c(1-diag(n)))
        V2 <- rowSums((SS^2)* c(1-diag(n))) 
        
        P[,i] <- as.numeric(V1 < V2)
        # Update Q (original SINDCLUS algo)
        T3 <- t(SS)
        X1 <- rowSums((T3 - t(t(rep(1,n)))%*%kronecker(W[i], P[,i]) )^2 * c(1-diag(n)))
        X2 <- rowSums((T3^2)* c(1-diag(n)))
        Q[,i] <- as.numeric(X1 < X2)
      }
      else{
        # Update P (original SINDCLUS algo)
        V1 <- rowSums((SS - t(t(rep(1,n)))%*%kronecker(W[i], Q[,i]) )^2)
        V2 <- rowSums((SS^2)) #check
        
        P[,i] <- as.numeric(V1 < V2)
        # Update Q (original SINDCLUS algo)
        T3 <- t(SS)
        X1 <- rowSums((T3 - t(t(rep(1,n)))%*%kronecker(W[i], P[,i]) )^2)
        X2 <- rowSums((T3^2))
        Q[,i] <- as.numeric(X1 < X2)
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

fit_sindclus <- function(S, 
                         r, 
                         init_method, 
                         conv_tol = 1e-6,
                         seed = 6, 
                         additive_term = TRUE, 
                         off_diagonal = TRUE){
  sind_obj_init <- sindclus_init(S, r, init_method, seed,
                                 additive_term, off_diagonal)
  sind_obj <- sindclus_update(sind_obj_init, S, r, off_diagonal, additive_term, 
                              conv_tol)
  return(sind_obj)
}
