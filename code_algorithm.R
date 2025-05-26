#------------
# Implementation of the spatial PCA estimator 
#------------



eval_grad <- function(x, v){
  n <- nrow(x)
  
  x_minus <- -1*sweep(x, 2, v, "-")
  x_plus <- sweep(x, 2, v, "+")
  
  norm_minus <- apply(x_minus, 1, function(z) sqrt(sum(z^2)))
  norm_plus <- apply(x_plus, 1, function(z) sqrt(sum(z^2)))
  
  xv_dist <- apply(cbind(norm_minus, norm_plus), 1, max)
  
  n_set <- (1:n)[which(xv_dist > 1e-6)]
  
  grad_now <- colMeans(sweep(x_minus[n_set, ], 1, norm_plus[n_set]/norm_minus[n_set], "*") + sweep(x_plus[n_set, ], 1, norm_minus[n_set]/norm_plus[n_set], "*"))
  
  grad_now
}


estim_v <- function(x, maxiter = 100){
  
  S <- cov(x)
  
  v0 <- eigen(S)$vectors[, 1]
  grad0 <- eval_grad(x, v0)
  step <- 0.1
  
  v1 <-  v0 - step*grad0
  grad1 <- eval_grad(x, v1)
  
  crit <- 1
  iter <- 1
  
  while(crit > 1e-6 & iter < maxiter){
    step <- abs(sum((v1 - v0)*(grad1 - grad0)))/sum((grad1 - grad0)^2)
    
    v2 <-  v1 - step*grad1
    grad2 <- eval_grad(x, v2)
    
    crit <- min(sqrt(sum((v2 - v1)^2)), sqrt(sum((v2 + v1)^2)))

    grad0 <- grad1
    grad1 <- grad2
    v0 <- v1
    v1 <- v2
    
    iter <- iter + 1
  }
  
  list(v = v1, iter = iter)
}



# n <- 1000
# x <- matrix(rnorm(n*2), n, 2)
# x <- sweep(x, 2, c(2, 1), "*")
# 
# res <- estim_v(x, maxiter = 100)
# res



