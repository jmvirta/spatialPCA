#------------
# Implementation of the spatial PCA estimator 
#------------


# Computes the gradient of the objective function
#
# x = data as n x p matrix
# v = the point in R^p where the gradient is evaluated
eval_grad <- function(x, v) {
  # Vectorized subtraction and addition
  x_minus <- sweep(x, 2, v, "-")
  x_plus  <- sweep(x, 2, v, "+")
  
  # Efficient Euclidean norms
  norm_minus <- sqrt(rowSums(x_minus^2))
  norm_plus  <- sqrt(rowSums(x_plus^2))
  
  # Filter small distances
  xv_dist <- pmin(norm_minus, norm_plus)
  n_set <- which(xv_dist > 1e-6)
  
  # Compute ratios only once
  ratio1 <- norm_plus[n_set] / norm_minus[n_set]
  ratio2 <- norm_minus[n_set] / norm_plus[n_set]
  
  x_m <- x_minus[n_set, ]
  x_p <- x_plus[n_set, ]
  
  # Element-wise scaling (vector * matrix row-wise)
  scaled_minus <- x_m * ratio1
  scaled_plus  <- x_p * ratio2
  
  # Compute gradient
  grad_now <- colMeans(scaled_minus + scaled_plus)
  
  return(grad_now)
}


# Evaluates the objective function
#
# x = data as n x p matrix
# v = the point in R^p where the objective function is evaluated
eval_f <- function(x, v) {
  # Vectorized subtraction and addition
  x_minus <- sweep(x, 2, v, "-")
  x_plus  <- sweep(x, 2, v, "+")
  
  # Fast row-wise Euclidean norms
  norm_minus <- sqrt(rowSums(x_minus^2))
  norm_plus  <- sqrt(rowSums(x_plus^2))
  
  # Compute mean of element-wise product
  f_now <- mean(norm_minus * norm_plus)
  
  return(f_now)
}



# Evaluates the objective function in a parametrization v = u*lambda
#
# x = data as n x p matrix
# u = the direction of the point in R^p where the objective function is evaluated
# lambda = the magnitude of the point in R^p where the objective function is evaluated
eval_f_lambda <- function(x, u, lambda) {
  v <- u * lambda
  
  # Vectorized subtraction and addition
  x_minus <- sweep(x, 2, v, "-")
  x_plus  <- sweep(x, 2, v, "+")
  
  # Vectorized row-wise Euclidean norms
  norm_minus <- sqrt(rowSums(x_minus^2))
  norm_plus  <- sqrt(rowSums(x_plus^2))
  
  # Compute mean of element-wise product
  f_now <- mean(norm_minus * norm_plus)
  
  return(f_now)
}


# The update function in Algorithm 1 when v is not equal to plus minus one of the sample points
#
# x = data as n x p matrix
# v = the current iterate (p-vector)
T_iter <- function(x, v) {
  n <- nrow(x)
  
  x_minus <- -1 * sweep(x, 2, v, "-")
  x_plus  <- sweep(x, 2, v, "+")
  
  norm_minus <- sqrt(rowSums(x_minus^2))
  norm_plus  <- sqrt(rowSums(x_plus^2))
  
  xv_dist <- pmin(norm_minus, norm_plus)
  n_set <- which(xv_dist > 1e-6)
  
  nm <- norm_minus[n_set]
  np <- norm_plus[n_set]
  
  np_nm <- np * nm  
  
  c_val <- mean((np^2 + nm^2) / np_nm)
  temp_ci <- (np^2 - nm^2) / np_nm
  
  v_new <- colSums(x[n_set, ] * temp_ci) / (length(n_set) * c_val)
  
  print(min(sum((v - v_new)^2), sum((v + v_new)^2)))
  
  return(list(v = v_new, f = mean(np_nm)))
}


# The update function in Algorithm 1 when v is equal to plus minus one of the sample points
#
# x = data as n x p matrix
# v = the current iterate (p-vector)
Tk_iter <- function(x, v) {
  # Create a logical vector indicating which rows match vector v
  matching_rows <- apply(x, 1, function(row) all(row == v))
  
  # Subset the matrix to remove the matching row(s)
  xk <- x[!matching_rows, , drop = FALSE]
  
  n <- nrow(xk)
  
  
  x_minus <- -1 * sweep(xk, 2, v, "-")
  x_plus  <- sweep(xk, 2, v, "+")
  
  nm <- sqrt(rowSums(x_minus^2))
  np  <- sqrt(rowSums(x_plus^2))
  
  np_nm <- nm * np
  
  c_val <- mean((np^2 + nm^2) / np_nm)
  temp_ci <- (np^2 - nm^2) / np_nm
  
  v_new <- colSums(xk * temp_ci) / (n * c_val)
  
  print(min(sum((v - v_new)^2), sum((v + v_new)^2)))
  
  return(list(v = v_new, f = eval_f(x,v_new)))
}



# The actual optimization function
#
# x = data as n x p matrix
# v = the initial value (p-vector)
# maxiter = maximal number of iterations before manual termination
# maxtol = convergence tolerance of the algorithm 
v_optimize_iterative <- function(x, v = NULL, maxiter = 50, maxtol = 10^(-6)){
  iter <- 0;
  tol <- 1;
  f_evals <- c();
  
  if(is.null(v)){
    v <- eigen(cov(x))$vector[,1];
    lambda <- optim(1, eval_f_lambda, x=x, u = v,  method = "L-BFGS-B")$par;
    v <- lambda*v;
  }
  while (iter < maxiter && tol > maxtol) {
    res <- T_iter( x, v );
    v1 <- res$v;
    f_evals <- c( f_evals, res$f );
    tol <- min(sum((v - v1)^2), sum((v + v1)^2));
    v <- v1;
    iter <- iter + 1;
  }
  lambda <- optim(sqrt(sum(v^2)), eval_f_lambda, x=x, u = v,  method = "L-BFGS-B")$par;
  v <- lambda*v;
  return(list(v = v, f_evals = f_evals))
}



# Example of running the code
#
# n <- 1000
# p <- 100
# sigma = 1
#
# x <- sigma*matrix(rt(n*p, 2), n, p)%*%diag(c(50, rep(1, p-1)))
#
# v0 <- rnorm(p)
#
# system.time({
# res <- v_optimize_iterative(x)
# })
#
# plot(res$f_evals, type="o")
# eval_f(x,res$v)
#
# res$v
# res$v/sqrt(sum(res$v^2))

