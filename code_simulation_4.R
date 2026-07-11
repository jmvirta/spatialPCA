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
  v <- u * sqrt(lambda)
  
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
T_iter <- function(x, v, print=FALSE) {
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
  
  if(print){print(min(sum((v - v_new)^2), sum((v + v_new)^2)))}
  
  return(list(v = v_new, f = mean(np_nm)))
}


# The update function in Algorithm 1 when v is equal to plus minus one of the sample points
#
# x = data as n x p matrix
# v = the current iterate (p-vector)
Tk_iter <- function(x, v, print=FALSE) {
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
  
  if(print){print(min(sum((v - v_new)^2), sum((v + v_new)^2)))}
  
  return(list(v = v_new, f = eval_f(x,v_new)))
}



# The actual optimization function
#
# x = data as n x p matrix
# v = the initial value (p-vector)
# maxiter = maximal number of iterations before manual termination
# maxtol = convergence tolerance of the algorithm 
estim_v <- function(x, v = NULL, maxiter = 100, maxtol = 10^(-6), lambda_optimize = FALSE){
  iter <- 0;
  tol <- 1;
  f_evals <- c();
  
  if(is.null(v)){
    v <- eigen(cov(x))$vector[,1];
    #lambda <- optim(1, eval_f_lambda, x=x, u = v,  method = "L-BFGS-B")$par;
    #v <- lambda*v;
  }
  while (iter < maxiter && tol > maxtol) {
    matching_rows <- apply(x, 1, function(row) all(row == v))
    ifelse(sum(matching_rows) == 0, res <- T_iter( x, v ), res <- Tk_iter( x, v ));
    
    v1 <- res$v;
    f_evals <- c( f_evals, res$f );
    tol <- min(sqrt(sum((v - v1)^2)), sqrt(sum((v + v1)^2)));
    v <- v1;
    iter <- iter + 1;
  }
  ifelse(lambda_optimize == TRUE, lambda <- (optim(sqrt(sum(v^2)), eval_f_lambda, x=x, u = v,  method = "L-BFGS-B")$par)^2, lambda <- 1);
  #v <- lambda*v;
  return(list(v = v, f_evals = f_evals))
}


norm_2 <- function(x){return(norm(x,type="2"))}


#############################################################################################
#--------------------------Sensitivity study-----------------------------------------------
#############################################################################################

# #test1 
# n <- 1000
# p <- 3
# sigma = 1
# x <- sigma*matrix(rnorm(n*p), n, p)%*%diag(c(4,1,1))
# res <- estim_v(x, v = c(1,0,0), maxiter = 100)

mean_sd_init_inner_products <- function(vs) {
  # vs: n_init x p matrix
  G <- abs(vs %*% t(vs))
  list(G=c(G),mu = mean(G[upper.tri(G)]), stdiv = sd(G[upper.tri(G)]))
}

#Function takes fixed seed, dimension p and number of random initializations is inputs
#if normalize = TRUE, focus on unit directions
#returns the mean and sd of the inner
#products between different initializations as estimate of sensitivity


run_one_simu<-function(x, seed=sample(1:10000000,1), p, n_init=25, normalize=TRUE){
  set.seed(seed); #fix seed within one multiple initialization setting
  u_pca <- eigen(cov(x))$vectors[,1]
  us<-rnorm(p*n_init); dim(us)<-c(n_init,p); #us is a matrix with initial approximations in rows 
  us<-rbind(us,u_pca)
  vs<-us; #here collect the solutions;
  for (i in 1:(n_init+1)) {
    res_rand <- estim_v(x=x, v = us[i,]/norm_2(us[i,]), maxiter = 100)$v
    ifelse(normalize, vs[i,] <- res_rand/norm_2(res_rand), vs[i,] <- res_rand)
  }
  return(mean_sd_init_inner_products(vs))
}


replicate_simu <- function(n_rep, n, p, lambda, n_init=25, normalize=TRUE, seed=42) {
  seeds <- seed + seq_len(n_rep)
  
  cl <- parallel::makeCluster(100)
  on.exit(parallel::stopCluster(cl))
  
  parallel::clusterExport(cl, c(
    "run_one_simu", "estim_v", "T_iter", "Tk_iter",
    "eval_f", "eval_f_lambda", "eval_grad",
    "norm_2", "mean_sd_init_inner_products",
    "n", "p", "lambda", "df", "n_init", "normalize"
  ), envir = environment())  # <-- this is the fix
  parallel::clusterEvalQ(cl, library(mvtnorm))
  
  results <- parallel::parLapply(cl, seeds, function(s) {
    set.seed(s)
    #x <- mvtnorm::rmvnorm(n=n, sigma=diag(c(lambda, rep(1, p-1))))
    x <- mvtnorm::rmvt(n=n, sigma=diag(c(lambda, rep(1, p-1))), df=5)
    run_one_simu(x=x, p=p, n_init=n_init, normalize=TRUE)
  })
  
  all_G  <- unlist(lapply(results, `[[`, "G"))
  all_mu <- sapply(results, `[[`, "mu")
  
  list(
    grand_mean = mean(all_G),
    grand_sd   = sd(all_G),
    mean_mu    = mean(all_mu) #sanity check
  )
}


## test
# replicate_simu(n_rep=25,n=1000,p=10,lambda=2,n_init=25)

df <- 5

settings <- expand.grid(
  lambda = c(2, 5, 10, 25, 40),
  n      = c(500, 1000),
  p      = c(3, 10, 50)
  # lambda = c(5, 10),
  # n      = c(50, 100),
  # p      = c(3, 5)
)

all_results <- Map(function(n, p, lambda) {
  replicate_simu(n_rep=100, n=n, p=p, lambda=lambda, n_init=25, normalize=TRUE, seed=42)
}, settings$n, settings$p, settings$lambda)

settings$grand_mean <- sapply(all_results, `[[`, "grand_mean")
settings$grand_sd   <- sapply(all_results, `[[`, "grand_sd")
settings$mean_mu    <- sapply(all_results, `[[`, "mean_mu")

print(settings)

write.table(settings, file="sensitivity_simulation_results.txt", row.names=FALSE, sep="\t", quote=FALSE)



