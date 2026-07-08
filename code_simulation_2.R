#-------------
# Simulation study #2
#-------------


library(tidyverse)


#
# Run the methods from code_algorithm.R
#


norm_2 <- function(x){return(norm(x,type="2"))}



n <- 1000
p <- 3
sigma = 1
x <- sigma*matrix(rnorm(n*3), n, 3)%*%diag(c(4,1,1))
res <- estim_v(x, v = c(1,0,0), maxiter = 100)

#mean(sqrt(apply((x-matrix(rep(res$v,n),byrow=T,ncol=3) )^2,1,sum))/sqrt(apply((x+matrix(rep(res$v,n),byrow=T,ncol=3))^2,1,sum)))


#x <- sweep(x, 1, c(2, 1), "*")

#res <- estim_v(x, maxiter = 100)
#res

(v <- res$v)

#Sigma_eigens <- eigen(cov(x))$values
#sigma1 <- Sigma_eigens[1]
#sigma2 <- Sigma_eigens[2]

sigma1 <- (sigma*4)^2;
sigma2 <- sigma^2;

psi <- norm(res$v,"2")

o1 <- c(1,rep(0,p-1))
o2 <- c(0,1,rep(0,p-2))

v_matrix <- matrix(rep(psi*o1,n),nrow=n,byrow=TRUE )

x_minus <- -1 * (x-v_matrix)
x_plus <- x+v_matrix

y_minus <- apply(x_minus,1,norm_2)
y_plus <- apply(x_plus,1,norm_2)


##Hessian part

a1 <- mean( y_minus/y_plus+y_plus/y_minus )

a2 <- mean( y_plus / y_minus^3 * (tcrossprod(o1, x_minus) )^2 )

a3 <- mean( y_minus / y_plus^3 * (tcrossprod(o1, x_plus) )^2 )

a4 <- mean( 1 / (y_minus * y_plus ) * tcrossprod(o1, x_plus) * tcrossprod(o1, x_minus) )


b1 <- a1

b2 <- mean( y_plus / y_minus^3 * (tcrossprod(o2, x_minus) )^2 )

b3 <- mean( y_minus / y_plus^3 * (tcrossprod(o2, x_plus) )^2 )

b4 <- mean( 1 / (y_minus * y_plus ) * tcrossprod(o2, x_plus) * tcrossprod(o2, x_minus) )

a_hessian <- (a1 - ( a2 + a3 ) + 2*a4)
b_hessian <- (b1 - ( b2 + b3 ) + 2*b4)


## Gradient part


c11 <- mean( y_plus^2 / y_minus^2 * (tcrossprod(o1, x_minus) )^2 )
c12 <- mean( y_minus^2 / y_plus^2 * (tcrossprod(o1, x_plus) )^2 )

c1 <- ( c11 + c12 )

d11 <- mean( y_plus^2 / y_minus^2 * (tcrossprod(o2, x_minus) )^2 )
d12 <- mean( y_minus^2 / y_plus^2 * (tcrossprod(o2, x_plus) )^2 )

d1 <- ( d11 + d12 )

c_grad <- c1 + 2*psi^2 - 2*sigma1
d_grad <- d1 - 2*sigma2

c <- c_grad / a_hessian^2
d <- d_grad / b_hessian^2

limit_covariance <- c * tcrossprod(o1,o1) + d * ( diag(p) - tcrossprod(o1,o1))

limit_covariance_unit_v <-  d / psi^2 * ( diag(p) - tcrossprod(o1,o1))


#############
############# Testing if it is true

library(MASS)  # For matrix operations (if needed)
library(doParallel)
library(foreach)
library(abind)
# # Define parameters
# n <- 1000
# p <- 3
# 
# num_cores <- 100  # Number of cores
# reps_per_core <- 50  # Number of repetitions per core
# total_reps <- num_cores * reps_per_core  # Total replications


# Define function for one replication
run_simulation <- function() {
  x <- matrix(rnorm(n * p), n, p) %*% diag(c(4, 1, 1))
  res <- estim_v(x, maxiter = 500)
  
  # Sign correction
  ind <- which.min(c(sum((res$v / sqrt(sum((res$v)^2)) - c(1, 0, 0))^2),
                     sum((res$v / sqrt(sum((res$v)^2)) + c(1, 0, 0))^2)))
  
  # Compute the correct residual
  resid <- res$v + (-1)^ind * psi* c(1, 0, 0)
  
  # Compute final term
  final_stat <- sqrt(n) * resid
  
  return(final_stat)  # Return result for one replication
}

# Define function for one replication
run_simulation_unit_v <- function(n,p,lambda,df,seed = set.seed()) {
  x <- mvtnorm::rmvt(n = n, sigma = diag(c(lambda, rep(1,p-1) ) ), df = df )# * sqrt( (df-2)/df );
  
  res <- estim_v(x, maxiter = 500)
  ifelse(norm_2(res$v)>0, v <- res$v / norm_2(res$v), v <- res$v )
  
  # Sign correction
  resid_plus <- v + c(1, rep(0,p-1));
  resid_minus <- v - c(1, rep(0,p-1));
  
  ifelse(norm_2(resid_plus) < norm_2(resid_minus), resid_spat <- resid_plus, resid_spat <- resid_minus)
  
  ##PCA
  
  Sig <- cov(x)
  v_pca <- eigen(Sig)$vector[,1]
  
  # Sign correction
  resid_plus_pca <- v_pca + c(1, rep(0,p-1));
  resid_minus_pca <- v_pca - c(1, rep(0,p-1));
  
  ifelse(norm_2(resid_plus_pca) < norm_2(resid_minus_pca), resid_pca <- resid_plus_pca, resid_pca <- resid_minus_pca)
  
  return(cbind(spca = sqrt(n) * resid_spat, pca = sqrt(n) * resid_pca ))  # Return result for one replication
}





# 
# 
# ############ Run a small simulation in paralel
# 
# #num_cores <- detectCores()-1
# 
# # Set up parallel backend
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)
# 
# # Run parallel computation for whole v
# results_list <- foreach(i = 1:num_cores, .combine = cbind, .packages = "MASS") %dopar% {
#   replicate(reps_per_core, run_simulation_unit_v(), simplify = "array")
# }
# 
# # # Run parallel computation for unit v
# # results_list <- foreach(i = 1:num_cores, .combine = cbind, .packages = "MASS") %dopar% {
# #   replicate(reps_per_core, run_simulation_unit_v(), simplify = "array")
# # }
# 
# # Stop the parallel cluster
# stopCluster(cl)
# 
# # Convert results to matrix format
# results_matrix <- t(results_list)
# 
# # Print first few rows of results
# head(results_matrix)
# 
# round(cov(results_matrix),4)
# limit_covariance
# 
# limit_covariance/cov(results_matrix)
# 
# 
# # limit_covariance_unit_v
# # limit_covariance_unit_v/cov(results_matrix)
# 
# 
# 
# qqnorm(c(results_matrix[,1]))
# plot(data.frame(results_matrix))
# 
# 
# 




#######################################################
#---------------FIGURE 4------------------------------#
#######################################################
ps <- c(3, 5, 10)
lambdas <- seq(2,40)
ns <- c(500, 1000, 10000)
dfs <- c(3, 5, 10)
#dfs <- c(3)
# Set up parallel backend
num_cores <- 100
cl <- makeCluster(num_cores)
registerDoParallel(cl)
set.seed(1243)
p <- 3
res_3 <- numeric(0)

for (df in dfs) {
  for (lambda in lambdas) {
    for (n in ns) {
      
      # Use foreach to parallelize over m
      results_list <- foreach(m = 1:num_cores, .combine = abind, .multicombine = TRUE, .packages = c("MASS","mvtnorm")) %dopar% {
        replicate(10, run_simulation_unit_v(n, p, lambda, df, m * p * lambda * n), simplify = "array")
      }
      
      # Separate SPCA and PCA results
      res_spca_temp <- cbind(t(results_list[, 1, ]), n = n, lambda = lambda, p = 3, df = df, method = "SPCA")
      res_pca_temp  <- cbind(t(results_list[, 2, ]), n = n, lambda = lambda, p = 3, df = df, method = "PCA")
      
      # Combine results
      res_3 <- rbind(res_3, res_pca_temp, res_spca_temp)
    }
  }
}

# Stop the cluster
stopCluster(cl)

# Save to CSV
write.csv2(res_3, "res_p3_noscaling.csv")


##############################################################
##############################################################
# Set up parallel backend
num_cores <- 100
cl <- makeCluster(num_cores)
registerDoParallel(cl)

p <- 5
res_3 <- numeric(0)

for (df in dfs) {
  for (lambda in lambdas) {
    for (n in ns) {
      
      # Use foreach to parallelize over m
      results_list <- foreach(m = 1:num_cores, .combine = abind, .multicombine = TRUE, .packages = c("MASS","mvtnorm")) %dopar% {
        replicate(10, run_simulation_unit_v(n, p, lambda, df, m * p * lambda * n), simplify = "array")
      }
      
      # Separate SPCA and PCA results
      res_spca_temp <- cbind(t(results_list[, 1, ]), n = n, lambda = lambda, p = 5, df = df, method = "SPCA")
      res_pca_temp  <- cbind(t(results_list[, 2, ]), n = n, lambda = lambda, p = 5, df = df, method = "PCA")
      
      # Combine results
      res_3 <- rbind(res_3, res_pca_temp, res_spca_temp)
    }
  }
}

# Stop the cluster
stopCluster(cl)

# Save to CSV
write.csv2(res_3, "res_p5_noscaling.csv")


##############################################################
##############################################################
# Set up parallel backend
num_cores <- 100
cl <- makeCluster(num_cores)
registerDoParallel(cl)

p <- 10
res_3 <- numeric(0)

for (df in dfs) {
  for (lambda in lambdas) {
    for (n in ns) {
      
      # Use foreach to parallelize over m
      results_list <- foreach(m = 1:num_cores, .combine = abind, .multicombine = TRUE, .packages = c("MASS","mvtnorm")) %dopar% {
        replicate(10, run_simulation_unit_v(n, p, lambda, df, m * p * lambda * n), simplify = "array")
      }
      
      # Separate SPCA and PCA results
      res_spca_temp <- cbind(t(results_list[, 1, ]), n = n, lambda = lambda, p = 10, df = df, method = "SPCA")
      res_pca_temp  <- cbind(t(results_list[, 2, ]), n = n, lambda = lambda, p = 10, df = df, method = "PCA")
      
      # Combine results
      res_3 <- rbind(res_3, res_pca_temp, res_spca_temp)
    }
  }
}

# Stop the cluster
stopCluster(cl)

# Save to CSV
write.csv2(res_3, "res_p10_noscaling.csv")


##############################################################
##############################################################
# Set up parallel backend
num_cores <- 100
set.seed(1243)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

p <- 50
res_3 <- numeric(0)

for (df in dfs) {
  for (lambda in lambdas) {
    for (n in ns) {
      
      # Use foreach to parallelize over m
      results_list <- foreach(m = 1:num_cores, .combine = abind, .multicombine = TRUE, .packages = c("MASS","mvtnorm")) %dopar% {
        replicate(10, run_simulation_unit_v(n, p, lambda, df, m * p * lambda * n), simplify = "array")
      }
      
      # Separate SPCA and PCA results
      res_spca_temp <- cbind(t(results_list[, 1, ]), n = n, lambda = lambda, p = 50, df = df, method = "SPCA")
      res_pca_temp  <- cbind(t(results_list[, 2, ]), n = n, lambda = lambda, p = 50, df = df, method = "PCA")
      
      # Combine results
      res_3 <- rbind(res_3, res_pca_temp, res_spca_temp)
    }
  }
}

# Stop the cluster
stopCluster(cl)

# Save to CSV
write.csv2(res_3, "res_p50_noscaling.csv")



##############################################################
##############################################################
# Set up parallel backend
num_cores <- 100
set.seed(1243)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

p <- 100
res_3 <- numeric(0)

for (df in dfs) {
  for (lambda in lambdas) {
    for (n in ns) {
      
      # Use foreach to parallelize over m
      results_list <- foreach(m = 1:num_cores, .combine = abind, .multicombine = TRUE, .packages = c("MASS","mvtnorm")) %dopar% {
        replicate(10, run_simulation_unit_v(n, p, lambda, df, m * p * lambda * n), simplify = "array")
      }
      
      # Separate SPCA and PCA results
      res_spca_temp <- cbind(t(results_list[, 1, ]), n = n, lambda = lambda, p = 100, df = df, method = "SPCA")
      res_pca_temp  <- cbind(t(results_list[, 2, ]), n = n, lambda = lambda, p = 100, df = df, method = "PCA")
      
      # Combine results
      res_3 <- rbind(res_3, res_pca_temp, res_spca_temp)
    }
  }
}

# Stop the cluster
stopCluster(cl)

# Save to CSV
write.csv2(res_3, "res_p100_noscaling.csv")



#Plots#


mse <- function(v){
  v <- as.numeric(unlist(v))
  p <- length(v)
  # Sign correction
  resid_plus <- v + c(1, rep(0,p-1));
  resid_minus <- v - c(1, rep(0,p-1));
  
  ifelse(norm_2(resid_plus) < norm_2(resid_minus), resid_spat <- norm_2(resid_plus), resid_spat <- norm_2(resid_minus))
  return(resid_spat)
}


p3_data <- read.csv2("res_p3_noscaling.csv")
p3_data <- as.data.frame(p3_data)

p5_data <- read.csv2("res_p5_noscaling.csv")
p5_data <- as.data.frame(p5_data)

p10_data <- read.csv2("res_p10_noscaling.csv")
p10_data <- as.data.frame(p10_data)
#p10_data$p <- 10

p50_data <- read.csv2("res_p50_noscaling.csv")
p50_data <- as.data.frame(p50_data)

p100_data <- read.csv2("res_p100_noscaling.csv")
p100_data <- as.data.frame(p100_data)

library(dplyr)
library(purrr)


df <- as.data.frame(p3_data)

# Ensure grouping variables are in proper types (e.g., numeric or factor)
df <- df %>%
  mutate(
    n = as.numeric(as.character(n)),
    lambda = as.numeric(as.character(lambda)),
    p = as.numeric(as.character(p)),
    df = as.numeric(as.character(df))
  )

# Group and calculate
result <- df %>%
  group_by(n, lambda, p, df,method) %>%
  group_modify(~ {
    # For each group .x, extract the first three columns as numeric vectors
    mat <- t(sapply(1:nrow(.x), function(i) as.numeric(unlist(.x[i, 2:4]))))
    cov_mat <- cov(mat)
    cov_2norm <- max(svd(cov_mat)$d)  # 2-norm
    tibble(cov_2norm = cov_2norm)
  }) %>%
  ungroup()


df <- as.data.frame(p5_data)

# Ensure grouping variables are in proper types (e.g., numeric or factor)
df <- df %>%
  mutate(
    n = as.numeric(as.character(n)),
    lambda = as.numeric(as.character(lambda)),
    p = as.numeric(as.character(p)),
    df = as.numeric(as.character(df))
  )

# Group and calculate
result_5 <- df %>%
  group_by(n, lambda, p, df,method) %>%
  group_modify(~ {
    # For each group .x, extract the first three columns as numeric vectors
    mat <- t(sapply(1:nrow(.x), function(i) as.numeric(unlist(.x[i, 2:6]))))
    cov_mat <- cov(mat)
    cov_2norm <- max(svd(cov_mat)$d)  # 2-norm
    tibble(cov_2norm = cov_2norm)
  }) %>%
  ungroup()

df <- as.data.frame(p10_data)

# Ensure grouping variables are in proper types (e.g., numeric or factor)
df <- df %>%
  mutate(
    n = as.numeric(as.character(n)),
    lambda = as.numeric(as.character(lambda)),
    p = as.numeric(as.character(p)),
    df = as.numeric(as.character(df))
  )

# Group and calculate
result_10 <- df %>%
  group_by(n, lambda, p, df,method) %>%
  group_modify(~ {
    # For each group .x, extract the first three columns as numeric vectors
    mat <- t(sapply(1:nrow(.x), function(i) as.numeric(unlist(.x[i, 2:11]))))
    cov_mat <- cov(mat)
    cov_2norm <- max(svd(cov_mat)$d)  # 2-norm
    tibble(cov_2norm = cov_2norm)
  }) %>%
  ungroup()


df <- as.data.frame(p50_data)

# Ensure grouping variables are in proper types (e.g., numeric or factor)
df <- df %>%
  mutate(
    n = as.numeric(as.character(n)),
    lambda = as.numeric(as.character(lambda)),
    p = as.numeric(as.character(p)),
    df = as.numeric(as.character(df))
  )

# Group and calculate
result_50 <- df %>%
  group_by(n, lambda, p, df,method) %>%
  group_modify(~ {
    # For each group .x, extract the first three columns as numeric vectors
    mat <- t(sapply(1:nrow(.x), function(i) as.numeric(unlist(.x[i, 2:4]))))
    cov_mat <- cov(mat)
    cov_2norm <- max(svd(cov_mat)$d)  # 2-norm
    tibble(cov_2norm = cov_2norm)
  }) %>%
  ungroup()



df <- as.data.frame(p100_data)

# Ensure grouping variables are in proper types (e.g., numeric or factor)
df <- df %>%
  mutate(
    n = as.numeric(as.character(n)),
    lambda = as.numeric(as.character(lambda)),
    p = as.numeric(as.character(p)),
    df = as.numeric(as.character(df))
  )

# Group and calculate
result_100 <- df %>%
  group_by(n, lambda, p, df,method) %>%
  group_modify(~ {
    # For each group .x, extract the first three columns as numeric vectors
    mat <- t(sapply(1:nrow(.x), function(i) as.numeric(unlist(.x[i, 2:4]))))
    cov_mat <- cov(mat)
    cov_2norm <- max(svd(cov_mat)$d)  # 2-norm
    tibble(cov_2norm = cov_2norm)
  }) %>%
  ungroup()




result_all <- rbind(result, result_5, result_10, result_50,result_100)
write.csv2(result_all, "res_cov.csv")

result_all <- (read.csv2("res_cov.csv"))[,-1]
result_all$n <- as.factor(result_all$n)
levels(result_all$n) <- c("n = 500", "n = 1000", "n = 10000")

result_all$p <- as.factor(result_all$p)
levels(result_all$p) <- c("p = 3", "p = 5", "p = 10", "p=50","p=100")

library(ggplot2)
ggplot(result_all, aes(x = lambda, y = log(cov_2norm), color = method, linetype = as.factor(df))) +
  geom_line() +
  facet_grid(p~n) +
  theme_minimal() +
  labs(
    #title = "MSE vs Lambda",
    x = expression(lambda),
    y = expression("log || ASCOV ||"[2]),#expression("log||" * Cov(w[n]) * "||"[2]),
    color = "",
    linetype = expression(nu)
  ) + 
  theme_bw()+
  theme(legend.position = "bottom")



############### Auxilliary plot - mean difference ################
###############

df <- as.data.frame(p3_data)

# Ensure grouping variables are in proper types (e.g., numeric or factor)
df <- df %>%
  mutate(
    n = as.numeric(as.character(n)),
    lambda = as.numeric(as.character(lambda)),
    p = as.numeric(as.character(p)),
    df = as.numeric(as.character(df))
  )

# Group and calculate
result <- df %>%
  group_by(n, lambda, p, df,method) %>%
  group_modify(~ {
    # For each group .x, extract the first three columns as numeric vectors
    mat <- t(sapply(1:nrow(.x), function(i) as.numeric(unlist(.x[i, 2:4]))))
    mean1 <- sum(mat^2)
    tibble(mean1 = mean1/n)
  }) %>%
  ungroup()


df <- as.data.frame(p5_data)

# Ensure grouping variables are in proper types (e.g., numeric or factor)
df <- df %>%
  mutate(
    n = as.numeric(as.character(n)),
    lambda = as.numeric(as.character(lambda)),
    p = as.numeric(as.character(p)),
    df = as.numeric(as.character(df))
  )

# Group and calculate
result_5 <- df %>%
  group_by(n, lambda, p, df,method) %>%
  group_modify(~ {
    # For each group .x, extract the first three columns as numeric vectors
    mat <- t(sapply(1:nrow(.x), function(i) as.numeric(unlist(.x[i, 2:6]))))
    mean1 <- sum(mat^2)
    tibble(mean1 = mean1/n)
  }) %>%
  ungroup()

df <- as.data.frame(p10_data)

# Ensure grouping variables are in proper types (e.g., numeric or factor)
df <- df %>%
  mutate(
    n = as.numeric(as.character(n)),
    lambda = as.numeric(as.character(lambda)),
    p = as.numeric(as.character(p)),
    df = as.numeric(as.character(df))
  )

# Group and calculate
result_10 <- df %>%
  group_by(n, lambda, p, df,method) %>%
  group_modify(~ {
    # For each group .x, extract the first three columns as numeric vectors
    mat <- t(sapply(1:nrow(.x), function(i) as.numeric(unlist(.x[i, 2:11]))))
    mean1 <- sum(mat^2)
    tibble(mean1 = mean1/n)
  }) %>%
  ungroup()


df <- as.data.frame(p50_data)

# Ensure grouping variables are in proper types (e.g., numeric or factor)
df <- df %>%
  mutate(
    n = as.numeric(as.character(n)),
    lambda = as.numeric(as.character(lambda)),
    p = as.numeric(as.character(p)),
    df = as.numeric(as.character(df))
  )

# Group and calculate
result_50 <- df %>%
  group_by(n, lambda, p, df,method) %>%
  group_modify(~ {
    # For each group .x, extract the first three columns as numeric vectors
    mat <- t(sapply(1:nrow(.x), function(i) as.numeric(unlist(.x[i, 2:4]))))
    mean1 <- sum(mat^2)
    tibble(mean1 = mean1/n)
  }) %>%
  ungroup()



df <- as.data.frame(p100_data)

# Ensure grouping variables are in proper types (e.g., numeric or factor)
df <- df %>%
  mutate(
    n = as.numeric(as.character(n)),
    lambda = as.numeric(as.character(lambda)),
    p = as.numeric(as.character(p)),
    df = as.numeric(as.character(df))
  )

# Group and calculate
result_100 <- df %>%
  group_by(n, lambda, p, df,method) %>%
  group_modify(~ {
    # For each group .x, extract the first three columns as numeric vectors
    mat <- t(sapply(1:nrow(.x), function(i) as.numeric(unlist(.x[i, 2:4]))))
    mean1 <- sum(mat^2)
    tibble(mean1 = mean1/n)
  }) %>%
  ungroup()




result_mean <- rbind(result, result_5, result_10, result_50,result_100)
write.csv2(result_mean, "res_mean.csv")

result_mean <- (read.csv2("res_mean.csv"))[,-1]

result_mean$n <- as.factor(result_mean$n)
levels(result_mean$n) <- c("n = 500", "n = 1000", "n = 10000")

result_all$p <- as.factor(result_mean$p)
levels(result_all$p) <- c("p = 3", "p = 5", "p = 10", "p=50","p=100")

library(ggplot2)
ggplot(result_mean, aes(x = lambda, y = log(mean1), color = method, linetype = as.factor(df))) +
  geom_line() +
  facet_grid(n~p,scales="fixed") +
  theme_minimal() +
  labs(
    #title = "MSE vs Lambda",
    x = expression(lambda),
    y = expression( "log( Mean ||"~w[n]^2~"||)" ),
    color = "",
    linetype = expression(nu)
  ) + 
  theme_bw()+
  theme(legend.position = "bottom")

















###-------------New plot--------------------------###

spatial_simu_asym <- function(n=1000,p=5,lambda=20,df=3,seed=set.seed()){
  
  x <- mvtnorm::rmvt(n = n, sigma = diag(c(lambda, rep(1,p-1) ) ), df = df );
  res <- estim_v(x, v=c(1,0,0), maxiter = 100)
  
  #mean(sqrt(apply((x-matrix(rep(res$v,n),byrow=T,ncol=3) )^2,1,sum))/sqrt(apply((x+matrix(rep(res$v,n),byrow=T,ncol=3))^2,1,sum)))
  
  
  #x <- sweep(x, 1, c(2, 1), "*")
  
  #res <- estim_v(x, maxiter = 100)
  #res
  
  v <- res$v
  
  #Sigma_eigens <- eigen(cov(x))$values
  #sigma1 <- Sigma_eigens[1]
  #sigma2 <- Sigma_eigens[2]
  
  sigma1 <- (sigma*4)^2;
  sigma2 <- sigma^2;
  
  psi <- norm(res$v,"2")
  
  o1 <- c(1,rep(0,p-1))
  o2 <- c(0,1,rep(0,p-2))
  
  v_matrix <- matrix(rep(psi*o1,n),nrow=n,byrow=TRUE )
  
  x_minus <- -1 * (x-v_matrix)
  x_plus <- x+v_matrix
  
  y_minus <- apply(x_minus,1,norm_2)
  y_plus <- apply(x_plus,1,norm_2)
  
  
  ##Hessian part
  
  a1 <- mean( y_minus/y_plus+y_plus/y_minus )
  
  a2 <- mean( y_plus / y_minus^3 * (tcrossprod(o1, x_minus) )^2 )
  
  a3 <- mean( y_minus / y_plus^3 * (tcrossprod(o1, x_plus) )^2 )
  
  a4 <- mean( 1 / (y_minus * y_plus ) * tcrossprod(o1, x_plus) * tcrossprod(o1, x_minus) )
  
  
  b1 <- a1
  
  b2 <- mean( y_plus / y_minus^3 * (tcrossprod(o2, x_minus) )^2 )
  
  b3 <- mean( y_minus / y_plus^3 * (tcrossprod(o2, x_plus) )^2 )
  
  b4 <- mean( 1 / (y_minus * y_plus ) * tcrossprod(o2, x_plus) * tcrossprod(o2, x_minus) )
  
  a_hessian <- (a1 - ( a2 + a3 ) + 2*a4)
  b_hessian <- (b1 - ( b2 + b3 ) + 2*b4)
  
  
  ## Gradient part
  
  
  c11 <- mean( y_plus^2 / y_minus^2 * (tcrossprod(o1, x_minus) )^2 )
  c12 <- mean( y_minus^2 / y_plus^2 * (tcrossprod(o1, x_plus) )^2 )
  
  c1 <- ( c11 + c12 )
  
  d11 <- mean( y_plus^2 / y_minus^2 * (tcrossprod(o2, x_minus) )^2 )
  d12 <- mean( y_minus^2 / y_plus^2 * (tcrossprod(o2, x_plus) )^2 )
  
  d1 <- ( d11 + d12 )
  
  c_grad <- c1 + 2*psi^2 - 2*sigma1
  d_grad <- d1 - 2*sigma2
  
  c <- c_grad / a_hessian^2
  d <- d_grad / b_hessian^2
  
  limit_covariance <- c * tcrossprod(o1,o1) + d * ( diag(p) - tcrossprod(o1,o1))
  
  limit_covariance_unit_v <-  d / psi^2 * ( diag(p) - tcrossprod(o1,o1))
  
  return(d/psi^2)
}

pca_simu_asym <- function(n=1000,p=5,lambda=20,df=3,seed=set.seed()){
  z <- mvtnorm::rmvt(n = n, sigma = diag(p), df = df );
  x <- z %*% diag(c(sqrt(lambda), rep(1,p-1) ) )
  z1 <- z[,1]
  z2 <- z[,2]
  
  const1 <- mean(z1^2*z2^2)/(mean(z1^2))^2*lambda/(lambda-1)^2
  
  return(const1)
}

ps <- c(3,5,10)
lambdas <- seq(2,40)
ns <- c(500, 1000, 50000)
dfs <- c(3, 5, 10)

pca_res <- matrix(ncol=6,nrow=3*3*3*length(lambdas));
pca_res <- as.data.frame(pca_res)
names(pca_res) <- c("SPCA","PCA","n","p","df","lambda")

i <- 0;

for (df in dfs) {
  for (lambda in lambdas) {
    for (n in ns) {
      for (p in ps) {
        i <- i+1;
        pca_temp <- pca_simu_asym(n=n,df=df,p=p,lambda=lambda)
        spca_temp <- spatial_simu_asym(n=n,df=df,p=p,lambda=lambda)
        
        pca_res[i,] <- c(spca_temp,pca_temp,n,p,df,lambda)
      }
    }
  }
}

pca_res_final <- cbind(const=c(pca_res$SPCA,pca_res$PCA),method=rep(c("SPCA","PCA"),each=dim(pca_res)[1]), rbind(pca_res[,-c(1:2)],pca_res[,-c(1:2)]))
write.csv2(pca_res_final,"Constant_simu_res.csv")

pca_res_final <- read.csv2("Constant_simu_res.csv")



ggplot(pca_res_final[pca_res_final$lambda>5,], aes(x = lambda, y = log(const), color = method, linetype = as.factor(df))) +
  geom_line() +
  facet_grid(p~n,   scales = "free") +
  theme_minimal() +
  labs(
    #title = "MSE vs Lambda",
    x = expression(lambda),
    y = expression("log||" * Cov(w[n]) * "||"[2]),
    color = "",
    linetype = expression(nu)
  ) + 
  theme_bw()+
  theme(legend.position = "bottom")







