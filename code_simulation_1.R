#-------------
# Simulation study #1
#-------------


library(tidyverse)


#
# Run the methods from code_algorithm.R
#



#-----------------------------------------------
# PART 1: Obtain the empirical estimates




# Function that generates a single sample from the model (theta^2, 1, 1) used in simulation study #1
# n = sample size
# my_sd = theta
# distribution = (1 = normal, 2 = uniform, 3 = Bernoulli)
single_simu <- function(n, my_sd, distribution){
  if(distribution == 1){
    x <- cbind(rnorm(n, 0, my_sd), rnorm(n), rnorm(n))
  }
  if(distribution == 2){
    x <- cbind(runif(n, -sqrt(3)*my_sd, sqrt(3)*my_sd), runif(n, -sqrt(3), sqrt(3)), runif(n, -sqrt(3), sqrt(3)))
  }
  if(distribution == 3){
    x <- cbind(sample(c(-my_sd, my_sd), n, replace = TRUE), sample(c(-1, 1), n, replace = TRUE), sample(c(-1, 1), n, replace = TRUE))
  }
  
  res <- tryCatch(estim_v(x, maxiter = 500)$v, error = function(e) NA)
  
  if(any(is.na(res))){
    my_v <- c(NA, NA, NA)
  }
  else{
    my_v <- res
  }
  
  data.frame(n = n, my_sd = my_sd, distribution = distribution, v1 = my_v[1], v2 = my_v[2], v3 = my_v[3])
}

# single_simu(100, 3, 2)



# Parameter grid
sd_set <- seq(1, 3, length.out = 20)
n_set <- c(20, 80, 320, 1280)
dist_set <- c(1, 2, 3)

par_grid <- expand.grid(n = n_set,
                        my_sd = sd_set,
                        distribution = dist_set)



# Run the below parallelized code twice, setting batch first to 1 then to 2
# Batch 1 = 500 reps
# Batch 2 = 500 reps

batch <- 1

set.seed(20250000 + batch)

res <- data.frame()



# Simu

library(parallel)

par_n <- nrow(par_grid)

reps <- 500


for(i in 1:reps){
  cl <- makeCluster(4)
  
  clusterExport(cl, c("single_simu",
                      "eval_grad",
                      "eval_f",
                      "eval_f_lambda",
                      "T_iter",
                      "Tk_iter",
                      "estim_v",
                      "par_grid"), envir = environment())
  
  # clusterExport(cl, ls(), envir = environment())
  
  temp <- parLapply(cl, 1:par_n, function(j) single_simu(par_grid$n[j], par_grid$my_sd[j], par_grid$distribution[j]))
  stopCluster(cl)
  
  temp <- do.call("rbind", temp)
  
  res <- rbind(res, temp)
  print(i)
}


write.table(res, paste0("FOLDER_HERE/simu_1_res/res_", batch, ".txt"))





# Function that generates a single sample from the model (theta^2, theta, 1) used in simulation study #1
# n = sample size
# my_sd = theta
# distribution = (1 = normal, 2 = uniform, 3 = Bernoulli)
single_simu <- function(n, my_sd, distribution){
  if(distribution == 1){
    x <- cbind(rnorm(n, 0, my_sd), rnorm(n, 0, sqrt(my_sd)), rnorm(n))
  }
  if(distribution == 2){
    x <- cbind(runif(n, -sqrt(3)*my_sd, sqrt(3)*my_sd), runif(n, -sqrt(3)*sqrt(my_sd), sqrt(3)*sqrt(my_sd)), runif(n, -sqrt(3), sqrt(3)))
  }
  if(distribution == 3){
    x <- cbind(sample(c(-my_sd, my_sd), n, replace = TRUE), sample(c(-sqrt(my_sd), sqrt(my_sd)), n, replace = TRUE), sample(c(-1, 1), n, replace = TRUE))
  }
  
  res <- tryCatch(estim_v(x, maxiter = 500)$v, error = function(e) NA)
  
  if(any(is.na(res))){
    my_v <- c(NA, NA, NA)
  }
  else{
    my_v <- res
  }
  
  data.frame(n = n, my_sd = my_sd, distribution = distribution, v1 = my_v[1], v2 = my_v[2], v3 = my_v[3])
}

# single_simu(10, 2, 3)



# Parameter grid
sd_set <- seq(1, 3, length.out = 20)
n_set <- c(20, 80, 320, 1280)
dist_set <- c(1, 2, 3)

par_grid <- expand.grid(n = n_set,
                        my_sd = sd_set,
                        distribution = dist_set)


# Run the below parallelized code twice, setting batch first to 1 then to 2
# Batch 1 = 500 reps
# Batch 2 = 500 reps

batch <- 1

set.seed(20251000 + batch)

res <- data.frame()



# Simu

library(parallel)

par_n <- nrow(par_grid)

reps <- 500


for(i in 1:reps){
  cl <- makeCluster(4)
  
  clusterExport(cl, c("single_simu",
                      "eval_grad",
                      "eval_f",
                      "eval_f_lambda",
                      "T_iter",
                      "Tk_iter",
                      "estim_v",
                      "par_grid"), envir = environment())
  
  # clusterExport(cl, ls(), envir = environment())
  
  temp <- parLapply(cl, 1:par_n, function(j) single_simu(par_grid$n[j], par_grid$my_sd[j], par_grid$distribution[j]))
  stopCluster(cl)
  
  temp <- do.call("rbind", temp)
  
  res <- rbind(res, temp)
  print(i)
}


write.table(res, paste0("FOLDER_HERE/simu_2_res/res_", batch, ".txt"))






#-----------------------------------------------
# PART 2: Approximate the true values of psi





# Function that generates a single sample from the model (theta^2, 1, 1) used in simulation study #1
# n = sample size
# my_sd = theta
# distribution = (1 = normal, 2 = uniform, 3 = Bernoulli)
single_simu <- function(n, my_sd, distribution){
  if(distribution == 1){
    x <- cbind(rnorm(n, 0, my_sd), rnorm(n), rnorm(n))
  }
  if(distribution == 2){
    x <- cbind(runif(n, -sqrt(3)*my_sd, sqrt(3)*my_sd), runif(n, -sqrt(3), sqrt(3)), runif(n, -sqrt(3), sqrt(3)))
  }
  if(distribution == 3){
    x <- cbind(sample(c(-my_sd, my_sd), n, replace = TRUE), sample(c(-1, 1), n, replace = TRUE), sample(c(-1, 1), n, replace = TRUE))
  }
  
  res <- tryCatch(estim_v(x, maxiter = 500)$v, error = function(e) NA)
  
  if(any(is.na(res))){
    my_v <- c(NA, NA, NA)
  }
  else{
    my_v <- res
  }
  
  data.frame(n = n, my_sd = my_sd, distribution = distribution, v1 = my_v[1], v2 = my_v[2], v3 = my_v[3])
}

# single_simu(100, 3, 2)


# my_res <- single_simu(20000, 3, 2)
# sqrt(my_res$v1^2 + my_res$v3^2 + my_res$v3^2)



# Auxiliary function for estimating psi
estim_psi <- function(my_sd, distribution, n = 20000){
  my_res <- single_simu(n, my_sd, distribution)
  psi <- sqrt(my_res$v1^2 + my_res$v3^2 + my_res$v3^2)
  
  data.frame(my_sd = my_sd, distribution = distribution, psi = psi)
}






# Parameter grid
sd_set <- seq(1, 3, length.out = 20)
dist_set <- c(1, 2, 3)

par_grid <- expand.grid(my_sd = sd_set,
                        distribution = dist_set)



# Only 1 batch here
batch <- 1

set.seed(20260000 + batch)

res <- data.frame()



# Simu

library(parallel)

par_n <- nrow(par_grid)

reps <- 100


for(i in 1:reps){
  cl <- makeCluster(4)
  
  clusterExport(cl, c("single_simu",
                      "eval_grad",
                      "eval_f",
                      "eval_f_lambda",
                      "T_iter",
                      "Tk_iter",
                      "estim_v",
                      "par_grid",
                      "estim_psi"), envir = environment())
  
  # clusterExport(cl, ls(), envir = environment())
  
  temp <- parLapply(cl, 1:par_n, function(j) estim_psi(par_grid$my_sd[j], par_grid$distribution[j], 10000))
  stopCluster(cl)
  
  temp <- do.call("rbind", temp)
  
  res <- rbind(res, temp)
  print(i)
}


res_agg <- res %>%
  group_by(my_sd, distribution) %>%
  summarise(psi_avg = mean(psi, na.rm = TRUE))


write.table(res_agg, paste0("FOLDER_HERE/simu_1_res/psi_estim.txt"))







# Function that generates a single sample from the model (theta^2, theta, 1) used in simulation study #1
# n = sample size
# my_sd = theta
# distribution = (1 = normal, 2 = uniform, 3 = Bernoulli)
single_simu <- function(n, my_sd, distribution){
  if(distribution == 1){
    x <- cbind(rnorm(n, 0, my_sd), rnorm(n, 0, sqrt(my_sd)), rnorm(n))
  }
  if(distribution == 2){
    x <- cbind(runif(n, -sqrt(3)*my_sd, sqrt(3)*my_sd), runif(n, -sqrt(3)*sqrt(my_sd), sqrt(3)*sqrt(my_sd)), runif(n, -sqrt(3), sqrt(3)))
  }
  if(distribution == 3){
    x <- cbind(sample(c(-my_sd, my_sd), n, replace = TRUE), sample(c(-sqrt(my_sd), sqrt(my_sd)), n, replace = TRUE), sample(c(-1, 1), n, replace = TRUE))
  }
  
  res <- tryCatch(estim_v(x, maxiter = 500)$v, error = function(e) NA)
  
  if(any(is.na(res))){
    my_v <- c(NA, NA, NA)
  }
  else{
    my_v <- res
  }
  
  data.frame(n = n, my_sd = my_sd, distribution = distribution, v1 = my_v[1], v2 = my_v[2], v3 = my_v[3])
}

# single_simu(100, 3, 2)


# my_res <- single_simu(20000, 3, 2)
# sqrt(my_res$v1^2 + my_res$v3^2 + my_res$v3^2)


# Auxiliary function for estimating psi
estim_psi <- function(my_sd, distribution, n = 20000){
  my_res <- single_simu(n, my_sd, distribution)
  psi <- sqrt(my_res$v1^2 + my_res$v3^2 + my_res$v3^2)
  
  data.frame(my_sd = my_sd, distribution = distribution, psi = psi)
}


# Parameter grid
sd_set <- seq(1, 3, length.out = 20)
dist_set <- c(1, 2, 3)


par_grid <- expand.grid(my_sd = sd_set,
                        distribution = dist_set)


# Only 1 batch here
batch <- 1

set.seed(20260000 + batch)

res <- data.frame()



# Simu

library(parallel)

par_n <- nrow(par_grid)

reps <- 100


for(i in 1:reps){
  cl <- makeCluster(4)
  
  clusterExport(cl, c("single_simu",
                      "eval_grad",
                      "eval_f",
                      "eval_f_lambda",
                      "T_iter",
                      "Tk_iter",
                      "estim_v",
                      "par_grid",
                      "estim_psi"), envir = environment())
  
  # clusterExport(cl, ls(), envir = environment())
  
  temp <- parLapply(cl, 1:par_n, function(j) estim_psi(par_grid$my_sd[j], par_grid$distribution[j], 10000))
  stopCluster(cl)
  
  temp <- do.call("rbind", temp)
  
  res <- rbind(res, temp)
  print(i)
}


res_agg <- res %>%
  group_by(my_sd, distribution) %>%
  summarise(psi_avg = mean(psi, na.rm = TRUE))


write.table(res_agg, paste0("FOLDER_HERE/simu_2_res/psi_estim.txt"))





#-----------------------------------------------
# PART 3: Create Figure 2 (directions)



# Load results of first scenario
res <- NULL

for(i in 1:2){
  temp <- read.table(paste0("FOLDER_HERE/simu_1_res/res_", i, ".txt"))
  res <- bind_rows(res, temp)
}


res_agg <- res %>%
  mutate(distribution = factor(distribution, levels = c(1, 2, 3), labels = c("Normal", "Uniform", "Bernoulli"))) %>%
  group_by(n, my_sd, distribution) %>%
  summarise(me_q1 = median(abs(v1), na.rm = TRUE),
            me_q2 = median(sqrt(v2^2 + v3^2), na.rm = TRUE),
            lo_q1 = quantile(abs(v1), probs = 0.25, na.rm = TRUE),
            hi_q1 = quantile(abs(v1), probs = 0.75, na.rm = TRUE),
            lo_q2 = quantile(sqrt(v2^2 + v3^2), probs = 0.25, na.rm = TRUE),
            hi_q2 = quantile(sqrt(v2^2 + v3^2), probs = 0.75, na.rm = TRUE))

res_temp_1 <- res_agg %>%
  select(-me_q2, -lo_q2, -hi_q2) %>%
  rename(me = me_q1, lo = lo_q1, hi = hi_q1) %>%
  mutate(quant = 1)

res_temp_2 <- res_agg %>%
  select(-me_q1, -lo_q1, -hi_q1) %>%
  rename(me = me_q2, lo = lo_q2, hi = hi_q2) %>%
  mutate(quant = 2)

res_agg_final <- bind_rows(res_temp_1, res_temp_2) %>%
  mutate(quant = factor(quant, levels = c(1, 2), labels = c("alab", "blab")))

res_1 <- res_agg_final %>%
  mutate(model = 1)




# Load results of second scenario
res <- NULL

for(i in 1:2){
  temp <- read.table(paste0("FOLDER_HERE/simu_2_res/res_", i, ".txt"))
  res <- bind_rows(res, temp)
}


res_agg <- res %>%
  mutate(distribution = factor(distribution, levels = c(1, 2, 3), labels = c("Normal", "Uniform", "Bernoulli"))) %>%
  group_by(n, my_sd, distribution) %>%
  summarise(me_q1 = median(abs(v1), na.rm = TRUE),
            me_q2 = median(sqrt(v2^2 + v3^2), na.rm = TRUE),
            lo_q1 = quantile(abs(v1), probs = 0.25, na.rm = TRUE),
            hi_q1 = quantile(abs(v1), probs = 0.75, na.rm = TRUE),
            lo_q2 = quantile(sqrt(v2^2 + v3^2), probs = 0.25, na.rm = TRUE),
            hi_q2 = quantile(sqrt(v2^2 + v3^2), probs = 0.75, na.rm = TRUE))

res_temp_1 <- res_agg %>%
  select(-me_q2, -lo_q2, -hi_q2) %>%
  rename(me = me_q1, lo = lo_q1, hi = hi_q1) %>%
  mutate(quant = 1)

res_temp_2 <- res_agg %>%
  select(-me_q1, -lo_q1, -hi_q1) %>%
  rename(me = me_q2, lo = lo_q2, hi = hi_q2) %>%
  mutate(quant = 2)

res_agg_final <- bind_rows(res_temp_1, res_temp_2) %>%
  mutate(quant = factor(quant, levels = c(1, 2), labels = c("alab", "blab")))

res_2 <- res_agg_final %>%
  mutate(model = 2)




# Combine the results
res_final <- bind_rows(res_1, res_2) %>%
  mutate(model = factor(model,
                        levels = c(1, 2),
                        labels = c("Cov == diag(theta^2 * ', ' * 1 * ', ' * 1)",
                                   "Cov == diag(theta^2 * ', ' * theta * ', ' * 1)")))

lt_labels <- c("alab" = expression("|v"[1]*"|"), "blab" = expression(sqrt("v"[2]^2~"+v"[3]^2)))


# Figure 2
ggplot(res_final, aes(x = my_sd, y = me, col = as.factor(n), fill = as.factor(n))) +
  geom_ribbon(aes(ymin = lo, ymax = hi, linetype = quant), alpha = 0.25, colour = NA) +
  geom_line(aes(linetype = quant)) +
  facet_grid(model ~ distribution, labeller = labeller(model = label_parsed)) +
  labs(x = expression("Standard deviation"~theta~"of PC1"), y = "Median value", col = "Sample size", linetype = "Quantity") +
  theme_bw() +
  scale_linetype_discrete(labels = lt_labels) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 1.815271, col = "darkgrey") +
  scale_fill_discrete(guide = "none")





#-----------------------------------------------
# PART 4: Create Figure 3 (magnitudes)



# Load results of first scenario
res <- NULL

for(i in 1:2){
  temp <- read.table(paste0("FOLDER_HERE/simu_1_res/res_", i, ".txt"))
  res <- bind_rows(res, temp)
}

psi_res <- read.table("FOLDER_HERE/simu_1_res/psi_estim.txt")

psi_res <- psi_res %>%
  mutate(distribution = factor(distribution, levels = c(1, 2, 3), labels = c("Normal", "Uniform", "Bernoulli")))

res_2 <- res %>%
  mutate(distribution = factor(distribution, levels = c(1, 2, 3), labels = c("Normal", "Uniform", "Bernoulli"))) %>%
  mutate(psi_estim = sqrt(v1^2 + v2^2 + v3^2)) %>%
  select(-v1, -v2, -v3)

res_agg <- res_2 %>%
  group_by(n, my_sd, distribution) %>%
  summarise(me = median(psi_estim, na.rm = TRUE),
            lo = quantile(psi_estim, probs = 0.25, na.rm = TRUE),
            hi = quantile(psi_estim, probs = 0.75, na.rm = TRUE))

res_3 <- res_agg %>%
  left_join(psi_res, by = c("my_sd", "distribution"))


res_psi_1 <- res_3 %>%
  mutate(model = 1)





# Load results of second scenario
res <- NULL

for(i in 1:2){
  temp <- read.table(paste0("FOLDER_HERE/simu_2_res/res_", i, ".txt"))
  res <- bind_rows(res, temp)
}

psi_res <- read.table("FOLDER_HERE/simu_2_res/psi_estim.txt")

psi_res <- psi_res %>%
  mutate(distribution = factor(distribution, levels = c(1, 2, 3), labels = c("Normal", "Uniform", "Bernoulli")))

res_2 <- res %>%
  mutate(distribution = factor(distribution, levels = c(1, 2, 3), labels = c("Normal", "Uniform", "Bernoulli"))) %>%
  mutate(psi_estim = sqrt(v1^2 + v2^2 + v3^2)) %>%
  select(-v1, -v2, -v3)

res_agg <- res_2 %>%
  group_by(n, my_sd, distribution) %>%
  summarise(me = median(psi_estim, na.rm = TRUE),
            lo = quantile(psi_estim, probs = 0.25, na.rm = TRUE),
            hi = quantile(psi_estim, probs = 0.75, na.rm = TRUE))


res_3 <- res_agg %>%
  left_join(psi_res, by = c("my_sd", "distribution"))


res_psi_2 <- res_3 %>%
  mutate(model = 2)



# Combine the results
res_final <- bind_rows(res_psi_1, res_psi_2) %>%
  mutate(model = factor(model,
                        levels = c(1, 2),
                        labels = c("Cov == diag(theta^2 * ', ' * 1 * ', ' * 1)",
                                   "Cov == diag(theta^2 * ', ' * theta * ', ' * 1)")))


psi_avg_data <- res_final %>%
  distinct(my_sd, distribution, model, psi_avg)


# Figure 3
ggplot(res_final, aes(x = my_sd, y = me, col = as.factor(n), fill = as.factor(n))) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, colour = NA) +
  geom_line() +
  geom_line(data = psi_avg_data, 
            aes(x = my_sd, y = psi_avg, 
                group = interaction(model, distribution),
                col = "psi_avg"),
            linetype = "dashed", linewidth = 0.7, inherit.aes = FALSE) +
  facet_grid(model ~ distribution, labeller = labeller(model = label_parsed)) +
  labs(x = expression("Standard deviation"~theta~"of PC1"), y = "Median value", col = "Sample size") +
  theme_bw() +
  scale_colour_manual(
    values = c("20"   = "#F8766D", 
               "80"   = "#7CAE00", 
               "320"  = "#00BCD8", 
               "1280" = "#C77CFF", 
               "psi_avg" = "black"),
    labels = c("20", "80", "320", "1280", "Limit")
  ) +
  scale_fill_manual(
    values = c("20"   = "#F8766D", 
               "80"   = "#7CAE00", 
               "320"  = "#00BCD8", 
               "1280" = "#C77CFF"),
    guide = "none"
  ) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 1.815271, col = "darkgrey") +
  geom_abline(slope = 1, intercept = 0, col = "lightgrey", linewidth = 0.7, linetype = "dashed")




