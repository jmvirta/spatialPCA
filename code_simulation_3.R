#-------------
# Simulation study #3
#-------------


library(tidyverse)


#
# Run the methods from code_algorithm.R
#


#-----------------------------------------------
# PART 1: Obtain the empirical estimates



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
  
  res <- tryCatch({
    v1 <- estim_v(x, maxiter = 500)$v
    Q <- diag(3) - v1%*%t(v1)/sum(v1^2)
    x_proj <- x%*%Q
    estim_v(x_proj, maxiter = 500)$v
    }, error = function(e) NA)
  
  if(any(is.na(res))){
    my_v <- c(NA, NA, NA)
  }
  else{
    my_v <- res
  }
  
  data.frame(n = n, my_sd = my_sd, distribution = distribution, v1 = my_v[1], v2 = my_v[2], v3 = my_v[3])
}

# single_simu(100, 2, 1)


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

batch <- 2

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


write.table(res, paste0("FOLDER_HERE/second_PC_res/res_", batch, ".txt"))




#-----------------------------------------------
# PART 2: Create Figure 5 


# Load results

res <- NULL

for(i in 1:2){
  temp <- read.table(paste0("FOLDER_HERE/second_PC_res/res_", i, ".txt"))
  res <- bind_rows(res, temp)
}



# Process the raw results

res_agg <- res %>%
  mutate(distribution = factor(distribution, levels = c(1, 2, 3), labels = c("Normal", "Uniform", "Bernoulli"))) %>%
  group_by(n, my_sd, distribution) %>%
  summarise(me_q1 = median(abs(v2), na.rm = TRUE),
            me_q2 = median(sqrt(v1^2 + v3^2), na.rm = TRUE),
            lo_q1 = quantile(abs(v2), probs = 0.25, na.rm = TRUE),
            hi_q1 = quantile(abs(v2), probs = 0.75, na.rm = TRUE),
            lo_q2 = quantile(sqrt(v1^2 + v3^2), probs = 0.25, na.rm = TRUE),
            hi_q2 = quantile(sqrt(v1^2 + v3^2), probs = 0.75, na.rm = TRUE))

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

res_agg_final <- res_agg_final %>%
  mutate(model = 2)



res_final <- res_agg_final %>%
  mutate(model = factor(model,
                        levels = 2,
                        labels = c("Cov == diag(theta^2 * ', ' * theta * ', ' * 1)")))

lt_labels <- c("alab" = expression("|w"[2]*"|"), "blab" = expression(sqrt("w"[1]^2~"+w"[3]^2)))


# Figure 5
ggplot(res_final, aes(x = my_sd, y = me, col = as.factor(n), fill = as.factor(n))) +
  geom_ribbon(aes(ymin = lo, ymax = hi, linetype = quant), alpha = 0.25, colour = NA) +
  geom_line(aes(linetype = quant)) +
  facet_grid(model ~ distribution, labeller = labeller(model = label_parsed)) +
  labs(x = expression("Standard deviation"~theta~"of PC1"), y = "Median value", col = "Sample size", linetype = "Quantity") +
  theme_bw() +
  scale_linetype_discrete(labels = lt_labels) +
  theme(legend.position = "bottom") +
  scale_fill_discrete(guide = "none")



