library(abc)


toy_model <- function(x) {
  return(2 * x + 5 + rnorm(1, 0, 0.1))
}

true_param <- 0.75
sum_stat_obs <- 6.5  


set.seed(123)
n_samples <- 10000
prior_samples <- runif(n_samples, 0, 10)

sum_stats <- numeric(n_samples)
for (i in 1:n_samples) {
  sum_stats[i] <- toy_model(prior_samples[i])
}


sim_data <- data.frame(param = prior_samples, stat = sum_stats)


tolerance_levels <- c(0.1, 0.3, 0.5)  
rejection_results <- list()

for (i in 1:length(tolerance_levels)) {
  rejection_results[[i]] <- abc(target = sum_stat_obs, 
                                param = sim_data$param, 
                                sumstat = sim_data$stat, 
                                tol = tolerance_levels[i], 
                                method = "rejection")
}


loclinear_results <- list()

for (i in 1:length(tolerance_levels)) {
  loclinear_results[[i]] <- abc(target = sum_stat_obs, 
                                param = sim_data$param, 
                                sumstat = sim_data$stat, 
                                tol = tolerance_levels[i], 
                                method = "loclinear")
}


rejection_params <- list()
loclinear_params <- list()

for (i in 1:length(tolerance_levels)) {
  rejection_params[[i]] <- rejection_results[[i]]$unadj.values
  loclinear_params[[i]] <- loclinear_results[[i]]$adj.values
}


par(mfrow = c(1, 1))

# Boxplot for rejection method
boxplot(rejection_params, 
        main = "ABC Rejection Method",
        xlab = "Tolerance Level", 
        ylab = "Parameter Value",
        names = tolerance_levels,  # Will now show 0.1, 0.3, 0.5
        col = "lightblue")
abline(h = true_param, col = "red", lty = 2)
legend("topright", legend = "True Parameter Value", col = "red", lty = 2)

# Boxplot for local linear regression method
boxplot(loclinear_params, 
        main = "ABC with Local Linear Regression Adjustment",
        xlab = "Tolerance Level", 
        ylab = "Parameter Value",
        names = tolerance_levels,  # Will now show 0.1, 0.3, 0.5
        col = "lightgreen")
abline(h = true_param, col = "red", lty = 2)
legend("topright", legend = "True Parameter Value", col = "red", lty = 2)

# Calculate theoretical posterior parameters
posterior_sd <- 0.05
posterior_mean <- true_param

# Create density plots for all tolerance levels
par(mfrow = c(length(tolerance_levels), 1), mar = c(4, 4, 2, 1))

for (i in 1:length(tolerance_levels)) {

  rej_dens <- density(rejection_params[[i]])
  loc_dens <- density(loclinear_params[[i]])
  

  x_seq <- seq(0, 2, length.out = 1000)
  true_dens <- dnorm(x_seq, mean = posterior_mean, sd = posterior_sd)
  

  y_max <- max(c(max(rej_dens$y), max(loc_dens$y), max(true_dens))) * 1.1
  

  plot(NULL, xlim = c(0, 2), ylim = c(0, y_max),
       main = paste("Density Comparison (Tolerance =", tolerance_levels[i], ")"),
       xlab = "Parameter Value", ylab = "Density")
  

  lines(rej_dens, col = "blue", lwd = 2, lty = 1)
  lines(loc_dens, col = "green", lwd = 2, lty = 2)
  lines(x_seq, true_dens, col = "purple", lwd = 2, lty = 3)
  

  abline(v = true_param, col = "red", lty = 4, lwd = 1.5)
  

  legend("topright", 
         legend = c("ABC Rejection", "Local Linear Regression", "True Posterior", "True Value"),
         col = c("blue", "green", "purple", "red"), 
         lty = c(1, 2, 3, 4),
         lwd = c(2, 2, 2, 1.5))
}


result_summary <- data.frame(
  Tolerance = rep(tolerance_levels, each = 3),
  Method = rep(c("ABC Rejection", "Local Linear Regression", "True Posterior"), times = length(tolerance_levels)),
  Posterior_Mean = c(
    # Tolerance level 0.5
    mean(rejection_params[[1]]), 
    mean(loclinear_params[[1]]), 
    posterior_mean,
    # Tolerance level 0.3
    mean(rejection_params[[2]]), 
    mean(loclinear_params[[2]]), 
    posterior_mean,
    # Tolerance level 0.1
    mean(rejection_params[[3]]), 
    mean(loclinear_params[[3]]), 
    posterior_mean
  ),
  Posterior_SD = c(
    # Tolerance level 0.5
    sd(rejection_params[[1]]), 
    sd(loclinear_params[[1]]), 
    posterior_sd,
    # Tolerance level 0.3
    sd(rejection_params[[2]]), 
    sd(loclinear_params[[2]]), 
    posterior_sd,
    # Tolerance level 0.1
    sd(rejection_params[[3]]), 
    sd(loclinear_params[[3]]), 
    posterior_sd
  )
)

# Print results organized by tolerance level
for(tol in tolerance_levels) {
  cat("\n--- Results for Tolerance =", tol, "---\n")
  subset_results <- result_summary[result_summary$Tolerance == tol, ]
  print(subset_results[, c("Method", "Posterior_Mean", "Posterior_SD")])
}