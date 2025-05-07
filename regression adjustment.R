# === ABC Regression Adjustment Function ===
abc_regression_adjustment_v2 <- function(theta_samples, 
                                         summary_stats_list, 
                                         observed_summary_matrix,
                                         opt_indices = NULL,
                                         true_params = NULL,
                                         bandwidth = NULL) {
  
  # Input validation and diagnostics
  cat("Checking regression inputs:\n")
  cat(sprintf("Number of samples: %d\n", nrow(theta_samples)))
  cat(sprintf("Number of summary statistics per group: %d\n", nrow(observed_summary_matrix)))
  
  n_samples <- nrow(theta_samples)
  n_dist_groups <- nrow(observed_summary_matrix)
  n_params_opt <- ncol(observed_summary_matrix)
  
  if (is.null(opt_indices)) {
    opt_indices <- 1:ncol(theta_samples)
  }
  
  theta_adj <- theta_samples  # Initialize adjusted result matrix
  
  # Compute distances between simulated summaries and observed summary
  distances <- numeric(n_samples)
  for (i in 1:n_samples) {
    diff_matrix <- as.matrix(summary_stats_list[[i]] - observed_summary_matrix)
    distances[i] <- sqrt(sum(diff_matrix^2))
  }
  
  # Normalize distances
  scaled_distances <- scale(distances)
  
  # Estimate bandwidth if not provided
  if (is.null(bandwidth)) {
    bandwidth <- quantile(abs(scaled_distances), probs = 1)
  }
  
  # Compute Epanechnikov kernel weights
  kernel_weights <- ifelse(abs(scaled_distances)/bandwidth <= 1,
                           3/4 * bandwidth^(-1) * (1 - (abs(scaled_distances)/bandwidth)^2),
                           0)
  
  cat(sprintf("\nBandwidth: %.6f\n", bandwidth))
  cat(sprintf("Non-zero weights: %d\n", sum(kernel_weights > 0)))
  cat(sprintf("Mean kernel weight: %.6f\n", mean(kernel_weights)))
  cat(sprintf("Weight range: [%.6f, %.6f]\n", 
              min(kernel_weights[kernel_weights > 0]), max(kernel_weights)))
  
  # Perform weighted regression adjustment for each parameter
  for (i in seq_along(opt_indices)) {
    j <- opt_indices[i]
    
    cat(sprintf("\nProcessing parameter %d (original index %d):\n", i, j))
    
    X <- matrix(0, nrow = n_samples, ncol = n_dist_groups)
    for (k in 1:n_samples) {
      diff <- summary_stats_list[[k]][, i] - observed_summary_matrix[, i]
      X[k, ] <- scale(diff)
    }
    
    model <- lm(theta_samples[, j] ~ X, weights = kernel_weights)
    
    cat(sprintf("R-squared: %.4f\n", summary(model)$r.squared))
    
    adjustment <- predict(model) - model$coefficients[1]
    cat(sprintf("Adjustment range: [%.4f, %.4f]\n", min(adjustment), max(adjustment)))
    
    theta_adj[, j] <- theta_samples[, j] - adjustment
  }
  
  return(theta_adj)
}


# === Plotting Function for Posterior Densities ===
plot_adjustment_v2 <- function(theta_original, theta_adjusted, opt_indices, true_params = NULL) {
  par(mfrow = c(2, 2))
  
  for (i in seq_along(opt_indices)) {
    param_idx <- opt_indices[i]
    
    d_orig <- density(theta_original[, param_idx])
    d_adj <- density(theta_adjusted[, param_idx])
    
    ylim <- range(c(d_orig$y, d_adj$y))
    
    plot(d_orig, main = paste("Parameter", param_idx),
         col = "gray", ylim = ylim, lwd = 2,
         xlab = "Parameter Value", ylab = "Density")
    lines(d_adj, col = "red", lwd = 2)
    
    if (!is.null(true_params)) {
      abline(v = true_params[param_idx], lty = 2, col = "black", lwd = 2)
    }
    
    legend("topright", 
           legend = c("Original", "Adjusted", if (!is.null(true_params)) "True Value"),
           col = c("gray", "red", if (!is.null(true_params)) "black"),
           lty = c(1, 1, if (!is.null(true_params)) 2),
           lwd = 2)
  }
  
  par(mfrow = c(1, 1))
}


# === Summary Statistics Printer ===
print_adjustment_summary_v2 <- function(theta_original, theta_adjusted, opt_indices, true_params = NULL) {
  cat("=== Mean and Standard Deviation Comparison ===\n")
  for (i in seq_along(opt_indices)) {
    param_idx <- opt_indices[i]
    cat(sprintf("\nParameter %d:\n", param_idx))
    
    if (!is.null(true_params)) {
      cat(sprintf("True value: %.4f\n", true_params[param_idx]))
    }
    
    cat(sprintf("Original mean (sd): %.4f (%.4f)\n",
                mean(theta_original[, param_idx]),
                sd(theta_original[, param_idx])))
    
    cat(sprintf("Adjusted mean (sd): %.4f (%.4f)\n",
                mean(theta_adjusted[, param_idx]),
                sd(theta_adjusted[, param_idx])))
    
    if (!is.null(true_params)) {
      orig_dist <- abs(mean(theta_original[, param_idx]) - true_params[param_idx])
      adj_dist <- abs(mean(theta_adjusted[, param_idx]) - true_params[param_idx])
      improvement <- (orig_dist - adj_dist) / orig_dist * 100
      cat(sprintf("Distance improvement: %.2f%%\n", improvement))
    }
  }
}


# === Boxplot Comparison Function ===
plot_boxplot_comparison <- function(theta_original, theta_adjusted, opt_indices, true_params = NULL) {
  par(mfrow = c(ceiling(length(opt_indices) / 2), 2))
  
  for (i in seq_along(opt_indices)) {
    param_idx <- opt_indices[i]
    
    data <- data.frame(
      Value = c(theta_original[, param_idx], theta_adjusted[, param_idx]),
      Group = rep(c("Original", "Adjusted"), each = nrow(theta_original))
    )
    
    boxplot(Value ~ Group, data = data, main = paste("Parameter", param_idx),
            xlab = "Group", ylab = "Parameter Value",
            col = c("gray", "red"), border = "black", outline = TRUE)
    
    if (!is.null(true_params)) {
      abline(h = true_params[param_idx], col = "blue", lty = 2, lwd = 2)
    }
  }
  
  par(mfrow = c(1, 1))
}


# === Batch Processing Function ===
process_abc_adjustments <- function(
    result_names,
    results_list,
    true_params,
    plot = TRUE,
    print_summary = TRUE
) {
  if (!is.vector(result_names) || !is.character(result_names)) stop("result_names must be a character vector")
  if (!is.list(results_list)) stop("results_list must be a list")
  if (!all(result_names %in% names(results_list))) stop("All result_names must exist in results_list")
  
  adjusted_results <- list()
  
  for (res in result_names) {
    cat(sprintf("\nProcessing %s...\n", res))
    
    current_result <- results_list[[res]]
    
    adjusted_results[[res]] <- abc_regression_adjustment_v2(
      theta_samples = current_result$final_theta,
      summary_stats_list = current_result$final_sim_scores,
      observed_summary_matrix = current_result$observed_scores,
      opt_indices = current_result$opt_indices,
      true_params = true_params
    )
  }
  
  cat("\nAll adjustments completed!\n")
  
  if (plot) {
    for (res in result_names) {
      cat(sprintf("\nPlotting %s...\n", res))
      
      current_result <- results_list[[res]]
      
      plot_adjustment_v2(
        theta_original = current_result$final_theta,
        theta_adjusted = adjusted_results[[res]],
        opt_indices = current_result$opt_indices,
        true_params = true_params
      )
    }
  }
  
  if (print_summary) {
    for (res in result_names) {
      cat(sprintf("\nSummary for %s:\n", res))
      
      current_result <- results_list[[res]]
      
      print_adjustment_summary_v2(
        theta_original = current_result$final_theta,
        theta_adjusted = adjusted_results[[res]],
        opt_indices = current_result$opt_indices,
        true_params = true_params
      )
    }
  }
  
  return(adjusted_results)
}


adjusted_theta <- abc_regression_adjustment_v2(
  theta_samples = result$final_theta,
  summary_stats_list = result$final_sim_scores,
  observed_summary_matrix = result$observed_scores,
  opt_indices = result$opt_indices,
  true_params = hat.param
)

plot_adjustment_v2(
  theta_original = result$final_theta,
  theta_adjusted = adjusted_theta,
  opt_indices = result$opt_indices,
  true_params = hat.param
)

print_adjustment_summary_v2(
  theta_original = result$final_theta,
  theta_adjusted = adjusted_theta,
  opt_indices = result$opt_indices,
  true_params = hat.param
)
