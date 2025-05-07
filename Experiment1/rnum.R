source("supplementary_function.R")
source("observed_data.R")

APMC_adj_q <- function(N, alf, paccmin, sim_coords,  observed_data, param,  
                       m,  
                       k,  
                       fixed_params = NULL, df = 5) {
  

  final_theta <- NULL      
  final_sim_scores <- NULL
  final_distances <- NULL  
 
  hatCholJ <- solve(t(chol(mod.spe1$var.score)))
  
  w_raw <- rep(1, N)
  diagnostics_data <- list()
  

  if (!is.numeric(N) || N <= 0 || N != floor(N)) 
    stop("N must be a positive integer")
  if (!is.numeric(alf) || alf <= 0 || alf >= 1) 
    stop("alf must be between 0 and 1")
  if (!is.numeric(paccmin) || paccmin <= 0 || paccmin >= 1) 
    stop("paccmin must be between 0 and 1")
  if (!is.matrix(sim_coords)) 
    stop("sim_coords must be a matrix")
  
  all_params <- 1:length(param$mean)
  fixed_indices <- as.numeric(names(fixed_params))
  opt_indices <- setdiff(all_params, fixed_indices)
  

  all_params <- 1:length(param$mean)
  if (!is.null(fixed_params)) {
    fixed_indices <- as.numeric(names(fixed_params))
    opt_indices <- as.integer(setdiff(all_params, fixed_indices))  
    cpp_indices <- opt_indices - 1L  
  } else {
    opt_indices <- as.integer(all_params) 
    cpp_indices <- opt_indices - 1L
  }
  
  
  
  make.positive.definite <- function(matrix, tol = 1e-6) {
    eigen_decomp <- eigen(matrix)
    values <- eigen_decomp$values
    vectors <- eigen_decomp$vectors
    values[values < tol] <- tol
    result <- vectors %*% diag(values) %*% t(vectors)
    result <- (result + t(result)) / 2
    return(result)
  }
  
  
 
  generate_new_particles <- function(thetao, sigma, param_scales, N, df = 5) {

    thetao_subset <- thetao[, opt_indices, drop = FALSE]
    

    thetao_scaled <- sweep(thetao_subset, 2, param_scales, "/")
    
    
    thetan_scaled <- t(sapply(1:nrow(thetao_scaled), function(i) {
      rmvt(1, sigma = 2*sigma, df = df, delta = thetao_scaled[i,])
    }))
    

    thetan <- matrix(0, nrow = nrow(thetan_scaled), ncol = length(param$mean))
    thetan[, opt_indices] <- sweep(thetan_scaled, 2, param_scales, "*")

    if (!is.null(fixed_params)) {
      for (idx in names(fixed_params)) {
        thetan[, as.numeric(idx)] <- fixed_params[[idx]]
      }
    }
    
    return(thetan)
  }
  
  
  
  
  

  compute_log_weights <- function(theta_new, theta_old, w_old, sigma, param_scales, df = 5) {    
    
    log_prior_new <- log_uniform_prior(theta_new)
    if (!is.finite(log_prior_new)) return(0)
    

    theta_new_scaled <- theta_new[opt_indices] / param_scales
    theta_old_scaled <- sweep(theta_old[, opt_indices, drop = FALSE], 2, param_scales, "/")
    
    w_normalized <- w_old / sum(w_old)
    

    log_densities <- numeric(nrow(theta_old_scaled))
    for(j in 1:nrow(theta_old_scaled)) {
      log_densities[j] <- dmvt(
        x     = theta_new_scaled,
        delta = theta_old_scaled[j,],
        sigma = 2 * sigma,
        df    = df,
        log   = TRUE
      )
    }
    
    log_weight <- log_prior_new - log(sum(w_normalized * exp(log_densities)))
    return(exp(log_weight))  
  }
  
  log_uniform_prior <- function(theta) {

    if(length(theta) == length(opt_indices)) {

      full_theta <- numeric(10)
      full_theta[opt_indices] <- theta
      if (!is.null(fixed_params)) {
        for (idx in names(fixed_params)) {
          full_theta[as.numeric(idx)] <- fixed_params[[idx]]
        }
      }
    } else {

      full_theta <- theta
    }
    

    if (full_theta[1] <= 0     || full_theta[1] >= 1000) return(-Inf)
    if (full_theta[2] <= -300  || full_theta[2] >= 300)  return(-Inf)
    if (full_theta[3] <= 0     || full_theta[3] >= 1000) return(-Inf)
    if (abs(full_theta[2]) >= sqrt(full_theta[1] * full_theta[3])) return(-Inf)
    if (full_theta[10] <= 0)                             return(-Inf)
    
    return(-log(1000) - log(600) - log(1000))
  }

  p <- round(N * alf)  
  max_iter <- 30 
  

  EPS <- numeric(max_iter)
  ACC <- numeric(max_iter)
  THETA <- vector("list", max_iter)
  W <- vector("list", max_iter)
  COUNTSIM <- numeric(max_iter)
  

  pb <- progress_bar$new(
    format = "  progress [:bar] :percent | time: :elapsed | estimated rest time: :eta",
    total = 4,
    clear = FALSE,
    width = 60
  )
  
 

  
  theta_init <- matrix(0, nrow = N, ncol = length(param$mean))
  

  if (!is.null(fixed_params)) {
    for (idx in names(fixed_params)) {
      theta_init[, as.numeric(idx)] <- fixed_params[[idx]]
    }
  }
  
  
  

  theta_list <- replicate(N, {
    repeat {
      theta <- rmvt(1, 
                    sigma = param$vcov[opt_indices, opt_indices, drop = FALSE], 
                    df = df, 
                    delta = param$mean[opt_indices])
      if(is.finite(log_uniform_prior(theta))) break
    }
    theta
  }, simplify = FALSE)
  
  theta_matrix <- do.call(rbind, theta_list)
  theta_init[, opt_indices] <- theta_matrix
  
  
  
  theta <- theta_init
  w <- rep(1, N)
  pb$tick()
  
  
  tryCatch({
    if(any(is.na(sim_coords))) stop("sim_coords contains NA values")

    dist_matrix <- as.matrix(dist(sim_coords))
    
   
    observed_batched_scores <- compute_distance_based_gradient_q(
      params = param$mean,
      data = observed_data,
      coords = sim_coords,
      dist_matrix = dist_matrix,
      m = m,  
      k = k, 
      target_indices = cpp_indices
    )
    

    if(is.null(observed_batched_scores)) stop("compute_adaptive_batched_gradient returned NULL")
    if(any(is.na(observed_batched_scores))) stop("observed_batched_scores contains NA values")
    
    score_scale <- apply(observed_batched_scores, 2, function(x) {
      if(any(is.na(x))) return(1)
      s <- sd(x, na.rm = TRUE)
      if(is.na(s) || s == 0) return(1)
      return(s)
    })
    
    message("Score scales computed successfully")
    
  }, error = function(e) {
    stop(sprintf(" observed summary statistics error: %s", e$message))
  })
  pb$tick()

  x <- lapply(1:N, function(i) {
    
    
    
    tryCatch({
      sim_data <- simData(theta[i, ])
      if(is.null(sim_data) || any(is.na(sim_data))) {
        message(sprintf(" %d simulation invalid", i))
        return(NULL)
      }
      return(sim_data)
    }, error = function(e) {
      message(sprintf(" %d simulation invalid: %s", i, e$message))
      return(NULL)
    })
  })
  
  sim_scores_list <- lapply(x, function(sim_data) {
    if(is.null(sim_data)) return(NULL)
    tryCatch({
      compute_distance_based_gradient_q(
        params = param$mean,
        data = sim_data,
        coords = sim_coords,
        dist_matrix = dist_matrix,
        m = m,  
        k = k,  
        target_indices = cpp_indices
      )
    }, error = function(e) {

      return(NULL)
    })
  })
  scores_old <- sim_scores_list
  gc()
  pb$tick()
  

  score_scales <- apply(observed_batched_scores, 2, function(col) {
    sd_val <- sd(col, na.rm = TRUE)
    if (is.na(sd_val) || sd_val < 1e-10) return(1)
    return(sd_val)
  })
  
  

  distances <- sapply(sim_scores_list, function(sim_score) {
    if(is.null(sim_score) || !all(is.finite(sim_score))) {
      return(Inf)
    }
    

    total_diff <- 0
    for(i in 1:nrow(sim_score)) {

      obs_row <- as.numeric(observed_batched_scores[i,])
      sim_row <- as.numeric(sim_score[i,])
      
  
      obs_row <- matrix(obs_row, ncol=1)
      sim_row <- matrix(sim_row, ncol=1)
      
   
      hatCholJ_subset <- hatCholJ[opt_indices, opt_indices]

      whitened_obs <- hatCholJ_subset %*% obs_row
      whitened_sim <- hatCholJ_subset %*% sim_row
      
   
      diff_vec <- whitened_obs - whitened_sim
      total_diff <- total_diff + sum(diff_vec^2)
    }
    
    return(sqrt(total_diff))
  })
  

  

  valid_indices <- which(!is.na(distances)) 
  if (length(valid_indices) < p) {
    stop("no enough particles")
  }
  

  distances <- distances[valid_indices]
  theta <- theta[valid_indices, , drop = FALSE]
  x <- x[valid_indices]
  w <- w[valid_indices]
  
  
  order_indices <- order(distances)
  theta <- theta[order_indices[1:p], , drop = FALSE]
  x <- x[order_indices[1:p]]
  distances <- distances[order_indices[1:p]]
  w <- w[order_indices[1:p]]
  

  param_scales <- apply(theta[, opt_indices, drop = FALSE], 2, sd)
  
  param_scales[param_scales < 1e-10] <- 1
  

  theta_scaled <- sweep(theta[, opt_indices, drop = FALSE], 2, param_scales, "/")
  sigma <- make.positive.definite(cov.wt(theta_scaled, wt = w/sum(w))$cov)
  sigma <- sigma + diag(0.01, ncol(sigma))
  
  eps <- distances[p]
  current_iter <- 1
  EPS[current_iter] <- eps
  THETA[[current_iter]] <- theta
  W[[current_iter]] <- w
  ACC[current_iter] <- 1
  COUNTSIM[current_iter] <- N
  
  pb$tick()
  

  main_pb <- progress_bar$new(
    format = "  main iteration [:bar] :percent | time used: :elapsed | estimated rest time: :eta",
    total = max_iter,
    clear = FALSE,
    width = 60
  )
  

  sd_old <- NULL
  sd_tol <- 0.001
  

  while (current_iter <= max_iter && ACC[current_iter] > paccmin) {
    current_iter <- current_iter + 1
    message(sprintf("\n========== %d iteration ==========", current_iter - 1))
    
    t <- sample(length(w), N-p, replace = TRUE, prob = w)
    xo <- x[t]
    thetao <- theta[t, , drop = FALSE]
    distances0 <- distances[t]

    tryCatch({
      thetan <- generate_new_particles(thetao, sigma, param_scales, N - p)
      
    }, error = function(e) {
      stop(sprintf("error: %s", e$message))
    })
    
    
    # 生成新的模拟数据
    xn <- lapply(1:(N - p), function(i) {
      tryCatch({
        sim_data <- simData(thetan[i, ])
        if(is.null(sim_data) || any(is.na(sim_data))) {
  
          return(NULL)
        }
        return(sim_data)
      }, error = function(e) {
     
        return(NULL)
      })
    })
    

    sim_scores_new <- lapply(xn, function(sim_data) {
      if(is.null(sim_data)) return(NULL)
      tryCatch({
        compute_distance_based_gradient_q(
          params = param$mean,
          data = sim_data,
          coords = sim_coords,
          dist_matrix = dist_matrix,
          m = m,
          k = k,
          target_indices = cpp_indices
        )
      }, error = function(e) {
      
        return(NULL)
      })
    })

    distances_new <- sapply(sim_scores_new, function(sim_score) {
      if (is.null(sim_score) || !all(is.finite(sim_score))) {
        return(Inf)
      }
      
   
      total_diff <- 0
      for(i in 1:nrow(sim_score)) {

        obs_row <- as.numeric(observed_batched_scores[i,])
        sim_row <- as.numeric(sim_score[i,])
        

        obs_row <- matrix(obs_row, ncol=1)
        sim_row <- matrix(sim_row, ncol=1)
       
        hatCholJ_subset <- hatCholJ[opt_indices, opt_indices]
        

        whitened_obs <- hatCholJ_subset %*% obs_row
        whitened_sim <- hatCholJ_subset %*% sim_row
        

        diff_vec <- whitened_obs - whitened_sim
        total_diff <- total_diff + sum(diff_vec^2)
      }
      
      return(sqrt(total_diff))
    })
    

    valid_before_eps <- !(distances_new == Inf | is.na(distances_new) | 
                            is.null(distances_new) | !is.finite(distances_new))
 
    

    invalid_particles <- !valid_before_eps | distances_new > eps
    

    invalid_particles <- distances_new == Inf | is.na(distances_new) | 
      is.null(distances_new) | !is.finite(distances_new)| distances_new > eps

    wt_new <- numeric(N - p)
    for(i in 1:(N - p)) {
      if(!invalid_particles[i]) {
        wt_new[i] <- compute_log_weights(thetan[i,], theta, w, sigma, param_scales)
      }
    }
    
    thetat <- rbind(theta, matrix(0, nrow = N-p, ncol = ncol(theta)))
    for(i in 1:(N-p)) {
      if(invalid_particles[i]) {
        thetat[p+i,] <- thetao[i,]
      } else {
        thetat[p+i,] <- thetan[i,]
      }
    }
    

    xt <- c(x, vector("list", N-p))
    for(i in 1:(N-p)) {
      if(invalid_particles[i]) {
        xt[[p+i]] <- xo[[i]]
      } else {
        xt[[p+i]] <- xn[[i]]
      }
    }
    message("Debug: Types of distances:")
    

    scores_combined <- c(scores_old, vector("list", N-p))
    for(i in 1:(N-p)) {
      if(invalid_particles[i]) {
        scores_combined[[p+i]] <- scores_old[[t[i]]]  
      } else {
        scores_combined[[p+i]] <- sim_scores_new[[i]] 
      }
    }
    
 
    distances_combined <- numeric(length(distances0))
    
    for (i in seq_along(invalid_particles)) {
      if (invalid_particles[i]) {
        distances_combined[i] <- distances0[i]
      } else {
        distances_combined[i] <- distances_new[i]
      }
    }
    
   
    distancest <- c(distances, distances_combined)
    
  
    wt <- c(w, wt_new)
    
    weight_indices <- which(wt > 0)
    
    
    valid_indices <- weight_indices[!sapply(scores_combined[weight_indices], is.null)]
    
    
    if(length(valid_indices) < p) {
  
      break
    }
    
    temp_order <- order(distancest[valid_indices])
    selected_indices <- valid_indices[temp_order[1:p]]
    
    
    
    theta <- thetat[selected_indices, , drop = FALSE]
    x <- xt[selected_indices]
    distances <- distancest[selected_indices]
    w <- wt[selected_indices]
    

    scores_old <- scores_combined[selected_indices]

    

    eps <- distances[p]
    param_scales <- apply(theta[, opt_indices, drop = FALSE], 2, sd)
    param_scales[param_scales < 1e-10] <- 1
    
  
    theta_scaled <- sweep(theta[, opt_indices, drop = FALSE], 2, param_scales, "/")
    sigma <- make.positive.definite(cov.wt(theta_scaled, wt = w/sum(w))$cov)
    


    EPS[current_iter] <- eps
    THETA[[current_iter]] <- theta
    W[[current_iter]] <- w
    accepted_new_particles <- sum(selected_indices > p)
    ACC[current_iter] <- accepted_new_particles / (N - p)
    COUNTSIM[current_iter] <- COUNTSIM[current_iter - 1] + (N - p)
 
    final_theta <- theta
    final_sim_scores <- scores_old  

    message(sprintf("acceptance rate: %.4f", ACC[current_iter]))
    message(sprintf("eps: %.4f", eps))
    message(sprintf("valid particles numbers: %d", sum(!invalid_particles)))
    message("mean:")
    print(colMeans(theta[, opt_indices, drop = FALSE]))
    message("sd:")
    print(apply(theta[, opt_indices, drop = FALSE], 2, sd))
    

    sd_new <- apply(theta[, opt_indices, drop = FALSE], 2, sd)
    if (!is.null(sd_old)) {
      rel_change <- max(abs(sd_new - sd_old) / sd_old)
      if (rel_change < sd_tol) {
        message("sd converges, end")
        break
      }
    }
    sd_old <- sd_new
    
 
    diagnostics_data[[current_iter]] <- list(
      iteration = current_iter,
      paramSD = apply(theta[, opt_indices, drop = FALSE], 2, sd),
      paramMean = colMeans(theta[, opt_indices, drop = FALSE]),
      weightMax = max(w),
      weightMin = min(w[w > 0]),
      ess = 1/sum((w/sum(w))^2),
      acceptance = ACC[current_iter],
      tolerance = eps,
      sigmaCondNum = kappa(sigma),
      sigmaEigenvals = eigen(sigma)$values
    )
    

    if(current_iter %% 5 == 0) {
      cat("\n==== iteration", current_iter, "information ====\n")
      cat("sd change rate:", 
          if(current_iter > 1) {
            (diagnostics_data[[current_iter]]$paramSD - 
               diagnostics_data[[max(1, current_iter-1)]]$paramSD) / 
              diagnostics_data[[max(1, current_iter-1)]]$paramSD
          } else {
            NA
          }, "\n")
      cat("ESS:", diagnostics_data[[current_iter]]$ess, "\n")
      cat("acceptance:", diagnostics_data[[current_iter]]$acceptance, "\n")
   
    }
    
    gc()
    main_pb$tick()
  }
  

  valid_iterations <- 1:current_iter
  gc()
  
  

  return(list(
    THETA = THETA[valid_iterations],
    ACC = ACC[valid_iterations],
    W = W[valid_iterations],
    EPS = EPS[valid_iterations],
    COUNTSIM = COUNTSIM[valid_iterations],
    diagnostics = diagnostics_data,
    fixed_params = fixed_params,
    opt_indices = opt_indices,
    final_theta = final_theta,
    final_sim_scores = final_sim_scores,
    final_distances = final_distances,
    observed_scores = observed_batched_scores
    
  ))
}


fixed_params <- list(
  # "1" = 332,       # Fix parameter 1 to 332
  # "2" = 70,        # Fix parameter 2 to 70
  # "3" = 185,       # Fix parameter 3 to 185
  "4" = 20.65,      # Fix parameter 4 to 20.65
  "5" = 0.064,      # Fix parameter 5 to 0.064
  "6" = -0.155,     # Fix parameter 6 to -0.155
  "7" = 3.54,       # Fix parameter 7 to 3.54
  "8" = 0.023,      # Fix parameter 8 to 0.023
  "9" = -0.039,     # Fix parameter 9 to -0.039
  "10" = 0.192      # Fix parameter 10 to 0.192
)


num102 <- APMC_adj_q(
  N = 1000,
  alf = 0.3,
  paccmin = 0.01,
  sim_coords = sim_coords,
  
  observed_data = observed_data,
  param = is.param,
  m = 10,
  k = 2,  
  fixed_params = fixed_params
)

