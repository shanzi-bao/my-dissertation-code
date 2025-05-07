library(progress)
library(mvtnorm)  
library(Matrix)
library(corpcor)
library(SpatialExtremes)
library(MASS)  
library(Rcpp)

sourceCpp('summary_statistics_calculation.cpp')



simData <- function(theta_row) {
  cov11 <- theta_row[1]
  cov12 <- theta_row[2]
  cov22 <- theta_row[3]
  beta.loc <- theta_row[4:6]
  beta.scale <- theta_row[7:9]
  shape <- theta_row[10]
  

  cov_matrix <- matrix(c(cov11, cov12, cov12, cov22), ncol = 2)
  cov_matrix <- cov_matrix + diag(1e-6, nrow = 2) 
  

  determinant_value <- cov11 * cov22 - cov12^2
  

  if (determinant_value <= 0) {
    
    return(NA)
  }
  

  if (!is.positive.definite(cov_matrix)) {

    return(NA)
  }

  result <- simData_spExt(47, 79, sim_coords, cov11, cov12, cov22, beta.loc, beta.scale, shape)

  return(result)
}


compute_distances_and_weights <- function(coords, d) {

  n_sites <- nrow(coords)
  

  dist_matrix <- matrix(0, nrow = n_sites, ncol = n_sites)
  weights <- matrix(0, nrow = n_sites, ncol = n_sites)

  for (i in 1:(n_sites - 1)) {
    for (j in (i + 1):n_sites) {
    
      dist_matrix[i, j] <- sqrt(sum((coords[i, ] - coords[j, ])^2))
      dist_matrix[j, i] <- dist_matrix[i, j] 
      

      if (dist_matrix[i, j] < d) {
        weights[i, j] <- 1
        weights[j, i] <- 1 
      }
    }
  }
  
  return(list(dist_matrix = dist_matrix, weights = weights))
}
