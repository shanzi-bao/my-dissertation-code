#include <Rcpp.h>
using namespace Rcpp;


struct SitePair {
  int site_i;
  int site_j;
  double distance;
};

// log bivariate density function
double grplik_singlPairR_cpp(NumericVector params, NumericMatrix data, NumericMatrix coords, 
                             int obs_index, int site_index_i, int site_index_j) {

  double cov11 = params[0];
  double cov12 = params[1];
  double cov22 = params[2];
  
  NumericVector beta_mu = params[Range(3, 5)];
  NumericVector beta_la = params[Range(6, 8)];
  double xi = params[9];
  

  double y_i = data(obs_index, site_index_i);
  double y_j = data(obs_index, site_index_j);
  
  NumericVector coord_i = coords.row(site_index_i);
  NumericVector coord_j = coords.row(site_index_j);

  NumericVector h = coord_j - coord_i;
  
  //  a(h)
  double idet = 1.0 / (cov11 * cov22 - cov12 * cov12);
  double a = sqrt((cov11 * h[1] * h[1] - 2 * cov12 * h[0] * h[1] + cov22 * h[0] * h[0]) * idet);
  
  // z_i 和 z_j
  NumericVector X_i = NumericVector::create(1.0, coord_i[0], coord_i[1]);
  NumericVector X_j = NumericVector::create(1.0, coord_j[0], coord_j[1]);
  
  double Xbeta_mui = sum(X_i * beta_mu);
  double Xbeta_muj = sum(X_j * beta_mu);
  double Xbeta_lai = sum(X_i * beta_la);
  double Xbeta_laj = sum(X_j * beta_la);
  

  double lambda_i = (Xbeta_lai);
  double lambda_j = (Xbeta_laj);
  

  if (lambda_i <= 0.0 || lambda_j <= 0.0) {
    return R_NegInf; 
  }
  
 
  double zero_i = std::max(1.0 + xi * (y_i - Xbeta_mui) / lambda_i, 0.0);
  double zero_j = std::max(1.0 + xi * (y_j - Xbeta_muj) / lambda_j, 0.0);
  
  double z_i = pow(zero_i, 1.0 / xi);
  double z_j = pow(zero_j, 1.0 / xi);
  

  double w = a * 0.5 + log(z_j / z_i) / a;
  double v = a - w;

  double Phiw = R::pnorm(w, 0.0, 1.0, true, false);
  double Phiv = R::pnorm(v, 0.0, 1.0, true, false);
  double phiw = R::dnorm(w, 0.0, 1.0, false);
  double phiv = R::dnorm(v, 0.0, 1.0, false);
  

  double A = -Phiw/z_i - Phiv/z_j;
  

  double z_i2 = z_i * z_i;
  double z_j2 = z_j * z_j;
  double B = Phiw / z_i2 + phiw / (z_i2 * a) - phiv / (a * z_j * z_i);
  double C = Phiv / z_j2 + phiv / (z_j2 * a) - phiw / (a * z_i * z_j);
  

  double D = v * phiw / (a * a * z_i2 * z_j) + w * phiv / (a * a * z_j2 * z_i);
  

  double BCpD = B * C + D;
  

  double E = log(1.0 / (Xbeta_lai * Xbeta_laj)) + 
    log(std::max(1.0 + xi * (y_i - Xbeta_mui) / Xbeta_lai, 0.0) / Xbeta_lai) * (1.0 / xi - 1) +
    log(std::max(1.0 + xi * (y_j - Xbeta_muj) / Xbeta_laj, 0.0) / Xbeta_laj) * (1.0 / xi - 1);

  double log_likelihood = A + log(BCpD) + E;

  return log_likelihood;
}



NumericVector compute_partial_gradient(NumericVector params, NumericMatrix data, NumericMatrix coords, 
                                       int obs_index, int site_index_i, int site_index_j, 
                                       IntegerVector target_indices) {
  
  double h = sqrt(std::numeric_limits<double>::epsilon());
  int target_size = target_indices.size();
  NumericVector partial_gradient(target_size);
  
  for (int idx = 0; idx < target_size; idx++) {
    int i = target_indices[idx];
    NumericVector ei(params.size(), 0.0);
    ei[i] = 1.0;
    
    double log_likelihood_forward = grplik_singlPairR_cpp(params + h * ei, data, coords, 
                                                          obs_index, site_index_i, site_index_j);
    double log_likelihood_backward = grplik_singlPairR_cpp(params - h * ei, data, coords, 
                                                           obs_index, site_index_i, site_index_j);
    
    partial_gradient[idx] = (log_likelihood_forward - log_likelihood_backward) / (2 * h);

    
  }
  
  return partial_gradient;
}

// [[Rcpp::export]]

NumericMatrix compute_distance_based_gradient_q(
    NumericVector params,
    NumericMatrix data,
    NumericMatrix coords,
    NumericMatrix dist_matrix,
    int m,                 
    int k,                 
    IntegerVector target_indices) {
  

  std::vector<double> distances;
  for (int i = 0; i < dist_matrix.nrow(); i++) {
    for (int j = i + 1; j < dist_matrix.ncol(); j++) {
      distances.push_back(dist_matrix(i, j));
    }
  }
  

  std::sort(distances.begin(), distances.end());
  

  std::vector<double> quantile_breaks(m + 1);
  quantile_breaks[0] = 0.0; 
  
  for (int i = 1; i < m; i++) {
    int index = (i * distances.size()) / m;
    quantile_breaks[i] = distances[index];
  }
  
  quantile_breaks[m] = distances.back() + 1e-10; 
  

  std::vector<std::vector<SitePair>> distance_groups(k);
  

  for (int i = 0; i < dist_matrix.nrow(); i++) {
    for (int j = i + 1; j < dist_matrix.ncol(); j++) {
      double dist = dist_matrix(i, j);
      

      int group_idx = 0;
      for (int q = 1; q <= m; q++) {
        if (dist < quantile_breaks[q]) {
          group_idx = q - 1;
          break;
        }
      }
      

      if (group_idx < k) {
        SitePair pair = {i, j, dist};
        distance_groups[group_idx].push_back(pair);
      }
    }
  }
  

  int p = target_indices.size();
  NumericMatrix batched_gradients_sum(k, p);
  

  for (int obs = 0; obs < data.nrow(); obs++) {

    for (int batch = 0; batch < k; batch++) {
   
      for (size_t pair_idx = 0; pair_idx < distance_groups[batch].size(); pair_idx++) {
        SitePair& pair = distance_groups[batch][pair_idx];
        
        NumericVector partial_gradient = compute_partial_gradient(
          params, data, coords, obs, pair.site_i, pair.site_j, target_indices
        );
        

        for (int p_idx = 0; p_idx < p; p_idx++) {
          batched_gradients_sum(batch, p_idx) += partial_gradient[p_idx];
        }
      }
    }
  }
  

  int n_obs = data.nrow();
  if (n_obs > 0) {
    for (int batch = 0; batch < k; batch++) {
      for (int p_idx = 0; p_idx < p; p_idx++) {
        batched_gradients_sum(batch, p_idx) /= n_obs;
      }
    }
  }
  
  return batched_gradients_sum;
}

