functions {
  real sparse_car_lpdf(vector phi, 
  // real tau,
  real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] +
        phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] +
        phi[W_sparse[i, 1]];
      }
      
    
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (sum(ldet_terms) - 
      (phit_D * phi - alpha * (phit_W * phi)));
  }
  real corr_val(vector x, vector y, int n){
    int n_vec;
    real num;
    real den_1;
    real den_2;
    n_vec = n;
    num = n_vec*sum(x .* y) - sum(x)*sum(y);
    den_1 = sqrt(n_vec*sum(x .* x) - pow(sum(x), 2));
    den_2 = sqrt(n_vec*sum(y .* y) - pow(sum(y), 2));
    return num/(den_1 * den_2);
    }
}
data {
  int<lower = 1> n;
  int<lower = 1> m;
  // Outcomes
  int<lower = 0> y1[n];
  int<lower = 0> y2[m];
  vector[n] log_offset1;
  vector[m] log_offset2;
  // Spatial matrices
  matrix<lower = 0, upper = 1>[n, n] W; // adjacency matrix
  int W_n; // number of adjacent region pairs
  // Multiple membership
  matrix[m,n] M_W;
}
transformed data {
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[n] D_sparse;     // diagonal of D (no of neigbors for each site)
  vector[n] lambda;       // eigenvalues of invsqrtD * W * invsqrtD
  // Adjacency
  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:n) D_sparse[i] = sum(W[i]);
  {
    vector[n] invsqrtD;  
    for (i in 1:n) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
}
parameters {
  vector[n] phi1;
  vector[n] phi2;
  vector[2] nought;
  real<lower = 0> tau1;
  real<lower = 0> tau2;
  real<lower = 0, upper = 1> alpha;
  real eta0;
  real<lower = 0> v_sig1;
  real<lower = 0> v_sig2;
}
transformed parameters{
  vector[n] rho_1;
  vector[n] c_phi2;
  vector[m] rho_2;
  c_phi2 = phi2*inv_sqrt(tau2);
  rho_1 = eta0*c_phi2 + phi1*inv_sqrt(tau1) + nought[1];
  rho_2 = nought[2] + M_W*c_phi2;
  
}
model {
  phi2 ~ sparse_car(alpha, W_sparse, D_sparse, lambda, n, W_n);
  phi1 ~ sparse_car(alpha, W_sparse, D_sparse, lambda, n, W_n);
  nought ~ normal(0, 5);
  eta0 ~ normal(0, 5);
  tau1 ~ normal(0, 5);
  tau2 ~ normal(0, 5);
  v_sig1 ~ gamma(1, .1);
  v_sig2 ~ gamma(1, .1);
  y2 ~ neg_binomial_2(exp(rho_2 + log_offset2), v_sig2);
  y1 ~ neg_binomial_2(exp(rho_1 + log_offset1), v_sig1);
}
generated quantities {
  real corr_risk;
  vector[n] log_lik1;
  vector[m] log_lik2;
  vector[n] log_lik_rep1;
  vector[m] log_lik_rep2;
  int<lower = 0> yrep1[n];
  int<lower = 0> yrep2[m];
  real sum_ll;
  real sum_ll_rep;
  real ppp1;
  real ppp2;
  real ppp;

  // Correlation
  corr_risk = corr_val(exp(rho_1), exp(c_phi2 + nought[2]), n);

  // Simulate from posterior
  for (i in 1:n) {
    // likelihood of the current parameter values (given the original data)
    log_lik1[i] = neg_binomial_2_lpmf(y1[i] | 
    exp(rho_1[i] + log_offset1[i]), v_sig1);
    // generate new data based on current parameter values
    yrep1[i] = neg_binomial_2_rng(exp(rho_1[i] + log_offset1[i]), v_sig1);
    // compute likelihood of the current parameter values
    // (given the new data)
    log_lik_rep1[i] = neg_binomial_2_lpmf(yrep1[i] | 
    exp(rho_1[i] + log_offset1[i]), v_sig1);
  }
  for (i in 1:m) {
    // likelihood of the current parameter values (given the original data)
    log_lik2[i] = neg_binomial_2_lpmf(y2[i] | 
    exp(rho_2[i] + log_offset2[i]), v_sig2);
    // generate new data based on current parameter values
    yrep2[i] = neg_binomial_2_rng(exp(rho_2[i] + log_offset2[i]), v_sig2);
    // compute likelihood of the current parameter values 
    // (given the new data)
    log_lik_rep2[i] = neg_binomial_2_lpmf(yrep2[i] | 
    exp(rho_2[i] + log_offset2[i]), v_sig2);
    }
    // sum up the likelihoods for all observations
    sum_ll = sum(log_lik1) + sum(log_lik2);
    sum_ll_rep = sum(log_lik_rep1) + sum(log_lik_rep2);
    // check which is higher
    ppp1 = sum(log_lik1) > sum(log_lik_rep1) ? 1 : 0;
    ppp2 = sum(log_lik2) > sum(log_lik_rep2) ? 1 : 0;
    ppp = sum_ll > sum_ll_rep ? 1 : 0;
}
