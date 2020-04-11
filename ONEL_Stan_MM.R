#### GENERAL SETUP ####

path <- file.path("....../StanCode")
setwd(path)

# Load package
library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Data load - North East: 4 boroughs coverage 90% (OLD)
load("./ONEL_Data.RData")
#

#### DATASET PREPARATION ####

# Covariates dataset and min-max normalization
X_cent <- apply(cbind(sa_pop = sa_pop_data,
                      imd = imd_data),
                2, function(x) (x-min(x))/diff(range(x)))

# Create datasets for RStan
sp_d_mm <- list(n = length(mortality), # number of MSOAs
                m = length(prevalence), # number of GP
                # Outcomes
                y1 = mortality, # observed number of cases 1
                y2 = prevalence, # observed number of cases 2#
                # Offsets
                log_offset1 = log(exp_mortality),
                log_offset2 = log(exp_prevalence),
                # Adjacecncy
                W_n = sum(W) / 2, # number of neighbor pairs
                W = W, # adjacency matrix
                # Multiple membership 
                M_W = weight_ext)

sp_d_mm_cov <- list(n = length(mortality), # number of MSOAs
                    m = length(prevalence), # number of GP
                    # Outcomes
                    y1 = mortality, # observed number of cases 1
                    y2 = prevalence, # observed number of cases 2#
                    # Covariates
                    k = ncol(X_cent),
                    X = X_cent,
                    # Offsets
                    log_offset1 = log(exp_mortality),
                    log_offset2 = log(exp_prevalence),
                    # Adjacecncy
                    W_n = sum(W) / 2, # number of neighbor pairs
                    W = W, # adjacency matrix
                    # Multiple membership
                    M_W = weight_ext)
#

#### MODELS COMPILATION ####

# GMCAR
## NOCOV
mod_GMCAR_MM_NegBin <- 
  stan_model(file = "Models/GMCAR/GMCAR_MM_NegBin.stan")
## COV
mod_GMCAR_MM_NegBin_cov <- 
  stan_model(file = "Models/GMCAR/GMCAR_MM_NegBin_cov.stan")

# MCAR
## NOCOV
mod_MVCAR_MM_NegBin <- 
  stan_model(file = "Models/MCAR/MCAR_Gelfand_MM.stan")
## COV
mod_MVCAR_MM_NegBin_cov <- 
  stan_model(file = "Models/MCAR/MCAR_Gelfand_MM_Cov.stan")
#

#### HMC PARAMETERS ####

niter <- 6E3
nchains <- 4
#

#### SAMPLING ####

# GMCAR
## NOCOV
fit_gmcar_nocov_negbin <- sampling(mod_GMCAR_MM_NegBin, data = sp_d_mm,
                                   iter = niter, chains = nchains,
                                   control = list(adapt_delta = .99, 
                                                  max_treedepth = 15))
## COV
fit_gmcar_cov_negbin <- sampling(mod_GMCAR_MM_NegBin_cov, 
                                 data = sp_d_mm_cov,
                                 iter = niter, chains = nchains,
                                 control = list(adapt_delta = .99, 
                                                max_treedepth = 15))

# MCAR
## NOCOV
fit_mvcar_nocov_negbin <- sampling(mod_MVCAR_MM_NegBin,
                                   data = sp_d_mm,
                                   iter = niter, chains = nchains,
                                   control = list(adapt_delta = .99, 
                                                  max_treedepth = 15))
## COV
fit_mvcar_cov_negbin <- sampling(mod_MVCAR_MM_NegBin_cov,
                                 data = sp_d_mm_cov,
                                 iter = niter, chains = nchains,
                                 control = list(adapt_delta = .99, 
                                                max_treedepth = 15))
#

#### PRINT OUTPUTS ####

# GMCAR
## NOCOV
pars_vec <- c('nought', 'alpha1', 'alpha2', 'eta0', 'eta1',
              'tau1', 'tau2', 'v_sig1', 'v_sig2', 'ppp1', 'ppp2')
print(fit_gmcar_nocov_negbin, pars = pars_vec, probs = c(.025,.5,.975))
## COV
pars_vec <- c('nought', 'alpha1', 'alpha2', 'eta0', 'eta1', 'beta1', 'beta2',
              'tau1', 'tau2', 'v_sig1', 'v_sig2', 'ppp1', 'ppp2')
print(fit_gmcar_cov_negbin, pars = pars_vec, probs = c(.025,.5,.975))

# MVCAR
## NOCOV
pars_vec <- c('nought', 'alpha', 'eta0', 
              'tau1', 'tau2', 'v_sig1', 'v_sig2', 'ppp1', 'ppp2')
print(fit_mvcar_nocov_negbin, pars = pars_vec, probs = c(.025,.5,.975))
## COV
pars_vec <- c('nought', 'alpha', 'eta0', 'beta1', 'beta2',
              'tau1', 'tau2', 'v_sig1', 'v_sig2', 'ppp1', 'ppp2')
print(fit_mvcar_cov_negbin, pars = pars_vec, probs = c(.025,.5,.975))
#

