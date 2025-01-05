library(rjags)
library(coda)
## install.packages('MCMCvis')
library(MCMCvis)
###--- Model Creation

library(rjags)
library(coda)
library(MCMCvis)
library(parallel)
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)


model_string2 <- "
model {

    #########################
    #Likelihood Specification
    #########################
    for (i in 1:n1) {
        for (j in 1:n4) {
            N[1,i,j] <- (rho[i,j]) * N0 + 1
            
            # Vectorized calculation for N[t,i,j]
            for (t in 2:n3) {
                N[t,i,j] <- (rho[i,j]^t) * N0 + (1 - (rho[i,j]^t)) / (1 - rho[i,j]);
            }
            
            for (k in 1:n2) {
                y[1,i,j,k] ~ dcat(p[1,i,j,k,1:m[i,j,k]]);
                
                # Vectorized calculation of p_temp and f
                for (g in 1:max(m)) {
                    p_temp[1,i,j,k,g] <- exp(lambda[i,j] * V0[i,j,k,g]);
                    V[1,i,j,k,g] <- ((phi[i,j]) * V0[i,j,k,g] * N0 + f[1,i,j,k,g]) / N[1,i,j];
                    f[1,i,j,k,g] <- (delta[i,j] + (1 - delta[i,j]) * equals(y[1,i,j,k], g)) * pi[i,j,k,1,g];
                }

                # Normalize p_temp to get p
                p[1,i,j,k,1:max(m)] <- p_temp[1,i,j,k,1:max(m)] / sum(p_temp[1,i,j,k,1:m[i,j,k]]);

                # Vectorized calculation for t > 1
                for (t in 2:n3) {
                    y[t,i,j,k] ~ dcat(p[t,i,j,k,1:m[i,j,k]]);
                    
                    for (g in 1:max(m)) {
                        p_temp[t,i,j,k,g] <- exp(lambda[i,j] * V[(t-1),i,j,k,g]);
                        
                        # Vectorized calculation of V and f
                        V[t,i,j,k,g] <- ((phi[i,j]^t) * V0[i,j,k,g] * N0 + 
                                         inprod(f[1:t,i,j,k,g], phi[i,j]^(t - c(1:t)))) / N[t,i,j];
                        
                        f[t,i,j,k,g] <- (delta[i,j] + (1 - delta[i,j]) * equals(y[t,i,j,k], g)) * pi[i,j,k,t,g];
                    }
                    
                    # Normalize p_temp to get p
                    p[t,i,j,k,1:max(m)] <- p_temp[t,i,j,k,1:max(m)] / sum(p_temp[t,i,j,k,1:m[i,j,k]]);
                }
            }
        }
    }

    ####################
    #Prior Distributions
    ####################
    for (i in 1:n1) {
        for (j in 1:n4) {
            phi[i,j] <- 1 / (1 + exp(-logit_phi[i,j]));
            logit_phi[i,j] ~ dnorm(mean_phi[i,j], sigma2_phi_inv);
            mean_phi[i,j] <- x[i,] %*% beta_phi + z[i,j,] %*% gamma_phi + RE[i, 1];

            delta[i,j] <- 1 / (1 + exp(-logit_delta[i,j]));
            logit_delta[i,j] ~ dnorm(mean_delta[i,j], sigma2_delta_inv);
            mean_delta[i,j] <- x[i,] %*% beta_delta + z[i,j,] %*% gamma_delta + RE[i, 2];

            rho[i,j] <- 1 / (1 + exp(-logit_rho[i,j]));
            logit_rho[i,j] ~ dnorm(mean_rho[i,j], sigma2_rho_inv);
            mean_rho[i,j] <- x[i,] %*% beta_rho + z[i,j,] %*% gamma_rho + RE[i, 3];

            lambda[i,j] <- exp(log_lambda[i,j]);
            log_lambda[i,j] ~ dnorm(mean_lambda[i,j], sigma2_lambda_inv);
            mean_lambda[i,j] <- x[i,] %*% beta_lambda + z[i,j,] %*% gamma_lambda + RE[i, 4];
        }
        RE[i, 1:4] ~ dmnorm(RE_mean, RE_Sigma_inv);
    }

    RE_mean[1:4] <- rep(0, 4);

    RE_Sigma_inv[1:4, 1:4] ~ dwish(Omega[1:4, 1:4], df);
    df <- 4 + 1;
    RE_Sigma[1:4, 1:4] <- inverse(RE_Sigma_inv[1:4, 1:4]);

    for (j in 1:p_x) {
        beta_phi[j] ~ dnorm(0, 0.0001);
        beta_delta[j] ~ dnorm(0, 0.0001);
        beta_rho[j] ~ dnorm(0, 0.0001);
        beta_lambda[j] ~ dnorm(0, 0.0001);
    }

    for (j in 1:p_z) {
        gamma_phi[j] ~ dnorm(0, 0.0001);
        gamma_delta[j] ~ dnorm(0, 0.0001);
        gamma_rho[j] ~ dnorm(0, 0.0001);
        gamma_lambda[j] ~ dnorm(0, 0.0001);
    }

    sigma2_phi_inv <- 1 / (sigma_phi * sigma_phi);
    sigma_phi ~ dunif(0.00, 1000.00);
    sigma2_delta_inv <- 1 / (sigma_delta * sigma_delta);
    sigma_delta ~ dunif(0.00, 1000.00);
    sigma2_rho_inv <- 1 / (sigma_rho * sigma_rho);
    sigma_rho ~ dunif(0.00, 1000.00);
    sigma2_lambda_inv <- 1 / (sigma_lambda * sigma_lambda);
    sigma_lambda ~ dunif(0.00, 1000.00);
}
"

model_jags_pet_m7 <- jags.model(data = list('y' = y_m7,
                                            'pi' = pi_m7,
                                            'x' = x_m7,
                                            'p_x' = p_x_m7,
                                            'z' = z_m7,
                                            'p_z' = p_z_m7,
                                            'n1' = n1_m7,
                                            'n2' = n2_m7,
                                            'n3' = n3_m7,
                                            'n4' = n4_m7,
                                            'm' = m_m7,
                                            'V0' = V0_m7,
                                            'N0' = N0,
                                            'Omega' = diag(4)),
                                textConnection(model_string2),
                                n.chains = n_cores)
update(model_jags_pet_m7, n.iter=500)

posterior_samples_pet_m7 <- coda.samples(model_jags_pet_m7,
                                         variable.names=c("phi", "beta_phi", "gamma_phi", "sigma_phi",
                                                          "delta", "beta_delta", "gamma_delta", "sigma_delta",
                                                          "rho", "beta_rho", "gamma_rho", "sigma_rho",
                                                          "lambda", "beta_lambda", "gamma_lambda", "sigma_lambda",
                                                          "RE", "RE_Sigma"),
                                         thin=10,
                                         n.iter=10000,
                                         n.chains=n_cores)
update(model_jags_pet_m7,
       n.iter=500)
###--- Outputs

geweke.diag(posterior_samples_pet_m7)


mean(MCMCsummary(posterior_samples_pet_m7, params = 'delta')$mean,na.rm=T)
sd(MCMCsummary(posterior_samples_pet_m7, params = 'delta')$mean,na.rm=T)

mean(MCMCsummary(posterior_samples_pet_m7, params = 'lambda')$mean,na.rm=T)
sd(MCMCsummary(posterior_samples_pet_m7, params = 'lambda')$mean,na.rm=T)

mean(MCMCsummary(posterior_samples_pet_m7, params = 'phi')$mean,na.rm=T)
sd(MCMCsummary(posterior_samples_pet_m7, params = 'phi')$mean,na.rm=T)

mean(MCMCsummary(posterior_samples_pet_m7, params = 'rho')$mean,na.rm=T)
sd(MCMCsummary(posterior_samples_pet_m7, params = 'rho')$mean,na.rm=T)

MCMCsummary(posterior_samples_pet_m7, params = c('beta_delta'), ISB = T, exact = F, round=2)

MCMCsummary(posterior_samples_pet_m7, params = c('beta_lambda'), ISB = T, exact = F, round=2)

MCMCsummary(posterior_samples_pet_m7, params = c('beta_phi'), ISB = T, exact = F, round=2)

MCMCsummary(posterior_samples_pet_m7, params = c('beta_rho'), ISB = T, exact = F, round=2)

MCMCsummary(posterior_samples_pet_m7, params = 'sigma_delta')

MCMCsummary(posterior_samples_pet_m7, params = 'sigma_lambda')

MCMCsummary(posterior_samples_pet_m7, params = 'sigma_phi')

MCMCsummary(posterior_samples_pet_m7, params = 'sigma_rho')

deltas <- MCMCsummary(posterior_samples_pet_m7, params = 'delta')
hist(deltas$mean)

lambdas <- MCMCsummary(posterior_samples_pet_m7, params = 'lambda')
hist(lambdas$mean)

phis <- MCMCsummary(posterior_samples_pet_m7, params = 'phi')
hist(phis$mean)

rhos <- MCMCsummary(posterior_samples_pet_m7, params = 'rho')
hist(rhos$mean)