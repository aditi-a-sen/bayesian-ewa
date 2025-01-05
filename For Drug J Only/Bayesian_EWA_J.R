##############
#Packages
##############

library(rjags)
library(coda)
## install.packages('MCMCvis')
library(MCMCvis)
# Changes to implement Parallel Processing (10/20)

#################################################################################################################
#Statistical Model
#################################################################################################################
model_string_J<-"

model{
 
     #########################
     #Likelihood Specification
     #########################
     for(i in 1:n1){
        for(j in 1:n2){
           y[1,i,j] ~ dcat(p[1,i,j,1:m[i,j]])
           for(k in 1:max(m)){
              p_temp[1,i,j,k]<-exp(lambda[i,j]*V0[i,j,k])
              V[1,i,j,k]<-((phi[i,j])*V0[i,j,k]*N0 + f[1,i,j,k])/N[1,i,j]
              f[1,i,j,k]<-(delta[i,j] + (1 - delta[i,j])*equals(y[1,i,j], k))*pi[i,j,1,k]
              }
           for(k in 1:max(m)){
              p[1,i,j,k]<-p_temp[1,i,j,k]/sum(p_temp[1,i,j,1:m[i,j]])
              }        
           N[1,i,j]<-(rho[i,j])*N0 + 1
           #model_fit[1,i,j]<-p[1,i,j,y[1,i,j]]
        
           for(t in 2:n3){
              y[t,i,j] ~ dcat(p[t,i,j,1:m[i,j]])              
              for(k in 1:max(m)){
                 p_temp[t,i,j,k]<-exp(lambda[i,j]*V[(t-1),i,j,k])
                 V[t,i,j,k]<-((phi[i,j]^t)*V0[i,j,k]*N0 + inprod(f[1:t,i,j,k], (phi[i,j]^(t - c(1:t)))))/N[t,i,j]
                 f[t,i,j,k]<-(delta[i,j] + (1 - delta[i,j])*equals(y[t,i,j], k))*pi[i,j,t,k]
                 }
              for(k in 1:max(m)){
                 p[t,i,j,k]<-p_temp[t,i,j,k]/sum(p_temp[t,i,j,1:m[i,j]])
                 }
              N[t,i,j] <- (rho[i,j]^t)*N0 + (1 - (rho[i,j]^t))/(1 - rho[i,j])
              #model_fit[t,i,j]<-p[t,i,j,y[t,i,j]]
              }
           }
        }
       
     ####################
     #Prior Distributions
     ####################
     for(i in 1:n1){
        for(j in 1:n2){
           phi[i,j]<-1/(1 + exp(-logit_phi[i,j]))
           logit_phi[i,j] ~ dnorm(mean_phi[i,j], sigma2_phi_inv)
           mean_phi[i,j]<-x[i,]%*%beta_phi + z[i,j,]%*%gamma_phi + RE[i,1]  

           delta[i,j]<-1/(1 + exp(-logit_delta[i,j]))
           logit_delta[i,j] ~ dnorm(mean_delta[i,j], sigma2_delta_inv)
           mean_delta[i,j]<-x[i,]%*%beta_delta + z[i,j,]%*%gamma_delta + RE[i,2]  

           rho[i,j]<-1/(1 + exp(-logit_rho[i,j]))
           logit_rho[i,j] ~ dnorm(mean_rho[i,j], sigma2_rho_inv)
           mean_rho[i,j]<-x[i,]%*%beta_rho + z[i,j,]%*%gamma_rho + RE[i,3]  

           lambda[i,j]<-exp(log_lambda[i,j])
           log_lambda[i,j] ~ dnorm(mean_lambda[i,j], sigma2_lambda_inv)
           mean_lambda[i,j]<-x[i,]%*%beta_lambda + z[i,j,]%*%gamma_lambda + RE[i,4] 
           } 

        RE[i,1:4] ~ dmnorm(RE_mean, RE_Sigma_inv)
        }

     RE_mean[1]<-0
     RE_mean[2]<-0
     RE_mean[3]<-0
     RE_mean[4]<-0

     RE_Sigma_inv[1:4,1:4] ~ dwish(Omega[1:4,1:4], df)
     df<-4 + 1
     RE_Sigma[1:4,1:4]<-inverse(RE_Sigma_inv[1:4,1:4])

     for(j in 1:p_x){
        beta_phi[j] ~ dnorm(0, 0.0001)
        beta_delta[j] ~ dnorm(0, 0.0001)
        beta_rho[j] ~ dnorm(0, 0.0001)
        beta_lambda[j] ~ dnorm(0, 0.0001)
        }
   
     for(j in 1:p_z){
        gamma_phi[j] ~ dnorm(0, 0.0001)
        gamma_delta[j] ~ dnorm(0, 0.0001)
        gamma_rho[j] ~ dnorm(0, 0.0001)
        gamma_lambda[j] ~ dnorm(0, 0.0001)
        }

     sigma2_phi_inv<-1/(sigma_phi*sigma_phi)
     sigma_phi ~ dunif(0.00, 1000.00)
     sigma2_delta_inv<-1/(sigma_delta*sigma_delta)
     sigma_delta ~ dunif(0.00, 1000.00)
     sigma2_rho_inv<-1/(sigma_rho*sigma_rho)
     sigma_rho ~ dunif(0.00, 1000.00)
     sigma2_lambda_inv<-1/(sigma_lambda*sigma_lambda)
     sigma_lambda ~ dunif(0.00, 1000.00)
     }
" 

####################################################
#Model Organization
####################################################

model_jags_J<-jags.model(data=list('y' = y_J, 
                                   'pi' = pi_array,
                                   'x' = x_J, 
                                   'p_x' = p_x_J,
                                   'z' = z_J,
                                   'p_z' = p_z_J,
                                   'n1' = n1_J,
                                   'n2' = n2_J,
                                   'n3' = n3_J,
                                   'm' = m_J, 
                                   'V0' = V0_J,
                                   'N0' = N0,
                                   'Omega' = diag(4)),
                         textConnection(model_string_J),
                         n.chains=1)

#model_jags_2<-jags.model(data=list('y' = y_2, 
# 'pi' = pi_2,
# 'x' = x_2, 
#'p_x' = p_x_2,
# 'z' = z_2,
# 'p_z' = p_z_2,
# 'n1' = n1_2,
# 'n2' = n2_2,
# 'n3' = n3_2,
# 'm' = m_2, 
# 'V0' = V0_2,
# 'N0' = N0,
# 'Omega' = diag(4)),
# textConnection(model_string),
# n.chains=1)

################################################################
#Posterior Sampling
################################################################
#update(model_jags_2, 
#       n.iter=100000)  #Burnin

update(model_jags_J, 
       n.iter=500) #Reduced Burnin

posterior_samples_J<-coda.samples(model_jags_J, 
                                  variable.names=c("phi",
                                                   "beta_phi",
                                                   "gamma_phi",
                                                   "sigma_phi",
                                                   "delta",
                                                   "beta_delta",
                                                   "gamma_delta",
                                                   "sigma_delta",
                                                   "rho",
                                                   "beta_rho",
                                                   "gamma_rho",
                                                   "sigma_rho",
                                                   "lambda",
                                                   "beta_lambda",
                                                   "gamma_lambda",
                                                   "sigma_lambda",
                                                   "RE",
                                                   "RE_Sigma"),
                                  thin=10,
                                  n.iter=10000)

#posterior_samples_1_delta <- coda.samples(model_jags_1, 
#variable.names=c("delta"),
#thin=50,
#n.iter=1000)

##########################
#Posterior Inference
##########################

par(ask=FALSE)


par(mar=c(1,1,1,1))


plot(posterior_samples_J)



summary(posterior_samples_J)

sum = summary(posterior_samples_J)

################################
#MCMC Diagnostics
################################
effectiveSize(posterior_samples_J)
geweke.diag(posterior_samples_J)

MCMCplot(posterior_samples_J, params=c("delta"))
deltas_J <- MCMCsummary(posterior_samples_J, params = 'delta')
hist(deltas_J$mean)

MCMCplot(posterior_samples_J, params=c("lambda"))
lambdas_J <- MCMCsummary(posterior_samples_J, params = 'lambda')
hist(lambdas_J$mean)

MCMCplot(posterior_samples_J, params=c("phi"))
phis_J <- MCMCsummary(posterior_samples_J, params = 'phi')
hist(phis_J$mean)

MCMCplot(posterior_samples_J, params=c("rho"))
rhos_J <- MCMCsummary(posterior_samples_J, params = 'rho')
hist(rhos_J$mean)

MCMCsummary(posterior_samples_J, params = 'sigma_delta')
MCMCsummary(posterior_samples_J, params = 'sigma_lambda')
MCMCsummary(posterior_samples_J, params = 'sigma_phi')
MCMCsummary(posterior_samples_J, params = 'sigma_rho')


MCMCsummary(posterior_samples_J, params = c('beta_delta'), ISB = T, exact = F, round=2)
MCMCsummary(posterior_samples_J, params = c('beta_lambda'), ISB = T, exact = F, round=2)
MCMCsummary(posterior_samples_J, params = c('beta_phi'), ISB = T, exact = F, round=2)
MCMCsummary(posterior_samples_J, params = c('beta_rho'), ISB = T, exact = F, round=2)

mean(MCMCsummary(posterior_samples_J, params = 'delta')$mean,na.rm=T)
sd(MCMCsummary(posterior_samples_J, params = 'delta')$mean,na.rm=T)

mean(MCMCsummary(posterior_samples_J, params = 'lambda')$mean,na.rm=T)
sd(MCMCsummary(posterior_samples_J, params = 'lambda')$mean,na.rm=T)

mean(MCMCsummary(posterior_samples_J, params = 'phi')$mean,na.rm=T)
sd(MCMCsummary(posterior_samples_J, params = 'phi')$mean,na.rm=T)

mean(MCMCsummary(posterior_samples_J, params = 'rho')$mean,na.rm=T)
sd(MCMCsummary(posterior_samples_J, params = 'rho')$mean,na.rm=T)

write.csv(dvr_demo_df)
print("done")
