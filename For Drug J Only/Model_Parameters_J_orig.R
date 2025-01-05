#---- Packages
library(dplyr)
library(stats)
library(tidyr)

# Load data frames
patent_race_J_only= read.csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 1/patent_race_J_only.csv")
dvr_demo_df = read.csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 1/dvr_demo_df.csv")
dopamine_master_demographics_dec_23 = read.csv("~/Downloads/bayesian-ewa/attempt 1/dopamine_master_demographics_dec_23.csv")

#---- Data Cleaning and Filtering

# Remove missing values
Analysis_Data_J = na.omit(patent_race_J_only[,1:12])

# Subset by responded ID's
Analysis_Data_J <- subset(Analysis_Data_J, UCB_ID %in% dvr_demo_df$UCB_ID)

# Remove bad ID's
ids_to_remove <- c(6022, 6100, 6103, 6144, 6171)

Analysis_Data_J <- Analysis_Data_J %>%
  filter(!(UCB_ID %in% ids_to_remove))

Analysis_Data_J <- Analysis_Data_J %>%
  filter(!(UCB_ID==6156))

#---- Summary Statistics

# Number of unique subjects
n1_J = length(unique(Analysis_Data_J$UCB_ID))

# Number of unique runs
n2_J = length(unique(Analysis_Data_J$Run))

# Number of unique rounds
n3_J = length(unique(Analysis_Data_J$Round))

#---- Variables and Calculations

# Categorize conditions as 1 or 0 (strong or weak)
Analysis_Data_J$Strong_Role_Indicator = ifelse(Analysis_Data_J$Condition=='s',1,
                                               ifelse(Analysis_Data_J$Condition=='w',0,NA))
# Matrices to hold calculated values for different subjects and sessions
m_J<-matrix(NA, nrow=n1_J, ncol=n2_J)
w_J<-matrix(NA, nrow=n1_J, ncol=n2_J)
r_J<-matrix(NA, nrow=n1_J, ncol=n2_J)

# Add to matrices based on Strong_Role_Indicator
for(i in 1:n1_J){
  for(j in 1:n2_J){
    m_J[i,j]<- 5 + mean(Analysis_Data_J$Strong_Role_Indicator[Analysis_Data_J$UCB_ID==sort(unique(Analysis_Data_J$UCB_ID))[i] &
                                                                Analysis_Data_J$Run==sort(unique(Analysis_Data_J$Run))[j]])
    
    w_J[i,j]<-mean(Analysis_Data_J$Run[Analysis_Data_J$UCB_ID==sort(unique(Analysis_Data_J$UCB_ID))[i] &
                                         Analysis_Data_J$Run==sort(unique(Analysis_Data_J$Run))[j]])
    
    r_J[i,j]<-mean(Analysis_Data_J$Strong_Role_Indicator[Analysis_Data_J$UCB_ID==sort(unique(Analysis_Data_J$UCB_ID))[i] &
                                                           Analysis_Data_J$Run==sort(unique(Analysis_Data_J$Run))[j]])
  }
}

#---- V0

# Store initial values based on role and session
V0_J<-array(NA, c(n1_J, n2_J, max(m_J,na.rm=T)))

for(i in 1:n1_J){
  for(j in 1:n2_J){
    
    #Weak Role, Session 1
    if((j==1) & (m_J[i,j]==5)){
      V0_J[i,j,]<-c(8.230893, 5.303274, 5, 6.666961, 8.828354, 0)
    }
    
    #Strong Role, Session 1
    if(j==1 & m_J[i,j]==6){
      V0_J[i,j,]<-c(6.363687, 4.89804, 5, 10.57402, 11.624766, 12.53241)
    }
    
    #Weak Role, Session 2
    if((j==2) & (m_J[i,j]==5)){
      V0_J[i,j,]<-c(8.58959, 7.161393, 5, 7.161393, 9.041251, 0)
    }
    
    #Strong Role, Session 2
    if(j==2 & m_J[i,j]==6){
      V0_J[i,j,]<-c(7.464666, 7.161393, 5, 7.959099, 9.717582, 11.081269)
    }
  }
}

#---- y (investment array)

# Arrays to hold MyInvestment and OpponentInvest
y_J<-array(NA, c(n3_J, n1_J, n2_J))
y_star_J<-array(NA, c(n3_J, n1_J, n2_J))
for(t in 1:n3_J){
  for(i in 1:n1_J){
    for(j in 1:n2_J){
      y_J[t,i,j]<-Analysis_Data_J$MyInvestment[Analysis_Data_J$Round==sort(unique(Analysis_Data_J$Round))[t] &
                                                 Analysis_Data_J$UCB_ID==sort(unique(Analysis_Data_J$UCB_ID))[i] &
                                                 Analysis_Data_J$Run==sort(unique(Analysis_Data_J$Run))[j]]
      
      y_star_J[t,i,j]<-Analysis_Data_J$OpponentInvest[Analysis_Data_J$Round==sort(unique(Analysis_Data_J$Round))[t] &
                                                        Analysis_Data_J$UCB_ID==sort(unique(Analysis_Data_J$UCB_ID))[i] &
                                                        Analysis_Data_J$Run==sort(unique(Analysis_Data_J$Run))[j]]
    }
  }
}


y_J<-y_J + 1
y_star_J<-y_star_J + 1

#---- pi

q<-10
N0<-1

# array to hold calculated probabilities based on m_J
pi_J<-array(NA, c(n1_J, n2_J, n3_J, 6))

for(i in 1:n1_J){
  for(j in 1:n2_J){
    for(t in 1:n3_J){
      for(k in 1:6){
        pi_J[i,j,t,k]<-((m_J[i,j] - 1) + q*as.numeric((k - 1) >= y_star_J[t,i,j]) - (k - 1))
      }
    }
  }
}
#---- x and p_x - covariate arrays, based on demographic information and the calculated role metrics

Analysis_Data_J <- Analysis_Data_J %>% arrange(UCB_ID)
dopamine_master_demographics_dec_23 <- dopamine_master_demographics_dec_23 %>% arrange(Participant.ID.)

x_J<-cbind(dopamine_master_demographics_dec_23$Age, as.numeric(dopamine_master_demographics_dec_23$Sex_assigned_at_birth=="Male"))
x_J[,1]<-scale(x_J[,1])

x_J<-cbind(rep(1, times=nrow(x_J)), x_J)
p_x_J<-ncol(x_J)

#---- z and p_z

# covariate array to hold information for each subject-session combination
z_J<-array(0, dim=c(n1_J, n2_J, 2))
for(i in 1:n1_J){
  for(j in 1:n2_J){
    z_J[i,j,]<-c((w_J[i,j] - 1), r_J[i,j])
  }
}

p_z_J<-2