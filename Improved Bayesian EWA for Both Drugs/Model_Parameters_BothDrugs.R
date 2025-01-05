library(dplyr)
library(stats)
library(tidyr)
library(readr)
library(rjags)
library(coda)
library(MCMCvis)
library(mcmcplots)
b_Analysis_Data = read_csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 2/patent_race_with_drug.csv")
dvr_demo_df = read_csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 2/dvr_demo_df.csv")
dopamine_master_demographics_dec_23 = read_csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 2/dopamine_master_demographics_dec_23.csv")
b_Analysis_Data <- subset(b_Analysis_Data, SID %in% dvr_demo_df$UCB_ID)
b_ids_to_remove <- c(6022, 6100, 6103, 6144, 6171)
b_Analysis_Data <- b_Analysis_Data %>% filter(!(SID %in% b_ids_to_remove))
## Check if Strong Role
b_Analysis_Data$Strong_Role_Indicator = ifelse(b_Analysis_Data$Condition=='s',1,ifelse(b_Analysis_Data$Condition=='w',0,NA))
# Number of unique subjects
b_n1_m7 = length(unique(b_Analysis_Data$SID))
# Number of unique runs
b_n2_m7 = length(unique(b_Analysis_Data$Run))
# Number of unique rounds
b_n3_m7 = length(unique(b_Analysis_Data$Round))
# Number of unique sessions
b_n4_m7 = length(unique(b_Analysis_Data$Session))
# Create arrays from above variables
b_m_m7<-array(NA, dim = c(b_n1_m7,b_n4_m7,b_n2_m7))
b_w_m7<-array(NA, dim = c(b_n1_m7,b_n4_m7,b_n2_m7))
b_r_m7<-array(NA, dim = c(b_n1_m7,b_n4_m7,b_n2_m7))
b_drug_m7<-array(NA, dim = c(b_n1_m7,b_n4_m7,b_n2_m7))
b_w_m7_j<-array(NA, dim = c(b_n1_m7,b_n4_m7))
b_r_m7_j<-array(NA, dim = c(b_n1_m7,b_n4_m7))
b_drug_m7_j<-array(NA, dim = c(b_n1_m7,b_n4_m7))
#form 4 arrays (max investment = 6, strong role, session, and drug)
for(i in 1:b_n1_m7){
  for(j in 1:b_n4_m7){
    for (k in 1:b_n2_m7){
      b_m_m7[i,j,k]<- 5 + mean(b_Analysis_Data$Strong_Role_Indicator[b_Analysis_Data$SID==sort(unique(b_Analysis_Data$SID))[i] & b_Analysis_Data$Session==sort(unique(b_Analysis_Data$Session))[j] & b_Analysis_Data$Run==sort(unique(b_Analysis_Data$Run))[k]],na.rm = T)
      b_r_m7[i,j,k]<-mean(b_Analysis_Data$Strong_Role_Indicator[b_Analysis_Data$SID==sort(unique(b_Analysis_Data$SID))[i] & b_Analysis_Data$Session==sort(unique(b_Analysis_Data$Session))[j] & b_Analysis_Data$Run==sort(unique(b_Analysis_Data$Run))[k]],na.rm = T)
      b_w_m7[i,j,k]<-mean(b_Analysis_Data$Session[b_Analysis_Data$SID==sort(unique(b_Analysis_Data$SID))[i] & b_Analysis_Data$Session==sort(unique(b_Analysis_Data$Session))[j] & b_Analysis_Data$Run==sort(unique(b_Analysis_Data$Run))[k]],na.rm=T)
      b_drug_m7[i,j,k]<-b_Analysis_Data[match(unique(b_Analysis_Data$Drug[b_Analysis_Data$SID==sort(unique(b_Analysis_Data$SID))[i] & b_Analysis_Data$Session==sort(unique(b_Analysis_Data$Session))[j] & b_Analysis_Data$Run==sort(unique(b_Analysis_Data$Run))[k]]), b_Analysis_Data$Drug),]$Drug
    }
  }
}

for(i in 1:b_n1_m7){
  for(j in 1:b_n4_m7){
    b_w_m7_j[i,j]<-mean(b_Analysis_Data$Session[b_Analysis_Data$SID==sort(unique(b_Analysis_Data$SID))[i] & b_Analysis_Data$Session==sort(unique(b_Analysis_Data$Session))[j]],na.rm=T)
    b_drug_m7_j[i,j]<-b_Analysis_Data[match(unique(b_Analysis_Data$Drug[b_Analysis_Data$SID==sort(unique(b_Analysis_Data$SID))[i] & b_Analysis_Data$Session==sort(unique(b_Analysis_Data$Session))[j]]), b_Analysis_Data$Drug),]$Drug
    b_r_m7_j[i,j]<-mean(b_Analysis_Data$Strong_Role_Indicator[b_Analysis_Data$SID==sort(unique(b_Analysis_Data$SID))[i] & b_Analysis_Data$Session==sort(unique(b_Analysis_Data$Session))[j]],na.rm = T)
  }
}

b_drug_m7_num <- ifelse(b_drug_m7_j == "H", 1, 0)
b_V0_m7<-array(NA, c(b_n1_m7, b_n4_m7, b_n2_m7, max(b_m_m7,na.rm=T)))
for(i in 1:b_n1_m7){
  for(j in 1:b_n2_m7){
    for(k in 1:b_n4_m7){
      #Weak Role, Run 1, Session 1
      if((j==1) & (b_m_m7[i,k,j]==5) & (k==1)){
        b_V0_m7[i,k,j,]<-(c(1.82, 0, 0.956, 0.956, 1.89, 0) + 5)
      } #add 0 to end when including weak+strong together
      #Strong Role, Run 1, Session 1
      if(j==1 & b_m_m7[i,k,j]==6 & (k==1)){
        b_V0_m7[i,k,j,]<-(c(0.405, 0, 0.511, 1.15, 1.34, 1.70) + 5)
      }
      #Weak Role, Run 2, Session 1
      if((j==2) & (b_m_m7[i,k,j]==5)& (k==1)){
        b_V0_m7[i,k,j,]<-(c(2.73, 1.20, 0, 1.10, 2.37, 0) + 5)
      } #add 0 to end when including weak+strong together
      #Strong Role, Run 2, Session 1
      if(j==2 & b_m_m7[i,k,j]==6  & (k==1)){
        b_V0_m7[i,k,j,]<-(c(0, 0.560, 0.223, 1.10, 1.70, 2.42) + 5)
      }
      #Weak Role, Run 1, Session 2
      if((j==1) & (b_m_m7[i,k,j]==5) & (k==2)){
        b_V0_m7[i,k,j,]<-(c(2.56, .336, 0, 0.182, 0.875, 0) + 5)
      } #add 0 to end when including weak+strong together
      #Strong Role, Run 1, Session 2
      if(j==1 & b_m_m7[i,k,j]==6 & (k==2)){
        b_V0_m7[i,k,j,]<-(c(0, 2.71, 2.08, 2.89, 2.48, 3.83) + 5)
      }
      #Weak Role, Run 2, Session 2
      if((j==2) & (b_m_m7[i,k,j]==5)  & (k==2)){
        b_V0_m7[i,k,j,]<-(c(2.71, 0, 0.811, 0.560, 1.61, 0) + 5)
      } #add 0 to end when including weak+strong together
      #Strong Role, Run 2, Session 2
      if(j==2 & b_m_m7[i,k,j]==6  & (k==2)){
        b_V0_m7[i,k,j,]<-(c(0.154, 0.405, 0, 0.693, 1.43, 1.79) + 5)
      }
    }
  }
}

b_y_m7<-array(NA, c(b_n3_m7, b_n1_m7, b_n4_m7, b_n2_m7))
b_y_star_m7<-array(NA, c(b_n3_m7, b_n1_m7, b_n4_m7, b_n2_m7))
for(t in 1:b_n3_m7){
  for(i in 1:b_n1_m7){
    for(j in 1:b_n2_m7){
      for(k in 1:b_n4_m7){
        b_y_m7[t,i,k,j]<-b_Analysis_Data$MyInvestment[b_Analysis_Data$Round==sort(unique(b_Analysis_Data$Round))[t] & b_Analysis_Data$SID==sort(unique(b_Analysis_Data$SID))[i] & b_Analysis_Data$Session==sort(unique(b_Analysis_Data$Session))[k] & b_Analysis_Data$Run==sort(unique(b_Analysis_Data$Run))[j]]
        b_y_star_m7[t,i,k,j]<-b_Analysis_Data$OpponentInvest[b_Analysis_Data$Round==sort(unique(b_Analysis_Data$Round))[t] & b_Analysis_Data$SID==sort(unique(b_Analysis_Data$SID))[i] & b_Analysis_Data$Session==sort(unique(b_Analysis_Data$Session))[k] & b_Analysis_Data$Run==sort(unique(b_Analysis_Data$Run))[j]]
      }
    }
  }
}

b_y_m7<-b_y_m7 + 1
b_y_star_m7<-b_y_star_m7 + 1
q<-10
N0<-1
b_pi_m7<-array(NA, c(b_n1_m7, b_n4_m7, b_n2_m7, b_n3_m7, max(b_m_m7)))
for(i in 1:b_n1_m7){
  for(j in 1:b_n4_m7){
    for(k in 1:b_n2_m7){
      for(t in 1:b_n3_m7){
        for(g in 1:max(b_m_m7)){
          b_pi_m7[i,j,k,t,g]<-((b_m_m7[i,j,k] - 1) + q*as.numeric((g - 1) >= b_y_star_m7[t,i,j,k]) - (g - 1))
        }
      }
    }
  }
}

b_Analysis_Data <- b_Analysis_Data %>% arrange(SID)
dopamine_master_demographics_dec_23 <- dopamine_master_demographics_dec_23 %>% arrange(Participant.ID.)
b_x_m7<-cbind(dopamine_master_demographics_dec_23$Age, as.numeric(dopamine_master_demographics_dec_23$Sex_assigned_at_birth=="Male"))
b_x_m7[,1]<-scale(b_x_m7[,1])
b_x_m7<-cbind(rep(1, times=nrow(b_x_m7)), b_x_m7)
b_p_x_m7<-ncol(b_x_m7)

# z and p_z
b_z_m7<-array(0, dim=c(b_n1_m7, b_n4_m7, 2))
for(i in 1:b_n1_m7){
  for(j in 1:b_n4_m7){
    b_z_m7[i,j,]<-c((b_w_m7_j[i,j] - 1), b_drug_m7_num[i,j])
  }
}
b_p_z_m7<-2