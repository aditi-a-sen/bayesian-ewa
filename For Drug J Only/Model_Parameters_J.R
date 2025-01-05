# Packages
library(dplyr)
library(stats)
library(tidyr)

# Load data frames
patent_race_J_only= read.csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 1/patent_race_J_only.csv")
dvr_demo_df = read.csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 1/dvr_demo_df.csv")
dopamine_master_demographics_dec_23 = read.csv("~/Downloads/bayesian-ewa/attempt 1/dopamine_master_demographics_dec_23.csv")

# Data Cleaning and Filtering
Analysis_Data_J = patent_race_J_only[,1:12] %>%
  na.omit() %>%
  filter(UCB_ID %in% dvr_demo_df$UCB_ID) %>%
  filter(!UCB_ID %in% c(6022, 6100, 6103, 6144, 6171, 6156))

# Summary Statistics
n1_J = n_distinct(Analysis_Data_J$UCB_ID) 
n2_J = n_distinct(Analysis_Data_J$Run)
n3_J = n_distinct(Analysis_Data_J$Round)


# Categorize conditions as strong or weak
Analysis_Data_J <- Analysis_Data_J %>%
  mutate(Strong_Role_Indicator = ifelse(Condition == 's', 1, ifelse(Condition == 'w', 0, NA)))

# Calculate mean values using dplyr
mean_values <- Analysis_Data_J %>%
  group_by(UCB_ID, Run) %>%
  summarize(
    m_value = 5 + mean(Strong_Role_Indicator, na.rm = TRUE),
    w_value = mean(Run, na.rm = TRUE),
    r_value = mean(Strong_Role_Indicator, na.rm = TRUE),
    a_value = 5+Strong_Role_Indicator,
    b_value = Run,
    c_value = Strong_Role_Indicator
  ) %>%
  ungroup()

m_J <- matrix(mean_values$m_value, nrow = n1_J, ncol = n2_J)
w_J <- matrix(mean_values$w_value, nrow = n1_J, ncol = n2_J)
r_J <- matrix(mean_values$r_value, nrow = n1_J, ncol = n2_J)


# V0 Calculation
V0_values <- mean_values %>%
  mutate(V0_row = case_when(
    Run == 1 & m_value == 5 ~ list(c(8.230893, 5.303274, 5, 6.666961, 8.828354, 0)),
    Run == 1 & m_value == 6 ~ list(c(6.363687, 4.89804, 5, 10.57402, 11.624766, 12.53241)),
    Run == 2 & m_value == 5 ~ list(c(8.58959, 7.161393, 5, 7.161393, 9.041251, 0)),
    Run == 2 & m_value == 6 ~ list(c(7.464666, 7.161393, 5, 7.959099, 9.717582, 11.081269))
  )) 

V0_J <- array(unlist(V0_values$V0_row), dim=c(n1_J,n2_J,max(m_J)))
Analysis_Data_J <- Analysis_Data_J %>%
  arrange(Round, UCB_ID, Run)
# y and y_star Arrays
y_J <- array(Analysis_Data_J$MyInvestment + 1, dim=c(n3_J,n1_J,n2_J))
y_star_J <- array(Analysis_Data_J$OpponentInvest + 1, dim=c(n3_J,n1_J,n2_J))

# pi Calculation
q<-10
N0<-1
# Initialize an empty list to store results
pi_J <- vector("list", n1_J * n2_J * n3_J)

index <- 1
for (i in 1:n1_J) {
  for (j in 1:n2_J) {
    for (t in 1:n3_J) {
      pi_J[[index]] <- sapply(1:6, function(k) {
        (m_J[i,j] - 1) + q * as.numeric((k - 1) >= y_star_J[t,i,j]) - (k - 1)
      })
      index <- index + 1
    }
  }
}

# Convert list to array if needed
pi_array <- array(unlist(pi_J), dim = c(n1_J, n2_J, n3_J, 6))

Analysis_Data_J <- Analysis_Data_J %>% arrange(UCB_ID)
dopamine_master_demographics_dec_23 <- dopamine_master_demographics_dec_23 %>% arrange(Participant.ID.)

x_J<-cbind(dopamine_master_demographics_dec_23$Age, as.numeric(dopamine_master_demographics_dec_23$Sex_assigned_at_birth=="Male"))
x_J[,1]<-scale(x_J[,1])

x_J<-cbind(rep(1, times=nrow(x_J)), x_J)
p_x_J<-ncol(x_J)


z_J <- array(0, dim = c(n1_J, n2_J, 2))

# Vectorized operation to fill z_J
z_J[,,1] <- w_J - 1  # Subtracting 1 from each element in w_J
z_J[,,2] <- r_J      # Directly assigning r_J values

p_z_J <- 2