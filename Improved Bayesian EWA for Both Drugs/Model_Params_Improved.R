library(dplyr)
library(readr)
library(tidyr)

# Load data
Analysis_Data <- read_csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 2/patent_race_with_drug.csv")
dvr_demo_df <- read_csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 2/dvr_demo_df.csv")
dopamine_master_demographics_dec_23 <- read_csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 2/dopamine_master_demographics_dec_23.csv")

# Data cleaning and filtering
Analysis_Data <- Analysis_Data %>%
  filter(SID %in% dvr_demo_df$UCB_ID) %>%
  filter(!SID %in% c(6022, 6100, 6103, 6144, 6171)) %>%
  mutate(Strong_Role_Indicator = ifelse(Condition == 's', 1, ifelse(Condition == 'w', 0, NA)))

# Calculate dimensions
n1_m7 <- n_distinct(Analysis_Data$SID)
n2_m7 <- n_distinct(Analysis_Data$Run)
n3_m7 <- n_distinct(Analysis_Data$Round)
n4_m7 <- n_distinct(Analysis_Data$Session)
############################
# Ensure the data is properly ordered
Analysis_Data <- Analysis_Data %>%
  arrange(SID, Session, Run)

# Calculate mean values
mean_values <- Analysis_Data %>%
  group_by(SID, Session, Run) %>%
  summarize(m_value = 5 + mean(Strong_Role_Indicator, na.rm = TRUE)) %>%
  ungroup()

# Create a complete grid of all possible combinations
all_combinations <- expand.grid(
  SID = sort(unique(Analysis_Data$SID)),
  Session = 1:n4_m7,
  Run = 1:n2_m7
)

# Join the mean values with all combinations
filled_values <- all_combinations %>%
  left_join(mean_values, by = c("SID", "Session", "Run"))

# Handle missing values and flip roles
filled_values <- filled_values %>%
  group_by(SID, Session) %>%
  mutate(
    m_value = case_when(
      !is.na(m_value) ~ m_value,
      is.na(m_value) & any(!is.na(m_value)) ~ 11 - m_value[!is.na(m_value)],
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

# Create the m_m7 array
m_m7 <- array(
  filled_values$m_value,
  dim = c(n1_m7, n4_m7, n2_m7)
)
############################
#---------------------------
# Ensure the data is properly ordered
Analysis_Data <- Analysis_Data %>%
  arrange(SID, Session, Run)

# Calculate mean values for r_m7, w_m7, and get drug values
mean_values <- Analysis_Data %>%
  group_by(SID, Session, Run) %>%
  summarize(
    r_value = mean(Strong_Role_Indicator, na.rm = TRUE),
    w_value = first(Session),
    drug_value = first(Drug)
  ) %>%
  ungroup()

# Create a complete grid of all possible combinations
all_combinations <- expand.grid(
  SID = sort(unique(Analysis_Data$SID)),
  Session = 1:n4_m7,
  Run = 1:n2_m7
)

# Join the mean values with all combinations
filled_values <- all_combinations %>%
  left_join(mean_values, by = c("SID", "Session", "Run"))

# Handle missing values for r_m7
filled_values <- filled_values %>%
  group_by(SID, Session) %>%
  mutate(
    r_value = case_when(
      !is.na(r_value) ~ r_value,
      is.na(r_value) & any(!is.na(r_value)) ~ 1 - r_value[!is.na(r_value)],
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

# Create the r_m7, w_m7, and drug_m7 arrays
r_m7 <- array(filled_values$r_value, dim = c(n1_m7, n4_m7, n2_m7))
w_m7 <- array(filled_values$w_value, dim = c(n1_m7, n4_m7, n2_m7))
drug_m7 <- array(as.character(filled_values$drug_value), dim = c(n1_m7, n4_m7, n2_m7))
#---------------------------

# Ensure the data is properly ordered
Analysis_Data <- Analysis_Data %>%
  arrange(SID, Session, Run)


session_values <- Analysis_Data %>%
  group_by(SID) %>%
  summarize(
    w_value_1 = first(Session),
    w_value_2 = last(Session),
    drug_1 = first(Drug),
    drug_2 = last(Drug),
    r_value_1 = mean(Strong_Role_Indicator[Session == min(Session)], na.rm = TRUE),
    r_value_2 = mean(Strong_Role_Indicator[Session == max(Session)], na.rm = TRUE)
  ) %>%
  ungroup()

w_m7_j <- cbind(session_values$w_value_1, session_values$w_value_2)
drug_m7_j <- cbind(session_values$drug_1, session_values$drug_2)
r_m7_j <- cbind(session_values$r_value_1, session_values$r_value_2)

drug_m7_num <- ifelse(drug_m7_j == "H", 1, 0)
# &&&&&&&&&&&&&&&&&&&&&
V0_m7<-array(NA, c(n1_m7, n4_m7, n2_m7, max(m_m7,na.rm=T)))
for(i in 1:n1_m7){
  for(j in 1:n2_m7){
    for(k in 1:n4_m7){
      #Weak Role, Run 1, Session 1
      if((j==1) & (m_m7[i,k,j]==5) & (k==1)){
        V0_m7[i,k,j,]<-(c(1.82, 0, 0.956, 0.956, 1.89, 0) + 5)
      } #add 0 to end when including weak+strong together
      #Strong Role, Run 1, Session 1
      if(j==1 & m_m7[i,k,j]==6 & (k==1)){
        V0_m7[i,k,j,]<-(c(0.405, 0, 0.511, 1.15, 1.34, 1.70) + 5)
      }
      #Weak Role, Run 2, Session 1
      if((j==2) & (m_m7[i,k,j]==5)& (k==1)){
        V0_m7[i,k,j,]<-(c(2.73, 1.20, 0, 1.10, 2.37, 0) + 5)
      } #add 0 to end when including weak+strong together
      #Strong Role, Run 2, Session 1
      if(j==2 & m_m7[i,k,j]==6  & (k==1)){
        V0_m7[i,k,j,]<-(c(0, 0.560, 0.223, 1.10, 1.70, 2.42) + 5)
      }
      #Weak Role, Run 1, Session 2
      if((j==1) & (m_m7[i,k,j]==5) & (k==2)){
        V0_m7[i,k,j,]<-(c(2.56, .336, 0, 0.182, 0.875, 0) + 5)
      } #add 0 to end when including weak+strong together
      #Strong Role, Run 1, Session 2
      if(j==1 & m_m7[i,k,j]==6 & (k==2)){
        V0_m7[i,k,j,]<-(c(0, 2.71, 2.08, 2.89, 2.48, 3.83) + 5)
      }
      #Weak Role, Run 2, Session 2
      if((j==2) & (m_m7[i,k,j]==5)  & (k==2)){
        V0_m7[i,k,j,]<-(c(2.71, 0, 0.811, 0.560, 1.61, 0) + 5)
      } #add 0 to end when including weak+strong together
      #Strong Role, Run 2, Session 2
      if(j==2 & m_m7[i,k,j]==6  & (k==2)){
        V0_m7[i,k,j,]<-(c(0.154, 0.405, 0, 0.693, 1.43, 1.79) + 5)
      }
    }
  }
}

# y and y_star calculation
# Create a data frame with all combinations, maintaining the correct order
all_combinations <- expand_grid(
  Round = sort(unique(Analysis_Data$Round)),
  SID = sort(unique(Analysis_Data$SID)),
  Session = sort(unique(Analysis_Data$Session)),
  Run = sort(unique(Analysis_Data$Run))
)

# Join the combinations with the original data
result_df <- all_combinations %>%
  left_join(Analysis_Data, by = c("Round", "SID", "Session", "Run"))

# Create index columns for each dimension
result_df <- result_df %>%
  mutate(
    t = match(Round, sort(unique(Analysis_Data$Round))),
    i = match(SID, sort(unique(Analysis_Data$SID))),
    k = match(Session, sort(unique(Analysis_Data$Session))),
    j = match(Run, sort(unique(Analysis_Data$Run)))
  )

# Create empty arrays
y_m7 <- array(NA, dim = c(n3_m7, n1_m7, n4_m7, n2_m7))
y_star_m7 <- array(NA, dim = c(n3_m7, n1_m7, n4_m7, n2_m7))

# Fill the arrays using vectorized operations
y_m7[cbind(result_df$t, result_df$i, result_df$k, result_df$j)] <- result_df$MyInvestment
y_star_m7[cbind(result_df$t, result_df$i, result_df$k, result_df$j)] <- result_df$OpponentInvest

# Add 1 to both arrays
y_m7 <- y_m7 + 1
y_star_m7 <- y_star_m7 + 1
# ---------------------------------------
# pi calculation
q <- 10
N0 <- 1

pi_m7 <- array(NA, c(n1_m7, n4_m7, n2_m7, n3_m7, max(m_m7, na.rm = TRUE)))
for(i in 1:n1_m7) {
  for(j in 1:n4_m7) {
    for(k in 1:n2_m7) {
      pi_m7[i,j,k,,] <- outer(y_star_m7[,i,j,k], 1:max(m_m7, na.rm = TRUE) - 1, 
                              function(y, g) (m_m7[i,j,k] - 1) + q * (g >= y) - g)
    }
  }
}

# x and p_x calculation
dopamine_master_demographics_dec_23 <- dopamine_master_demographics_dec_23 %>% 
  arrange(Participant.ID.)

x_m7 <- cbind(1, scale(dopamine_master_demographics_dec_23$Age), 
              as.numeric(dopamine_master_demographics_dec_23$Sex_assigned_at_birth == "Male"))
p_x_m7 <- ncol(x_m7)

# z and p_z calculation
z_m7 <- array(0, dim = c(n1_m7, n4_m7, 2))
z_m7[,,1] <- w_m7_j - 1
z_m7[,,2] <- drug_m7_num
p_z_m7 <- 2