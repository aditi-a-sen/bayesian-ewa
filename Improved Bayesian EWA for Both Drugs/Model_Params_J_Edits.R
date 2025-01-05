library(dplyr)
library(tidyr)
library(purrr)
library(stats)

library(readr)
library(rjags)
library(coda)
library(MCMCvis)
library(mcmcplots)

dvr_demo_df = read_csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 2/dvr_demo_df.csv")
dopamine_master_demographics_dec_23 = read_csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 2/dopamine_master_demographics_dec_23.csv")

# Load and preprocess data
load_and_preprocess_data <- function(file_path, ids_to_remove) {
  read_csv(file_path) %>%
    filter(SID %in% dvr_demo_df$UCB_ID, !(SID %in% ids_to_remove)) %>%
    mutate(Strong_Role_Indicator = case_when(
      Condition == 's' ~ 1,
      Condition == 'w' ~ 0,
      TRUE ~ NA_real_
    ))
}

Analysis_Data <- load_and_preprocess_data("/Users/aditisen/Downloads/bayesian-ewa/attempt 2/patent_race_with_drug.csv", c(6022, 6100, 6103, 6144, 6171))

n1_m7 <- length(unique(Analysis_Data$SID))
n2_m7 <- length(unique(Analysis_Data$Run))
n3_m7 <- length(unique(Analysis_Data$Round))
n4_m7 <- length(unique(Analysis_Data$Session))

# Create summary data
summary_data <- Analysis_Data %>%
  group_by(SID, Session, Run) %>%
  summarise(
    m = 5 + mean(Strong_Role_Indicator, na.rm = TRUE),
    r = mean(Strong_Role_Indicator, na.rm = TRUE),
    w = mean(Session, na.rm = TRUE),
    drug = first(Drug),
    .groups = 'drop'
  )

# Create arrays
create_arrays <- function(data) {
  list(
    m = array(data$m, dim = c(n1_m7, n4_m7, n2_m7)),
    r = array(data$r, dim = c(n1_m7, n4_m7, n2_m7)),
    w = array(data$w, dim = c(n1_m7, n4_m7, n2_m7)),
    drug = array(data$drug, dim = c(n1_m7, n4_m7, n2_m7))
  )
}

arrays <- create_arrays(summary_data)
m_m7 <- arrays$m

# Create V0 array
# Create V0 array
V0_values <- list(
  c(1.82, 0, 0.956, 0.956, 1.89, 0),
  c(0.405, 0, 0.511, 1.15, 1.34, 1.70),
  c(2.73, 1.20, 0, 1.10, 2.37, 0),
  c(0, 0.560, 0.223, 1.10, 1.70, 2.42),
  c(2.56, .336, 0, 0.182, 0.875, 0),
  c(0, 2.71, 2.08, 2.89, 2.48, 3.83),
  c(2.71, 0, 0.811, 0.560, 1.61, 0),
  c(0.154, 0.405, 0, 0.693, 1.43, 1.79)
)

get_V0_index <- function(run, m, session) {
  (run - 1) * 4 + (m - 5) * 2 + session
}

# Initialize V0_m7 as a numeric array
V0_m7 <- array(0, dim = c(n1_m7, n4_m7, n2_m7, max(arrays$m, na.rm = TRUE)))

# Fill V0_m7 with the correct values
for (i in 1:n1_m7) {
  for (j in 1:n4_m7) {
    for (k in 1:n2_m7) {
      index <- get_V0_index(k, arrays$m[i,j,k], j)
      V0_m7[i,j,k,1:length(V0_values[[index]])] <- V0_values[[index]] + 5
    }
  }
}

# Create y and y_star arrays
create_y_arrays <- function(data) {
  data %>%
    mutate(
      MyInvestment = MyInvestment + 1,
      OpponentInvest = OpponentInvest + 1
    ) %>%
    arrange(Round, SID, Session, Run) %>%
    summarise(
      y = list(array(MyInvestment, dim = c(n3_m7, n1_m7, n4_m7, n2_m7))),
      y_star = list(array(OpponentInvest, dim = c(n3_m7, n1_m7, n4_m7, n2_m7)))
    )
}

y_arrays <- create_y_arrays(Analysis_Data)
y_m7 <- y_arrays$y[[1]]
y_star_m7 <- y_arrays$y_star[[1]]

q <- 10
N0 <- 1

# Create pi array
create_pi_array <- function(m, y_star) {
  array_dim <- c(dim(m), n3_m7, max(m))
  
  expand.grid(
    i = 1:n1_m7, j = 1:n4_m7, k = 1:n2_m7, t = 1:n3_m7, g = 1:max(m)
  ) %>%
    mutate(
      pi = ((m[cbind(i,j,k)] - 1) + q * as.numeric((g - 1) >= y_star[cbind(t,i,j,k)]) - (g - 1))
    ) %>%
    pull(pi) %>%
    array(dim = array_dim)
}

pi_m7 <- create_pi_array(arrays$m, y_star_m7)

# vectorize the remaining code for x_m7 and z_m7 
x_m7 <- dopamine_master_demographics_dec_23 %>%
  arrange(Participant.ID.) %>%
  transmute(
    intercept = 1,
    age = scale(Age),
    sex = as.numeric(Sex_assigned_at_birth == "Male")
  ) %>%
  as.matrix()

p_x_m7 <- ncol(x_m7)

# Create z_m7
z_m7 <- summary_data %>%
  select(SID, Session, w, drug) %>%
  mutate(
    session_minus_one = w - 1,
    drug_numeric = ifelse(drug == "H", 1, 0)
  ) %>%
  select(SID, Session, session_minus_one, drug_numeric) %>%
  pivot_wider(
    names_from = Session,
    values_from = c(session_minus_one, drug_numeric),
    names_glue = "Session{Session}_{.value}"
  ) %>%
  select(-SID) %>%
  as.matrix()

p_z_m7 <- 2

# Reshape z_m7 into a 3D array
z_m7 <- array(z_m7, dim = c(n1_m7, n4_m7, 2))

z_m7_data <- unlist(z_m7)

# Create the array with the correct dimensions
z_m7 <- array(z_m7_data, dim = c(n1_m7, n4_m7, 2))
m_m7 <- array(as.numeric(m_m7), dim = dim(m_m7))

# Extract numeric data from V0_m7 if it's a list
#V0_m7_data <- unlist(V0_m7)
# Create the array with the correct dimensions
# V0_m7 <- array(as.numeric(V0_m7_data), dim = V0_m7_dim)