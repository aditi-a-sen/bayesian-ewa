

# Documentation for Model Parameters

This documentation provides a comprehensive explanation of the code used to set up and prepare data for a Bayesian Experience-Weighted Attraction (EWA) model. The code is written in R and utilizes several libraries for data manipulation and analysis.

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code.

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Cmd+Shift+Enter*.

```{r}
print(1+1)
```

## 0. Packages Installation and Setup

First, load the required R packages and install dplyr if it's not already available:

```{r}
install.packages("dplyr",dependencies=TRUE)
library(dplyr)
library(stats)
library(tidyr)
library(readr)
```

This code loads three essential libraries:

dplyr: For data manipulation, readr: For reading CSV files, tidyr: For data tidying

Files and data frames used: patent_race , dvr_demo_df, dopamine_master_demographics_dec_23

Make sure to edit the code block below and replace the current pathname to reflect the location of the .csv files on your device.

```{r}
# Load data
Analysis_Data <- read_csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 2/patent_race_with_drug.csv")
dvr_demo_df <- read_csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 2/dvr_demo_df.csv")
dopamine_master_demographics_dec_23 <- read_csv("/Users/aditisen/Downloads/bayesian-ewa/attempt 2/dopamine_master_demographics_dec_23.csv")

```

It then imports three CSV files:

1.  patent_race_with_drug.csv: Contains the main analysis data
2.  dvr_demo_df.csv: Contains demographic information
3.  dopamine_master_demographics_dec_23.csv: Contains additional demographic data

### Data Cleaning

Ensure your dataset is properly cleaned and filtered before analysis. In this case, the dataset patent_race is filtered for relevant columns, and any missing data is removed using na.omit(). :

```{r}
# Data cleaning and filtering
Analysis_Data <- Analysis_Data %>%
  filter(SID %in% dvr_demo_df$UCB_ID) %>%
  filter(!SID %in% c(6022, 6100, 6103, 6144, 6171)) %>%
  mutate(Strong_Role_Indicator = ifelse(Condition == 's', 1, ifelse(Condition == 'w', 0, NA)))
```

This code performs the following operations:

1.  Filters Analysis_Data to include only subjects present in dvr_demo_df
2.  Removes specific subject IDs (6022, 6100, 6103, 6144, 6171)
3.  Creates a new column Strong_Role_Indicator:
    i.  1 for strong condition ('s')
    ii. 0 for weak condition ('w')
    iii. NA for any other condition

You can adjust this step for your experiment by removing specific subjects who are outliers or deemed unnecessary for the analysis and using the data frames with relevant information to the participant.

## 1. Model Parameter Set-up

We define key model parameters such as the number of subjects, runs, rounds, and sessions.

```{r}
# Calculate dimensions
n1_m7 <- n_distinct(Analysis_Data$SID)
n2_m7 <- n_distinct(Analysis_Data$Run)
n3_m7 <- n_distinct(Analysis_Data$Round)
n4_m7 <- n_distinct(Analysis_Data$Session)
```

These calculations determine:

1.  n1_m7: Number of unique subjects
2.  n2_m7: Number of unique runs
3.  n3_m7: Number of unique rounds
4.  n4_m7: Number of unique sessions

In this data set, there were 40 unique subjects, 2 runs, 2 sessions, and 80 rounds:

## 2. Array Creation for Model Parameters

Next, we set up matrices (m, w, r) to store calculated values based on subject performance and experimental conditions. This section creates arrays for various model parameters.

### 2.1 Calculating m_value

This creates a 3D array representing the number of choices each participant can make in each session and run. It's a fundamental input for the EWA model, providing the strategy space for each participant across all conditions.

Order the data by SID, Session, and Run

```{r}
# Ensure the data is properly ordered
Analysis_Data <- Analysis_Data %>%
  arrange(SID, Session, Run)

```

Calculate mean values for each SID, Session, and Run combination. This calculates the number of possible investments a participant can make. The m_value is essential for the EWA model to understand the range of choices available to each participant.

-   The code calculates mean_values for each unique combination of SID, Session, and Run.

-   m_value is set to 5 plus the mean of Strong_Role_Indicator. This represents the number of possible investments a participant can make.

```{r}
mean_values <- Analysis_Data %>%
  group_by(SID, Session, Run) %>%
  summarize(m_value = 5 + mean(Strong_Role_Indicator, na.rm = TRUE)) %>%
  ungroup()
```

Create a grid of all possible combinations. This ensures all potential combinations of Subject ID, Session, and Run are accounted for, even if they're not present in the original data. It's crucial for creating a complete data set without missing combinations, which could affect the model's accuracy.

```{r}
# Create a complete grid of all possible combinations
all_combinations <- expand.grid(
  SID = sort(unique(Analysis_Data$SID)),
  Session = 1:n4_m7,
  Run = 1:n2_m7
)
```

Join mean values with all combinations. This step fills in values for existing combinations and creates placeholders for missing ones.

```{r}
# Join the mean values with all combinations
filled_values <- all_combinations %>%
  left_join(mean_values, by = c("SID", "Session", "Run"))
```

Handle missing values and flips roles when necessary. This step handles missing values and flips roles when necessary. It's crucial for ensuring that the model has complete information about each participant's role and available strategies, even when data is missing. If an m_value is missing but others exist in the group, it's set to 11 minus the non-missing value. This effectively flips the role (strong to weak or vice versa). If all values are missing, NA is assigned.

```{r}
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
```

Create the m_m7 array with dimensions [subjects, sessions, runs]. Each cell m_m7[i,j,k] represents the number of choices participant i can make in session j, run k.

```{r}
# Create the m_m7 array
m_m7 <- array(
  filled_values$m_value,
  dim = c(n1_m7, n4_m7, n2_m7)
)
```

### 2.2 Calculating r_value, w_value, and drug_value

Similar to m_value, mean values are calculated for r_value (Strong_Role_Indicator), w_value (Session), and drug_value (Drug).

Order the data by SID, Session, and Run. This step ensures the data is properly ordered by Subject ID, Session, and Run. It's crucial for maintaining consistency in subsequent calculations and array creations, as the order of data affects how values are assigned in later steps.

```{r}
# Ensure the data is properly ordered
Analysis_Data <- Analysis_Data %>%
  arrange(SID, Session, Run)
```

Calculate mean values for r_m7, w_m7, and get drug values

```{r}
mean_values <- Analysis_Data %>%
  group_by(SID, Session, Run) %>%
  summarize(
    r_value = mean(Strong_Role_Indicator, na.rm = TRUE),
    w_value = first(Session),
    drug_value = first(Drug)
  ) %>%
  ungroup()

```

Create a complete grid of all possible combinations

```{r}

all_combinations <- expand.grid(
  SID = sort(unique(Analysis_Data$SID)),
  Session = 1:n4_m7,
  Run = 1:n2_m7
)

# Join the mean values with all combinations
filled_values <- all_combinations %>%
  left_join(mean_values, by = c("SID", "Session", "Run"))
```

Missing r_values are handled similarly to m_values, but flipped using 1 minus the non-missing value.

```{r}
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
```

Three arrays are created with dimensions [n1_m7, n4_m7, n2_m7].

-   r_m7[i,j,k] represents the Strong_Role_Indicator for participant i in session j, run k.

-   w_m7[i,j,k] represents the session number for participant i in session j, run k.

-   drug_m7[i,j,k] represents the drug condition for participant i in session j, run k.

```{r}
r_m7 <- array(filled_values$r_value, dim = c(n1_m7, n4_m7, n2_m7))
w_m7 <- array(filled_values$w_value, dim = c(n1_m7, n4_m7, n2_m7))
drug_m7 <- array(as.character(filled_values$drug_value), dim = c(n1_m7, n4_m7, n2_m7))
```

### 2.3 Calculating Session-Specific Values:

w_m7_j is a simplified version of w_m7, focusing only on the first and last sessions for each subject. This allows for easier comparison between initial and final sessions, which can be useful for analyzing changes over time.

```{r}
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
```

## 3. Initial Value Set-up for V0 (Investment Parameters)

This code creates the V0_m7 array with initial attraction values for each strategy, considering different combinations of roles, runs, and sessions. The values are predefined and a constant of 5 is added to each.

Here we set up the initial investment parameters (V0) for both weak and strong roles across two runs. V0 represents the initial attraction values for each strategy.

Recall that in this context, a strategy refers to the number of units a participant can choose to invest in a round. For example, a participant with a strong role has six strategies to choose from. They can invest anywhere from 0 to 5 tokens. Similarly, a weak participant can choose one of 5 strategies, investing anywhere from 0 to 4 tokens.

V0 values represent the initial attraction that a user can have to picking one strategy and this is represented as a decimal value in an array. 

> Run == 1 & m_value == 5 \~ list(c(8.230893, 5.303274, 5, 6.666961, 8.828354, 0))

For example, in the line above, if the Run is 1 and the participant has 5 strategies to choose from (weak role), the probability that they choose to invest 0 tokens is 8.230893, and the probability that they invest 5 tokens is 0 since they do not 5 tokens awarded.

We then attach these initial values to each individual game and put them in the array V0:

```{r}
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
```

More specifically, the value contained in V0[i,j,k]represents participant i’s initial attraction to pick strategy k in run j. This value is important because:

-   It influences the participant's initial probability of choosing that strategy.

-   It serves as a starting point for updating beliefs as the game progresses.

-   It can reflect prior experience or intuitions about the game.

## 5. Response Variables (y) Setup

The responses from the experiment are stored in two arrays: y for the subject's investment and y_star for the opponent's investment. The code below creates 3D arrays y and y_star from Analysis_Data\$MyInvestment and Analysis_Data\$OpponentInvest, then adds 1 to all values during array creation to specify the strategy number the participant picked. For example, if the participant invested 3 tokens, they picked strategy number 4 of all the choices available to them.

```{r}
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
```

## 5. Probability Distribution (pi) Calculation

A 4D array pi is initialized to store the calculated probabilities. The dimensions correspond to the number of subjects, runs, rounds, and the investment options (in this case, there are a maximum of 6 investment options). A nested loop iterates through each subject, run, and round to calculate the probabilities based on the model parameters (m and q). This utilizes the logic that influences investment decisions. The mathematics behind the logic are explained in the Attraction Updating section above. This calculates the payoff values for each strategy, which is a crucial component of the EWA model. 

\
The probability of participant i investing amount k during game t of run j =\> pi[i, j, t] = ](p(0), p(1)...p(6)): where m[i,j] -1 is the amount the participant invested + (reward \* (0 if the opponent invested more, if not, 1 since the participant i must have invested more) - the cost of playing (k units that were given up).

```{r}
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
```

## 7. Demographic Variables (x) Setup

This prepares covariates (age and sex) for the hierarchical model, allowing for the investigation of individual differences. The demographic data is extracted from the demographics data frame. The data includes age and sex, with sex being converted into a binary format (Male = 1 and Female = 0).

Next, the age column is scaled for 

-   Numerical stability: It helps prevent numerical instability in the model estimation process.

-   Interpretation: It puts Age on a comparable scale with other variables, making interpretation of coefficients easier.

-   Model convergence: It can improve model convergence and reduce computational time.

Scaling the age column involves these two steps:

-   Subtracting the mean age from the age in each row(centering)

-   Dividing by the standard deviation

After scaling, the column will have a mean of 0 and a standard deviation of 1. This process is also known as standardization or z-score normalization.

In the context of this model, scaling Age helps ensure that the magnitude of the Age variable doesn't overshadow other variables in the model and allows for more stable and interpretable parameter estimates in the Bayesian EWA model.

\
Lastly, an intercept term (a column of all 1’s) is added to the matrix to facilitate regression modeling. p_x is set to the number of columns in x, representing the number of covariates being considered.

```{r}
dopamine_master_demographics_dec_23 <- dopamine_master_demographics_dec_23 %>% 
  arrange(Participant.ID.)

x_m7 <- cbind(1, scale(dopamine_master_demographics_dec_23$Age), 
              as.numeric(dopamine_master_demographics_dec_23$Sex_assigned_at_birth == "Male"))
p_x_m7 <- ncol(x_m7)
```

## 8. Role Variables (z) Setup

The role variables represent run-specific information related to the investments. A 3D array z is initialized to store role information based on previous calculations (w and r). The number of role variables defined is the vector of participant and run-specific explanatory variables used in calculating the regression model.

Create a 3D array with the following dimensions

> z \<- array(0, dim = c(num_subjects, num_runs , 2))

Account for Session Effect:

> z [,,1] \<- w - 1

w_m7_j is a matrix indicating the session number (1 or 2) for each subject and run.

Subtracting 1 converts this to a 0/1 indicator, where:

-   0 represents the first session (reference category)

-   1 represents the second session

This allows the model to capture how behavior might change from the first to the second session.

Account for Role effect:

-   z [,,2] \<- r

-   r is a matrix indicating the role strength for each subject and run.

-   1 represents a "strong" role

-   0 represents a "weak" role

This allows the model to capture how behavior might differ between strong and weak roles.

By incorporating these role variables, the model can estimate how the session number and role strength affect the learning parameters, providing a more nuanced understanding of the participants' behavior across different experimental conditions.

```{r}
z_m7 <- array(0, dim = c(n1_m7, n4_m7, 2))
z_m7[,,1] <- w_m7_j - 1
z_m7[,,2] <- drug_m7_num
p_z_m7 <- 2
```
