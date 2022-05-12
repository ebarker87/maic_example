################################ Script details ################################
## Conference: R for HTA 
## Title: MAIC-ing it easy: R for matching-adjusted indirect comparisons
## Date: 20th May 2022
## Author: Erin Barker
## Organisation: York Health Economics Consortium
#==============================================================================#

##### Notes #####
# IPD = individual patient data; t1 = trial 1; t2 = trial 2

# Load the relevant libraries
library("tidyverse"); library("maic"); library("sandwich")


#==============================================================================#
#### IPD ####
#==============================================================================#

# Read in IPD
t1_ipd <- read.csv("ipd.csv")

# Generate unadjusted mean difference and CI for change in HbA1c from baseline
t1_t_test <- t.test(t1_ipd$hba1c_change[t1_ipd$trt== "Treatment A"], t1_ipd$hba1c_change[t1_ipd$trt== "Placebo"])

# Extract t1 outcomes
t1_outcomes <- data.frame(
  treatment = "Treatment A vs. Placebo (unadjusted)",
  md = as.numeric(t1_t_test$estimate[1]) - as.numeric(t1_t_test$estimate[2]),
  se = (t1_t_test$conf.int[2] -t1_t_test$conf.int[1])/3.92,
  lci = t1_t_test$conf.int[1],
  uci = t1_t_test$conf.int[2]
)
remove(t1_t_test)


#==============================================================================#
#### Aggregate data ####
#==============================================================================#

# Read in the t2 aggregate data
t2_agg_data <- read.csv("aggregate_data.csv")

#### Outcomes from t2 ####
t2_outcomes <- t2_agg_data %>% 
  mutate(treatment = "Treatment B vs. Placebo",
         md = hba1c_change,
         se = hba1c_change,
         lci = hba1c_change_lci,
         uci = hba1c_change_uci
         ) %>% 
  select(treatment, md, lci, uci, se)


#==============================================================================#
#### Baseline characteristics ####
#==============================================================================#

# Compare the baseline characteristics
rbind(t1_ipd %>% 
  summarise(
    age = mean(age),
    hba1c_baseline = mean(hba1c_baseline),
    fpg = mean(fpg)
  ),
  t2_agg_data %>% 
    select (age, hba1c_baseline, fpg)
  )


#==============================================================================#
#### Weights ####
#==============================================================================#

# Define the data dictionary
maic_dictionary <- data.frame(
  "match.id" = c(
    "age",
    "hba1c_baseline",
    "fpg"
  ),
  "target.variable" = c(
    "age",
    "hba1c_baseline",
    "fpg"
  ),
  "index.variable" = c(
    "age",
    "hba1c_baseline",
    "fpg"
  ),
  "match.type" = c(
    "mean",
    "mean",
    "mean"
  )
)

# Construct the MAIC input matrix
MAIC_matrix <- createMAICInput(
  index = t1_ipd,
  target = t2_agg_data,
  dictionary = maic_dictionary,
  matching.variables = c(
    "age",
    "hba1c_baseline",
    "fpg"
  )
)

# Extract the MAIC weights
MAIC_weights <- maicWeight(MAIC_matrix)


#==============================================================================#
#### Evaluate the weights  ####
#==============================================================================#

# Calculate the adjusted baseline characteristics
baseline_summary <- reportCovariates(
  index = t1_ipd,
  target = t2_agg_data,
  dictionary = maic_dictionary,
  matching.variables = c(
    "age",
    "hba1c_baseline",
    "fpg"),
  weights = MAIC_weights
)

# Rescale the weights
MAIC_weights_RS <- (MAIC_weights / sum(MAIC_weights)) * 300

# Summary of weights
summary(MAIC_weights_RS)

# Plot the distribution of the weights
plot(MAIC_weights_RS, main = "Hisotgram of rescaled weights")

# Calculate the effective sample size (ESS)
sum(MAIC_weights)^2/sum(MAIC_weights^2)
# ESS is 239.9466


#==============================================================================#
#### Calculate the adjusted outcome ####
#==============================================================================#

# Add the weights to the IPD
t1_ipd <- t1_ipd %>%
  mutate(weights = MAIC_weights)

# Run a linear model
model_1 <- lm(hba1c_change ~ trt, data= t1_ipd, weights=weights)

# Extract the adjusted outcome 
t1_outcomes <- t1_outcomes %>% 
  add_row(
    treatment = "Treatment A vs. Placebo (adjusted)",
    md = coef(model_1)["trtTreatment A"],
    se = sqrt(vcovHC(model_1)["trtTreatment A", "trtTreatment A"]),
    lci = md - 1.96*se,
    uci = md + 1.96*se
  )


#==============================================================================#
#### Results ####
#==============================================================================#

# Generate unadjusted (native) and adjusted outcome for treatment A vs. treatment B 
MAIC_results <- rbind(
  data.frame(
    treatment = "Treatment A vs. Treatment B (unadjusted)",
    md = t1_outcomes$md[t1_outcomes$treatment == "Treatment A vs. Placebo (unadjusted)"] - t2_outcomes$md,
    se = sqrt(t1_outcomes$se[t1_outcomes$treatment == "Treatment A vs. Placebo (unadjusted)"]^2+t2_outcomes$se^2)
    ) %>% 
    mutate(lci = md - 1.96*se,
           uci = md + 1.96*se),
  data.frame(
    treatment = "Treatment A vs. Treatment B (adjusted)",
    md = t1_outcomes$md[t1_outcomes$treatment == "Treatment A vs. Placebo (adjusted)"] - t2_outcomes$md,
    se = sqrt(t1_outcomes$se[t1_outcomes$treatment == "Treatment A vs. Placebo (adjusted)"]^2+t2_outcomes$se^2)
    ) %>% 
    mutate(lci = md - 1.96*se,
           uci = md + 1.96*se)
)

row.names(MAIC_results) <- NULL
