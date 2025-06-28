#################################################################################
###                                                                           ### 
###    Estimating Add-On Effects with the IPW Estimator                       ###
###                                                                           ### 
#################################################################################

# Estimand: Difference in total opioid dosage during follow-up
# Regimes under comparison: add-on-1 vs add-on-0
# Estimator: IPW estimator
# Follow-up: KK
# Treatment period: kappa

#################################################################################
###     Prerequisites                                                         ###
#################################################################################

# Load required packages
library(gfoRmula)
library(tidyverse)
library(purrr)
library(data.table)
library(survminer)
library(Hmisc)
library(survival)
library(boot)
library(splitstackshape)
library(speedglm)

# Disable scientific notation in output
options(scipen = 999)

# Import data
rm(list = ls())
setwd("N:/durable/IBMS Biostatistics/ACS/NTR+/Data")
data_pre <- readRDS("N:/durable/IBMS Biostatistics/ACS/NTR+/Data/G.A.2.5_catta.rds")

# Set predefined parameters
KK <- 21      # End of follow-up (time points 0 to KK)
kappa <- 1    # End of treatment period (0 to kappa), with kappa < KK
BB <- 200     # Number of bootstrap samples

#################################################################################
###     Prepare Data Set                                                      ###
#################################################################################

# Filter data to include only relevant time points
data_0 <- data_pre %>% filter(time < KK + 1)

# Create lagged versions of NSAID indicator
data_0 <- data_0 %>%
  group_by(id) %>%
  mutate(
    lag1_nsaid_ind = lag(nsaid_ind, default = 0),
    lag2_nsaid_ind = lag(lag1_nsaid_ind, default = 0)
  ) %>%
  ungroup()

# Restrict to individuals alive at time 0
id_death <- data_0 %>% filter(time == 0, death_ind == 0)
data_0 <- data_0 %>% filter(id %in% id_death$id)

# Number of eligible individuals
length(unique(data_0$id)) 

# Save outcome data for later use
dataY_0 <- data_0

#################################################################################
###     Add-on Effects                                                        ###
#################################################################################

# Data
data <- data_0
dataY <- dataY_0

# Fit logistic regression for propensity scores
model_pr <- speedglm(
  nsaid_ind ~
    # Baseline covariates:
    age + I(age^2) + sex + RHF + trm_cent + Kommuneindex +
    pt_asa_preinjury + acc_transport + acc_fall + acc_work +
    inj_mechanism + hosp_care_level + hosp_icu_days +
    hosp_los_days + res_gos_dischg + inntektGroup +
    ais_group + opioid_pre_dis + nsaid_pre_dis +
    # Time:
    time + I(time^2) +
    # Time-varying covariates:
    healthcare_use + opioid_omeq + B01A_ind + hosp_days +
    # Treatment:
    lag1_nsaid_ind + lag2_nsaid_ind,
  data = data,
  family = binomial()
)

# Predict propensity scores
pr <- predict(model_pr, newdata = data, type = "response")

data_pr <- data %>%
  add_column(pr_1 = pr) %>%
  mutate(
    pr_a = ifelse(nsaid_ind == 1, pr_1, 1 - pr_1),
    pr_0 = 1 - pr_1
  ) %>%
  group_by(id) %>%
  relocate(id, time, pr_a, pr_1, pr_0)
data <- data_pr

# Compute indicators and IP weights
data_pre_weights <- data %>%
  mutate(
    ind_a0_c0 = ifelse((opioid_omeq > 0 & nsaid_ind == 0) |
                         (opioid_omeq == 0 & nsaid_ind == 0), 1, 0),
    ind_a1_c0 = ifelse((opioid_omeq > 0 & nsaid_ind == 0) |
                         (opioid_omeq == 0 & nsaid_ind == 1), 1, 0),
    ind_a0_c1 = ifelse((opioid_omeq > 0 & nsaid_ind == 1) |
                         (opioid_omeq == 0 & nsaid_ind == 0), 1, 0),
    ind_a1_c1 = ifelse((opioid_omeq > 0 & nsaid_ind == 1) |
                         (opioid_omeq == 0 & nsaid_ind == 1), 1, 0)
  ) %>%
  mutate(
    w_a0_c0 = (pr_0 / pr_a) * ind_a0_c0,
    w_a1_c0 = (pr_1 / pr_a) * ind_a1_c0,
    w_a0_c1 = (pr_0 / pr_a) * ind_a0_c1,
    w_a1_c1 = (pr_1 / pr_a) * ind_a1_c1
  ) %>%
  dplyr::select(id, time, nsaid_ind, opioid_omeq,
                w_a0_c0, w_a1_c0, w_a0_c1, w_a1_c1)
data <- data_pre_weights

# Keep only time points within treatment period
data_time <- data %>% filter(time < kappa + 1)
data <- data_time

# Compute full time-dependent weights under add-on-0 and add-on-1 regimes

# Function to compute final weights for a single individual
compute_weights <- function(subset_data) {
  
  # Length of dataset
  mm <- dim(subset_data[, 1])[1]
  
  # Convert data frame to matrix for efficiency
  f_values_c0 <- as.matrix(subset_data[, c("w_a0_c0", "w_a1_c0")])
  f_values_c1 <- as.matrix(subset_data[, c("w_a0_c1", "w_a1_c1")])
  
  # Generate all combinations of 0/1 for mm time points
  combinations <- expand.grid(replicate(mm, c(0, 1), simplify = FALSE))
  
  # Compute product of weights for each combination
  compute_product_c0 <- function(row) {
    indices <- row + 1  # shift from 0/1 to 1/2 indexing
    prod(f_values_c0[cbind(1:mm, indices)])
  }
  compute_product_c1 <- function(row) {
    indices <- row + 1
    prod(f_values_c1[cbind(1:mm, indices)])
  }
  
  # Sum over all combinations
  result_c0 <- apply(combinations, 1, compute_product_c0)
  result_c1 <- apply(combinations, 1, compute_product_c1)
  
  c(sum(result_c0), sum(result_c1))
}

# Create matrix to store weights
data_w <- data.frame(
  id = data$id,
  time = data$time,
  Y = data$opioid_omeq,
  w_c0 = rep(NA, length(unique(data$id))),
  w_c1 = rep(NA, length(unique(data$id)))
)

# Compute weights for each individual and time
id_vec <- unique(data$id)
nn <- length(id_vec)

for (ii in 1:nn) {
  for (jj in 0:kappa) {
    subset_data <- data %>% filter(id == id_vec[ii] & time < jj + 1)
    data_w[jj + 1 + (ii - 1) * (kappa + 1), 4] <- compute_weights(subset_data)[1]
    data_w[jj + 1 + (ii - 1) * (kappa + 1), 5] <- compute_weights(subset_data)[2]
  }
  print(ii)
}

data_w_fix <- data_w %>% mutate(time = as.numeric(time))
data_w <- as_tibble(data_w_fix) %>% group_by(id) %>% ungroup()

# Create empty data frames to store estimates
results_0 <- data.frame(time = 1:KK, estimate = rep(NA, KK))
results_1 <- data.frame(time = 1:KK, estimate = rep(NA, KK))
results_diff <- data.frame(time = 1:KK, estimate = rep(NA, KK))

# Compute weighted expectations of Y_k for k = 1, ..., KK
for (jj in 1:KK) {
  
  # Filter data at time jj
  data_Y <- dataY %>%
    ungroup() %>%
    filter(time == jj) %>%
    mutate(Y = opioid_omeq) %>%
    dplyr::select(id, Y)
  
  data_merge <- data_w %>%
    dplyr::select(-Y) %>%
    left_join(data_Y, by = c("id")) %>%
    filter(time == min(jj - 1, kappa))
  
  # Save estimates
  results_0[jj, 2] <- weighted.mean(data_merge$Y, data_merge$w_c0)
  results_1[jj, 2] <- weighted.mean(data_merge$Y, data_merge$w_c1)
  results_diff[jj, 2] <- weighted.mean(data_merge$Y, data_merge$w_c1) -
    weighted.mean(data_merge$Y, data_merge$w_c0)
}

results_cum <- sum(results_diff$estimate[1:KK])  

# Save IPW estimates
addon_0 <- results_0
addon_1 <- results_1
addon_diff <- results_diff
addon_cum <- results_cum

#################################################################################
###     IPW with bootstrap                                                    ###
#################################################################################

# Data
data <- data_0

# Bootstrap specifications
id <- unique(dataY$id)

# Create matrices to store estimates
results_0_re <- data.frame(time = 1:KK, estimate = addon_0$estimate)
results_1_re <- data.frame(time = 1:KK, estimate = addon_1$estimate)
results_diff_re <- data.frame(time = 1:KK, estimate = addon_diff$estimate)
results_cum_re <- rep(NA, BB)

# Define function to compute weights (same as above)
compute_weights <- function(subset_data) {
  
  mm <- dim(subset_data[, 1])[1]
  
  f_values_c0 <- as.matrix(subset_data[, c("w_a0_c0", "w_a1_c0")])
  f_values_c1 <- as.matrix(subset_data[, c("w_a0_c1", "w_a1_c1")])
  
  combinations <- expand.grid(replicate(mm, c(0, 1), simplify = FALSE))
  
  compute_product_c0 <- function(row) {
    indices <- row + 1
    prod(f_values_c0[cbind(1:mm, indices)])
  }
  compute_product_c1 <- function(row) {
    indices <- row + 1
    prod(f_values_c1[cbind(1:mm, indices)])
  }
  
  result_c0 <- apply(combinations, 1, compute_product_c0)
  result_c1 <- apply(combinations, 1, compute_product_c1)
  
  c(sum(result_c0), sum(result_c1))
}

for (bb in 1:BB) {
  
  # Create resampled dataset
  id_resample <- sample(id, size = length(id), replace = TRUE)
  data_re <- data[data$id %in% id_resample, ]
  dataY_re <- data_re
  
  # Estimate propensity scores
  model_pr <- speedglm(
    nsaid_ind ~
      age + I(age^2) + sex + RHF + trm_cent + Kommuneindex +
      pt_asa_preinjury + acc_transport + acc_fall + acc_work +
      inj_mechanism + hosp_care_level + hosp_icu_days +
      hosp_los_days + res_gos_dischg + inntektGroup +
      ais_group + opioid_pre_dis + nsaid_pre_dis +
      time + I(time^2) +
      healthcare_use + opioid_omeq + B01A_ind + hosp_days +
      lag1_nsaid_ind + lag2_nsaid_ind,
    data = data_re,
    family = binomial()
  )
  
  # Add propensity scores and IP weights
  pr <- predict(model_pr, newdata = data_re, type = "response")
  
  data_pr_re <- data_re %>%
    add_column(pr_1 = pr) %>%
    mutate(
      pr_a = ifelse(nsaid_ind == 1, pr_1, 1 - pr_1),
      pr_0 = 1 - pr_1
    ) %>%
    group_by(id) %>%
    mutate(
      w_s0 = ifelse(cumsum(nsaid_ind) == 0,
                    cumprod(1 / pr_0), 0),
      w_s1 = ifelse(cumsum(nsaid_ind) == time + 1,
                    cumprod(1 / pr_1), 0)
    ) %>%
    relocate(id, time, pr_a, pr_1, pr_0, w_s0, w_s1)
  data_re <- data_pr_re
  
  data_pre_weights_re <- data_re %>%
    mutate(
      ind_a0_c0 = ifelse((opioid_omeq > 0 & nsaid_ind == 0) |
                           (opioid_omeq == 0 & nsaid_ind == 0), 1, 0),
      ind_a1_c0 = ifelse((opioid_omeq > 0 & nsaid_ind == 0) |
                           (opioid_omeq == 0 & nsaid_ind == 1), 1, 0),
      ind_a0_c1 = ifelse((opioid_omeq > 0 & nsaid_ind == 1) |
                           (opioid_omeq == 0 & nsaid_ind == 0), 1, 0),
      ind_a1_c1 = ifelse((opioid_omeq > 0 & nsaid_ind == 1) |
                           (opioid_omeq == 0 & nsaid_ind == 1), 1, 0)
    ) %>%
    mutate(
      w_a0_c0 = (pr_0 / pr_a) * ind_a0_c0,
      w_a1_c0 = (pr_1 / pr_a) * ind_a1_c0,
      w_a0_c1 = (pr_0 / pr_a) * ind_a0_c1,
      w_a1_c1 = (pr_1 / pr_a) * ind_a1_c1
    ) %>%
    dplyr::select(id, time, nsaid_ind, opioid_omeq,
                  w_s0, w_s1,
                  w_a0_c0, w_a1_c0, w_a0_c1, w_a1_c1)
  data_re <- data_pre_weights_re
  
  # Keep time points within treatment period
  data_time_re <- data_re %>% filter(time < kappa + 1)
  data_re <- data_time_re
  
  # Create matrix to store weights
  data_w_re <- data.frame(
    id = data_re$id,
    time = data_re$time,
    Y = data_re$opioid_omeq,
    w_c0 = rep(NA, length(unique(data_re$id))),
    w_c1 = rep(NA, length(unique(data_re$id)))
  )
  
  # Compute weights
  id_vec <- unique(data_re$id)
  nn <- length(id_vec)
  for (ii in 1:nn) {
    for (jj in 0:kappa) {
      subset_data <- data_re %>% filter(id == id_vec[ii] & time < jj + 1)
      data_w_re[jj + 1 + (ii - 1) * (kappa + 1), 4] <- compute_weights(subset_data)[1]
      data_w_re[jj + 1 + (ii - 1) * (kappa + 1), 5] <- compute_weights(subset_data)[2]
    }
  }
  
  # Estimate weighted expectations of Y_k
  for (jj in 1:KK) {
    data_Y_re <- dataY_re %>%
      ungroup() %>%
      filter(time == jj) %>%
      mutate(Y = opioid_omeq) %>%
      dplyr::select(id, Y)
    data_merge_re <- data_w_re %>%
      dplyr::select(-Y) %>%
      left_join(data_Y_re, by = c("id")) %>%
      filter(time == min(jj - 1, kappa))
    results_0_re[jj, bb + 2] <- weighted.mean(data_merge_re$Y, data_merge_re$w_c0)
    results_1_re[jj, bb + 2] <- weighted.mean(data_merge_re$Y, data_merge_re$w_c1)
    results_diff_re[jj, bb + 2] <- weighted.mean(data_merge_re$Y, data_merge_re$w_c1) -
      weighted.mean(data_merge_re$Y, data_merge_re$w_c0)
  }
  results_cum_re[bb] <- sum(results_diff_re[1:KK, bb + 2])
  print(bb)
}

# Save bootstrap results
addon_0_re <- results_0_re
addon_1_re <- results_1_re
addon_diff_re <- results_diff_re
addon_cum_re <- results_cum_re
