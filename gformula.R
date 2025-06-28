#################################################################################
###                                                                           ### 
###    Estimating Add-On Effects with the G-Formula Estimator                 ###
###                                                                           ### 
#################################################################################

# Estimand: Difference in total counterfactual opioid dosage during follow-up
# Regimes under comparison: add-on-1 vs. add-on-0
# Estimator: G-formula estimator
# Follow-up period: 0 to KK
# Treatment period: 0 to kappa (with kappa < KK)

#################################################################################
###     Prerequisites                                                         ###
#################################################################################

# Load required packages
library(tidyverse)
library(gfoRmula)
library(purrr)
library(data.table)
library(survminer)
library(Hmisc)
library(survival)
library(boot)
library(dplyr)
library(splitstackshape)
library(speedglm)

# Disable scientific notation in printed output
options(scipen = 999)

# Clear environment and set working directory
rm(list = ls())
setwd("N:/durable/IBMS Biostatistics/ACS/NTR+/Data")

# Import pre-processed data
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

# Subset data to treatment period (0 to kappa)
dataA <- data_0 %>% filter(time < kappa + 1)

# Full data for outcome
dataY <- dataY_0

# Create empty data frames to store results
results_0     <- data.frame(time = 1:KK, estimate = rep(NA, KK))  # Intervention 0
results_1     <- data.frame(time = 1:KK, estimate = rep(NA, KK))  # Intervention 1
results_diff  <- data.frame(time = 1:KK, estimate = rep(NA, KK))  # Difference

# Loop over each follow-up time point jj (1 to KK)
for (jj in 1:KK) {
  
  # Filter treatment data up to time min(jj, kappa + 1)
  dataA_j <- dataA %>% filter(time < min(jj, kappa + 1))
  
  # Filter and reshape outcome data to include Y only at time == jj
  dataY_j <- dataY %>%
    filter(time < (jj + 1)) %>%
    mutate(Y = ifelse(time == jj, opioid_omeq, NA)) %>%
    dplyr::select(id, time, Y) %>%
    mutate(time = time - max(1, jj - kappa)) %>%
    filter(time > -1)
  
  # Merge treatment and outcome data
  data <- dataA_j %>%
    left_join(dataY_j, by = c("id", "time")) %>%
    relocate(id, time, Y)
  
  # Define dynamic treatment strategies under comparison
  # Add-on-0 regime
  add_on_0 <- function(newdf, pool, intvar, intvals, time_name, t) {
    newdf[opioid_omeq > 0, (intvar) := 0]
  }
  # Add-on-1 regime 
  add_on_1 <- function(newdf, pool, intvar, intvals, time_name, t) {
    newdf[opioid_omeq > 0, (intvar) := 1]
  }
  
  # List of interventions to pass to gformula()
  interventions <- list(list(c(add_on_0)), list(c(add_on_1)))
  
  # Run parametric g-formula estimator (without bootstrap)
  g.results <- gformula(
    seed = 54514,
    
    # Dataset and identifiers
    obs_data = data,
    id = "id",
    
    # Time
    time_name = "time",
    time_points = jj,
    
    # Outcome
    outcome_name = "Y",
    outcome_type = "continuous_eof",
    ymodel = Y ~
      # Baseline covariates
      age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
      acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level +
      hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group +
      opioid_pre_dis + nsaid_pre_dis +
      # Treatment and time
      nsaid_ind * time + nsaid_ind * I(time^2) + lag1_nsaid_ind +
      # Time-varying covariates
      healthcare_use + opioid_omeq + B01A_ind + hosp_days,
    
    # Covariates
    basecovs = c(
      "age", "sex", "RHF", "trm_cent", "Kommuneindex", "pt_asa_preinjury",
      "acc_transport", "acc_fall", "acc_work", "inj_mechanism", "hosp_care_level",
      "hosp_icu_days", "hosp_los_days", "res_gos_dischg", "inntektGroup",
      "ais_group", "opioid_pre_dis", "nsaid_pre_dis"
    ),
    
    covnames = c("nsaid_ind", "B01A_ind", "healthcare_use", "opioid_omeq", "hosp_days"),
    covtypes = c("binary", "binary", "normal", "normal", "normal"),
    
    # Lagged covariates and history
    histories = c(lagged),
    histvars = list(c("nsaid_ind", "B01A_ind", "healthcare_use", "opioid_omeq", "hosp_days")),
    
    covparams = list(
      covmodels = c(
        nsaid_ind ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex +
          pt_asa_preinjury + acc_transport + acc_fall + acc_work + inj_mechanism +
          hosp_care_level + hosp_icu_days + hosp_los_days + res_gos_dischg +
          inntektGroup + ais_group + opioid_pre_dis + nsaid_pre_dis +
          time + I(time^2) +
          healthcare_use + opioid_omeq + B01A_ind + hosp_days +
          lag1_nsaid_ind + lag2_nsaid_ind,
        
        B01A_ind ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex +
          pt_asa_preinjury + acc_transport + acc_fall + acc_work + inj_mechanism +
          hosp_care_level + hosp_icu_days + hosp_los_days + res_gos_dischg +
          inntektGroup + ais_group + opioid_pre_dis + nsaid_pre_dis +
          time + I(time^2) +
          lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days +
          lag1_nsaid_ind,
        
        healthcare_use ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex +
          pt_asa_preinjury + acc_transport + acc_fall + acc_work + inj_mechanism +
          hosp_care_level + hosp_icu_days + hosp_los_days + res_gos_dischg +
          inntektGroup + ais_group + opioid_pre_dis + nsaid_pre_dis +
          time + I(time^2) +
          lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days +
          lag1_nsaid_ind,
        
        opioid_omeq ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex +
          pt_asa_preinjury + acc_transport + acc_fall + acc_work + inj_mechanism +
          hosp_care_level + hosp_icu_days + hosp_los_days + res_gos_dischg +
          inntektGroup + ais_group + opioid_pre_dis + nsaid_pre_dis +
          time + I(time^2) +
          lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days +
          lag1_nsaid_ind,
        
        hosp_days ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex +
          pt_asa_preinjury + acc_transport + acc_fall + acc_work + inj_mechanism +
          hosp_care_level + hosp_icu_days + hosp_los_days + res_gos_dischg +
          inntektGroup + ais_group + opioid_pre_dis + nsaid_pre_dis +
          time + I(time^2) +
          lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days +
          lag1_nsaid_ind
      )
    ),
    
    # Intervention definitions
    intvars = list("nsaid_ind", "nsaid_ind"),  # Intervening on NSAIDs in both arms
    interventions = interventions,             # Use the list of custom strategies
    int_descript = c("Add-on-0", "Add-on-1"),  # Description of strategies
    ref_int = 1,                               # Reference: intervention 1 (add-on 0 regime)
    
    # Bootstrap settings
    nsamples = 0  # No bootstrap used here
  )
  
  # Store estimates
  results_0[jj, 2]    <- data.frame(g.results$result)[2, 4]
  results_1[jj, 2]    <- data.frame(g.results$result)[3, 4]
  results_diff[jj, 2] <- data.frame(g.results$result)[3, 6]
}

# Calculate cumulative effect estimate
results_cum <- sum(results_diff$estimate[1:KK])  

# Save estimates
addon_0   <- results_0
addon_1   <- results_1
addon_diff <- results_diff # The estimated add-on effect at each time-point 
addon_cum <- results_cum 
    
#################################################################################
###     Add-on Effects (with bootstrap)                                       ###
#################################################################################

# Data 
dataA <- data_0 %>% filter(time < kappa + 1) 
dataY <- dataY_0

# Empty matrixes to store results
results_0_re <- data.frame(time = 1:KK, estimate = results_0$estimate, matrix(NA, nrow = KK, ncol = BB))
results_1_re <- data.frame(time = 1:KK, estimate = results_1$estimate, matrix(NA, nrow = KK, ncol = BB))
results_diff_re <- data.frame(time = 1:KK, estimate = results_diff$estimate, matrix(NA, nrow = KK, ncol = BB))
results_cum_re <- c(rep(NA, BB))

id <- unique(dataY$id)
set.seed(123)  # Set an initial seed for reproducibility
seeds <- sample(1:100000, BB, replace = FALSE)  # Generate BB unique seeds

for (bb in 1:BB){
  
  id_resample <- sample(id, size = length(id), replace = TRUE)
  dataA_re <- dataA[dataA$id %in% id_resample, ]
  dataY_re <- dataY[dataY$id %in% id_resample, ]
  boot_seed <- seeds[bb]  
  
  for (jj in 1:KK){ # jj representing outcome at month jj
    
    # Create data set  
    dataA_j <- dataA_re %>% filter(time < min(jj, kappa + 1) )
    dataY_j <- dataY_re %>% 
      filter(time < (jj + 1)) %>%
      mutate(Y = ifelse(time == jj, opioid_omeq, NA)) %>%
      dplyr::select(id, time, Y) %>%
      mutate(time = time - max(1, jj - kappa)) %>%
      filter(time > -1)
    data <- dataA_j %>% left_join(dataY_j, by = c("id", "time")) %>% relocate(id, time, Y)
    
    # Specifying the dynamic treatment strategies under comparison
    add_on_0 <- function(newdf, pool, intvar, intvals, time_name, t){newdf[opioid_omeq > 0, (intvar):=0]}
    add_on_1 <- function(newdf, pool, intvar, intvals, time_name, t){newdf[opioid_omeq > 0, (intvar):=1]}
    
    # Specify both of the custom functions above in the list of interventions
    interventions <- list(list(c(add_on_0)), list(c(add_on_1)))
    
    # Specifying parameters for implementing parametric g-formula w. bootstrap  
    g.results <- gformula(
      
      seed = boot_seed,
      
      # Dataset and participant identifiers
      obs_data = data,
      id = "id",
      
      # Time
      time_name = "time",
      time_points = jj, 
      
      # Outcome
      outcome_name = "Y",
      outcome_type = "continuous_eof",
      ymodel = Y ~ 
        # Baseline covariates: 
        age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury + 
        acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level + hosp_icu_days +
        hosp_los_days + res_gos_dischg + inntektGroup + ais_group + opioid_pre_dis + nsaid_pre_dis +
        # Treatment and time:
        nsaid_ind*time + nsaid_ind*I(time^2) + lag1_nsaid_ind +
        # Time-varying covariates:
        healthcare_use + opioid_omeq + B01A_ind + hosp_days,  
      
      # Covariates
      
      # Time-fixed 
      basecovs = c("age", "sex", "RHF", "trm_cent", "Kommuneindex", "pt_asa_preinjury", 
                   "acc_transport", "acc_fall", "acc_work", "inj_mechanism", "hosp_care_level", 
                   "hosp_icu_days", "hosp_los_days", "res_gos_dischg", "inntektGroup", "ais_group", 
                   "opioid_pre_dis", "nsaid_pre_dis"),
      
      # Time-varying   
      covnames = c("nsaid_ind", "B01A_ind", "healthcare_use", "opioid_omeq", "hosp_days"), 
      covtypes = c("binary", "binary", "normal", "normal", "normal"),
      
      histories = c(lagged),  # Custom lag function defined outside; consider documenting clearly
      
      histvars = list(c("nsaid_ind", "B01A_ind", "healthcare_use", "opioid_omeq", "hosp_days")),
      
      covparams = list(covmodels = c(
        # Consider checking if model coefficients are being estimated stably across all bootstrap samples
        nsaid_ind ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
          acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level + 
          hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group + 
          opioid_pre_dis + nsaid_pre_dis + time + I(time^2) + 
          healthcare_use + opioid_omeq + B01A_ind + hosp_days + 
          lag1_nsaid_ind + lag2_nsaid_ind,
        B01A_ind ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
          acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level + 
          hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group + 
          opioid_pre_dis + nsaid_pre_dis + time + I(time^2) + 
          lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days +
          lag1_nsaid_ind,
        healthcare_use ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
          acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level + 
          hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group + 
          opioid_pre_dis + nsaid_pre_dis + time + I(time^2) + 
          lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days +
          lag1_nsaid_ind,
        opioid_omeq ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
          acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level + 
          hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group + 
          opioid_pre_dis + nsaid_pre_dis + time + I(time^2) + 
          lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days +
          lag1_nsaid_ind,
        hosp_days ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
          acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level + 
          hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group + 
          opioid_pre_dis + nsaid_pre_dis + time + I(time^2) + 
          lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days +
          lag1_nsaid_ind
      )),
      
      # Intervention
      intvars = list("nsaid_ind", "nsaid_ind"),         # Ensure the order matches 'interventions' list
      interventions = interventions,
      int_descript = c("Add-on-0", "Add-on-1"),          
      ref_int = 1,                                      # Add-on-0 regime (first in list) is reference. 
      
      # Bootstrapping parameters
      nsamples = 0   
      
    )
    
    # Save results 
    results_0_re[jj, bb + 2] <- data.frame(g.results$result)[2, 4] 
    results_1_re[jj, bb + 2] <- data.frame(g.results$result)[3, 4] 
    results_diff_re[jj, bb + 2] <- data.frame(g.results$result)[3, 6]  
  }
  
  # Save cumulative contrast over time
  results_cum_re[bb] <- sum(results_diff_re[1:KK, bb + 2])  
  
}

### Save bootstrap estimates
addon_0_re <- results_0_re
addon_1_re <- results_1_re
addon_diff_re <- results_diff_re # The estimated add-on effect at each time-point 
addon_cum_re <- results_cum_re

#################################################################################
###     Effects under static regimes                                          ###
#################################################################################

# Data
dataA <- data_0 %>% filter(time < kappa + 1) 
dataY <- dataY_0

# Create empty matrixes to save results  
results_0 <- data.frame(time = c(1:KK), estimate = rep(NA, KK))  
results_1 <- data.frame(time = c(1:KK), estimate = rep(NA, KK))
results_diff <- data.frame(time = c(1:KK), estimate = rep(NA, KK))

for (jj in 1:KK){ # jj representing outcome at month jj
  
  # ---------------------------------------------- #
  # Prepare data up to time jj (truncated by kappa)
  # ---------------------------------------------- #
  
  dataA_j <- dataA %>% filter(time < min(jj, kappa + 1) )
  
  dataY_j <- dataY %>% 
    filter(time < (jj+1)) %>% 
    mutate(Y = ifelse(time == jj, opioid_omeq, NA)) %>% 
    dplyr::select(id, time, Y) %>%
    mutate(time = time - max(1, jj - kappa)) %>% 
    filter(time > -1)
  
  data <- dataA_j %>% 
    left_join(dataY_j, by = c("id", "time")) %>% 
    relocate(id, time, Y)
  
  # ----------------------------------------- #
  # Define static interventions for g-formula #
  # ----------------------------------------- #
  
  interventions <- list(
    list(c(static, rep(0, min(jj, kappa + 1)))), 
    list(c(static, rep(1, min(jj, kappa + 1))))
  )
  
  # ------------------------------------------------------ #
  # Run parametric g-formula for each intervention setting #
  # ------------------------------------------------------ #
  
  g.results <- gformula(
    
    seed = 54514,  # Set random seed for reproducibility
    
    # ----------------------------------------- #
    # -- Dataset and participant identifiers -- #
    # ----------------------------------------- #
    
    obs_data = data,
    id = "id",
    
    # ---------- #    
    # -- Time -- #    
    # ---------- # 
    
    time_name = "time",
    time_points = jj, 
    
    # ------------- #
    # -- Outcome -- #
    # ------------- #   
    
    outcome_name = "Y",
    outcome_type = "continuous_eof",
    ymodel = Y ~ 
      # Baseline covariates:
      age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury + 
      acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level + 
      hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group + 
      opioid_pre_dis + nsaid_pre_dis +
      # Treatment and time:
      nsaid_ind*time + nsaid_ind*I(time^2) + lag1_nsaid_ind +
      # Time-varying covariates:
      healthcare_use + opioid_omeq + B01A_ind + hosp_days,  
    
    # ---------------- #
    # -- Covariates -- #
    # ---------------- #  
    
    basecovs = c("age", "sex", "RHF", "trm_cent", "Kommuneindex", 
                 "pt_asa_preinjury", "acc_transport", "acc_fall", 
                 "acc_work", "inj_mechanism", "hosp_care_level", 
                 "hosp_icu_days", "hosp_los_days", "res_gos_dischg", 
                 "inntektGroup", "ais_group", "opioid_pre_dis", 
                 "nsaid_pre_dis"),
    
    covnames = c("nsaid_ind", "B01A_ind", "healthcare_use", "opioid_omeq", "hosp_days"), 
    covtypes = c("binary", "binary", "normal", "normal", "normal"),
    
    histories = c(lagged),
    histvars = list(c("nsaid_ind", "B01A_ind", "healthcare_use", "opioid_omeq", "hosp_days")),
    
    covparams = list(covmodels = c(
      nsaid_ind ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
        acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level + 
        hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group + 
        opioid_pre_dis + nsaid_pre_dis + time + I(time^2) + 
        healthcare_use + opioid_omeq + B01A_ind + hosp_days + 
        lag1_nsaid_ind + lag2_nsaid_ind,
      B01A_ind ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
        acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level + 
        hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group + 
        opioid_pre_dis + nsaid_pre_dis + time + I(time^2) + 
        lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days + 
        lag1_nsaid_ind,
      healthcare_use ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
        acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level + 
        hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group + 
        opioid_pre_dis + nsaid_pre_dis + time + I(time^2) + 
        lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days + 
        lag1_nsaid_ind,
      opioid_omeq ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
        acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level + 
        hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group + 
        opioid_pre_dis + nsaid_pre_dis + time + I(time^2) + 
        lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days + 
        lag1_nsaid_ind,
      hosp_days ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
        acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level + 
        hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group + 
        opioid_pre_dis + nsaid_pre_dis + time + I(time^2) + 
        lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days + 
        lag1_nsaid_ind
    )),
    
    # ------------------ #
    # -- Intervention -- #
    # ------------------ #
    
    intvars = list("nsaid_ind", "nsaid_ind"),         
    interventions = interventions,                    
    int_descript = c("Add-on-0", "Add-on-1"),         
    ref_int = 1,                                    
    
    # ------------------------------ #
    # -- Bootstrapping parameters -- #
    # ------------------------------ #
    
    nsamples = 0   # No bootstrap samples
    
  )
  
  # ------------------------ #
  # Extract and store result #
  # ------------------------ #
  
  results_0[jj, 2] <- data.frame(g.results$result)[2, 4] 
  results_1[jj, 2] <- data.frame(g.results$result)[3, 4] 
  results_diff[jj, 2] <- data.frame(g.results$result)[3, 6]  
  
}

# -------------------------------------------- #
# Calculate cumulative difference across months
# -------------------------------------------- #

results_cum <- sum(results_diff$estimate[1:KK])   

# Save estimates 
static_0 <- results_0
static_1 <- results_1
static_diff <- results_diff
static_cum <- results_cum

#################################################################################
###     Effects under static regimes (with bootstrap)                         ###
#################################################################################

# Load data
dataA <- data_0 %>% filter(time < kappa + 1)
dataY <- dataY_0

# Matrices to store results
results_0_re     <- data.frame(time = 1:KK, estimate = results_0$estimate, matrix(NA, nrow = KK, ncol = BB))
results_1_re     <- data.frame(time = 1:KK, estimate = results_1$estimate, matrix(NA, nrow = KK, ncol = BB))
results_diff_re  <- data.frame(time = 1:KK, estimate = results_diff$estimate, matrix(NA, nrow = KK, ncol = BB))
results_cum_re   <- rep(NA, BB)

# IDs and bootstrap seeds
id <- unique(dataY$id)
set.seed(123)
seeds <- sample(1:100000, BB, replace = FALSE)

# Bootstrap loop
for (bb in 1:BB) {
  
  # Resample individuals with replacement
  id_resample <- sample(id, size = length(id), replace = TRUE)
  dataA_re <- dataA %>% filter(id %in% id_resample)
  dataY_re <- dataY %>% filter(id %in% id_resample)
  boot_seed <- seeds[bb]
  
  # Loop over each time point
  for (jj in 1:KK) {
    
    # Prepare analysis dataset up to time jj
    dataA_j <- dataA_re %>% filter(time < min(jj, kappa + 1))
    
    dataY_j <- dataY_re %>%
      filter(time < (jj + 1)) %>%
      mutate(Y = ifelse(time == jj, opioid_omeq, NA)) %>%
      select(id, time, Y) %>%
      mutate(time = time - max(1, jj - kappa)) %>%
      filter(time > -1)
    
    data <- dataA_j %>%
      left_join(dataY_j, by = c("id", "time")) %>%
      relocate(id, time, Y)
    
    # Define interventions
    interventions <- list(
      list(c(static, rep(0, min(jj, kappa + 1)))),
      list(c(static, rep(1, min(jj, kappa + 1))))
    )
    
    # Run g-formula
    g.results <- gformula(
      seed = boot_seed,
      obs_data = data,
      id = "id",
      time_name = "time",
      time_points = jj,
      outcome_name = "Y",
      outcome_type = "continuous_eof",
      ymodel = Y ~
        age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
        acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level + 
        hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group +
        opioid_pre_dis + nsaid_pre_dis +
        nsaid_ind * time + nsaid_ind * I(time^2) + lag1_nsaid_ind +
        healthcare_use + opioid_omeq + B01A_ind + hosp_days,
      
      basecovs = c(
        "age", "sex", "RHF", "trm_cent", "Kommuneindex", "pt_asa_preinjury",
        "acc_transport", "acc_fall", "acc_work", "inj_mechanism", "hosp_care_level",
        "hosp_icu_days", "hosp_los_days", "res_gos_dischg", "inntektGroup", 
        "ais_group", "opioid_pre_dis", "nsaid_pre_dis"
      ),
      
      covnames = c("nsaid_ind", "B01A_ind", "healthcare_use", "opioid_omeq", "hosp_days"),
      covtypes = c("binary", "binary", "normal", "normal", "normal"),
      
      histories = c(lagged),
      histvars = list(c("nsaid_ind", "B01A_ind", "healthcare_use", "opioid_omeq", "hosp_days")),
      
      covparams = list(covmodels = c(
        nsaid_ind ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
          acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level +
          hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group +
          opioid_pre_dis + nsaid_pre_dis + time + I(time^2) +
          healthcare_use + opioid_omeq + B01A_ind + hosp_days +
          lag1_nsaid_ind + lag2_nsaid_ind,
        
        B01A_ind ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
          acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level +
          hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group +
          opioid_pre_dis + nsaid_pre_dis + time + I(time^2) +
          lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days +
          lag1_nsaid_ind,
        
        healthcare_use ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
          acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level +
          hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group +
          opioid_pre_dis + nsaid_pre_dis + time + I(time^2) +
          lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days +
          lag1_nsaid_ind,
        
        opioid_omeq ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
          acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level +
          hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group +
          opioid_pre_dis + nsaid_pre_dis + time + I(time^2) +
          lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days +
          lag1_nsaid_ind,
        
        hosp_days ~ age + I(age^2) + sex + RHF + trm_cent + Kommuneindex + pt_asa_preinjury +
          acc_transport + acc_fall + acc_work + inj_mechanism + hosp_care_level +
          hosp_icu_days + hosp_los_days + res_gos_dischg + inntektGroup + ais_group +
          opioid_pre_dis + nsaid_pre_dis + time + I(time^2) +
          lag1_healthcare_use + lag1_opioid_omeq + lag1_B01A_ind + lag1_hosp_days +
          lag1_nsaid_ind
      )),
      
      intvars = list("nsaid_ind", "nsaid_ind"),
      interventions = interventions,
      int_descript = c("Add-on-0", "Add-on-1"),
      ref_int = 1,
      nsamples = 0
    )
    
    # Store results
    results_0_re[jj, bb + 2]    <- data.frame(g.results$result)[2, 4]
    results_1_re[jj, bb + 2]    <- data.frame(g.results$result)[3, 4]
    results_diff_re[jj, bb + 2] <- data.frame(g.results$result)[3, 6]
  }
  
  # Cumulative effect for this bootstrap replicate
  results_cum_re[bb] <- sum(results_diff_re[1:KK, bb + 2])
  
  # Optional: progress printout
  print(bb)
}

# Save bootstrap estimates
static_0_re    <- results_0_re
static_1_re    <- results_1_re
static_diff_re <- results_diff_re
static_cum_re  <- results_cum_re
