# Estimating Add-On Effects

This R code is part of the supplementary material to the paper “Add-On Regimes and Their Relevance for Quantifying the Effects of Opioid-Sparing Treatments.”

The code was used to estimate opioid-sparing effects of supplementing opioid treatments with nonsteroidal anti-inflammatory drugs (NSAIDs), in a cohort of Norwegian trauma patients, using national registry data. 

More generally, the code can be used to estimate add-on effects—defined as the average causal effect of implementing one intervention (e.g., NSAID treatment) whenever another event occurs (e.g., opioid treatment) during a treatment period—on an outcome of interest (e.g., subsequent opioid use).

This folder contains two `R` scripts:
- `gformula.R`: Implements estimation of the add-on effect using a g-formula estimator. The results presented in Section 7 of the main article are based on this script.
- `IPW.R`: Implements estimation of the add-on effect using an inverse probability weighting (IPW) estimator.

## Data

The input data `data_pre` is in long format, where each row corresponds to one time point for one study subject. Below is a description of the key variables required:
- `id`: A unique identifier for each study subject.
- `time`: Discrete follow-up time, starting from 0. Each subject has one row per time point.
- `nsaid_ind`: A binary indicator of NSAID dispensing at each time point. This is the treatment indicator in the present example.
- `opioid_omeq`: The outcome variable in this example, representing opioid use measured in oral morphine equivalents. This variable is also used in defining the add-on regime (see note below).
- `death_ind`: A binary indicator equal to 1 if the subject has died by time k, and 0 otherwise.

Additional covariates: The dataset also includes a range of baseline and time-varying covariates used for adjustment in the estimation procedures. See Appendix C of the paper for a detailed description of all variables used in the analysis.

Note on censoring: In the data used in this specific example, there is no censoring during follow-up from month 0 to month 21. 

Note on `opioid_omeq`: In this specific example on opioid-sparing effects, the variable `opioid_omeq` play a dual role: it serves both as the outcome variable (since we are interested in the effect on opioid dose) and as a key time-varying covariate (since we are interested in the effect of adding NSAIDs when opioids are administered). In other applications, these roles may be represented by distinct variables. For example, when evaluating whether prescribing probiotics alongside antibiotics reduces antibiotic-associated diarrhea, a variable `diarrhea` may be the outcome variable of interest, indicating diarrhea at time $k$, while a separate time-varying variable, `antibiotic`, indicates antibiotic use at time $k$.

The input data `data_pre` stem from observational data from [NTRplus]([https://pages.github.com/](https://www.ous-research.no/home/ipot/Projects/20448)), which links the [Norwegian National Trauma Registry](https://www.ous-research.no/home/ipot/Projects/20448)
to several national databases: the [Norwegian Prescription Database](https://www.norpd.no/), the [Cause of Death Registry](https://www.fhi.no/en/ch/cause-of-death-registry/), the [Norwegian Patient Registry](https://helsedata.no/en/forvaltere/norwegian-institute-of-public-health/norwegian-patient-registry-npr/), and [Statistics Norway](https://www.ssb.no/en). This data is not publicly available. 


