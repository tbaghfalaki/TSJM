\dontrun{
# Remove all objects from the workspace to start with a clean environment
rm(list = ls())

# Load the data
data(dataLong)
data(dataSurv)

# Set the seed for reproducibility
set.seed(2)

# Create training and validation indices
INDTRAIN <- sample(dataSurv$id, 0.8 * (dim(dataSurv)[1]))
INDVALID <- dataSurv$id[-INDTRAIN]

# Subset data into training and validation sets based on indices
dataLong_t <- subset(dataLong, dataLong$id %in% INDTRAIN)
dataSurv_t <- subset(dataSurv, dataSurv$id %in% INDTRAIN)
dataLong_v <- subset(dataLong, dataLong$id %in% INDVALID)
dataSurv_v <- subset(dataSurv, dataSurv$id %in% INDVALID)

# Define the model formulas for fixed effects, random effects, group effects,
# and survival model
formFixed <- list(Y1 ~ obstime, Y2 ~ obstime, Y3 ~ obstime)
formRandom <- list(~obstime, ~obstime, ~obstime)
formGroup <- list(~id, ~id, ~id)
formSurv <- survival::Surv(survtime, death) ~ x1 + x2

# Specify the model types
model <- list("intercept", "linear", "quadratic")

# Fit the TSC (Time-Scale Continuous) model
TSC0 <- TSC(formFixed, formRandom, formGroup, formSurv,
  nmark = 3, K1 = 15, K2 = 15,
  model = model, n.chains1 = 1, n.iter1 = 500, n.burnin1 = 200,
  n.thin1 = 1,
  Obstime = "obstime", ncl = 2, Limp = 50,
  DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
)

##################################################################

# Merge survival data for modeling
surv_data <- survival::tmerge(dataSurv_t, dataSurv_t,
  id = id,
  endpt = event(survtime, death)
)

# Merge longitudinal data with survival data
long.data1 <- survival::tmerge(surv_data, dataLong_t,
  id = id, lY1 = tdc(obstime, TSC0$lPredY[, 1]),
  lY2 = tdc(obstime, TSC0$lPredY[, 2]),
  lY3 = tdc(obstime, TSC0$lPredY[, 3])
)

# Fit a Cox proportional hazards model using the joint model's longitudinal predictions
cox_model_tsjm <- coxph(
  Surv(time = tstart, time2 = tstop, endpt) ~ x1 + x2 +
    lY1 + lY2 + lY3,
  data = long.data1, id = id
)

####################################################

# Prepare data for multiple imputation
surv_data_new <- survival::tmerge(dataSurv_t, dataSurv_t,
  id = id,
  endpt = event(survtime, death)
)

Limp <- 50
SD2 <- matrix(0, Limp, length(cox_model_tsjm$coefficients))
TEM <- matrix(0, Limp, length(cox_model_tsjm$coefficients))

# Perform multiple imputation to estimate variability of coefficients
for (k in 1:Limp) {
  lPredY <- TSC0$LPredY[[k]]
  long.data1 <- survival::tmerge(surv_data_new, dataLong_t,
    id = id,
    lY1 = tdc(obstime, lPredY[, 1]),
    lY2 = tdc(obstime, lPredY[, 2]),
    lY3 = tdc(obstime, lPredY[, 3])
  )

  cox_model <- coxph(
    Surv(time = tstart, time2 = tstop, endpt) ~ x1 + x2 +
      lY1 + lY2 + lY3,
    data = long.data1, id = id
  )

  TEM[k, ] <- summary(cox_model)$coefficients[, 1]
  SD2[k, ] <- summary(cox_model)$coefficients[, 3]
}

# Compute the within-imputation variance (Wv) and between-imputation variance (Bv)
Wv <- apply(SD2^2, 2, mean)
Bv <- apply(TEM, 2, var)
sdnew <- sqrt(Wv + (Limp + 1) / Limp * Bv)

# Combine results into a single summary table
Res <- cbind(
  cox_model_tsjm$coefficients, sdnew,
  cox_model_tsjm$coefficients - qnorm(.95) * sdnew,
  cox_model_tsjm$coefficients + qnorm(.95) * sdnew
)

colnames(Res) <- c("coefficients", "sd", "L_CI", "U_CI")
print(Res)

################## Dynamic Prediction (DP) ##############

# Perform dynamic predictions using the fitted TSC model
LP_v0 <- LP_v(TSC0,
  s = 0.5, t = 0.5, n.chains = 1, n.iter = 2000, n.burnin = 1000,
  n.thin = 1, DIC = TRUE, quiet = FALSE, dataLong = dataLong_v,
  dataSurv = dataSurv_v
)

# Time points for prediction
s <- 0.5
Dt <- 0.5

# Load necessary libraries
library(dplyr)
library(DPCri)

# Filter longitudinal data up to time point s
data_Long_s <- dataLong_v[dataLong_v$obstime <= s, ]

# Add longitudinal predictions to the validation data
long.data_v <- data_Long_s %>%
  mutate(
    lY1 = LP_v0$lPredY[, 1],
    lY2 = LP_v0$lPredY[, 2],
    lY3 = LP_v0$lPredY[, 3]
  )

# Prepare survival data for validation
temp_v <- subset(dataSurv_v, select = c(id, survtime, death))
pbc21_v <- tmerge(temp_v, temp_v, id = id, endpt = event(survtime, death))
long.data1_v <- tmerge(pbc21_v, long.data_v,
  id = id, lY1 = tdc(obstime, lY1),
  lY2 = tdc(obstime, lY2), lY3 = tdc(obstime, lY3)
)

# Extract the baseline hazard function from the Cox model
base_hazard <- basehaz(cox_model_tsjm, centered = FALSE)

# Estimate cumulative hazard at time s for each individual
cum_hazard_s <- rep(0, nrow(dataSurv_v))
for (i in 1:nrow(dataSurv_v)) {
  individual_data <-
    long.data1_v[long.data1_v$id == dataSurv_v$id[i] & long.data1_v$tstop <= s, ]
  X1 <- c(dataSurv_v$x1[i], dataSurv_v$x2[i])
  if (nrow(individual_data) > 0) {
    linear_predictor <- sum(coef(cox_model_tsjm) *
      c((X1), unlist(tail(individual_data[, c("lY1", "lY2", "lY3")], 1))))
    baseline_hazard_s <- base_hazard$hazard[base_hazard$time <= s]
    cum_hazard_s[i] <- tail(baseline_hazard_s, 1) * exp(linear_predictor)
  }
}

# Estimate survival probability at time s
surv_prob_s <- exp(-cum_hazard_s)

# Estimate cumulative hazard at time s + Dt for each individual
cum_hazard_sDt <- rep(0, nrow(dataSurv_v))
for (i in 1:nrow(dataSurv_v)) {
  individual_data <-
    long.data1_v[long.data1_v$id == dataSurv_v$id[i] & long.data1_v$tstop <= (s + Dt), ]
  X1 <- c(dataSurv_v$x1[i], dataSurv_v$x2[i])
  if (nrow(individual_data) > 0) {
    linear_predictor <- sum(coef(cox_model_tsjm) * c(
      (X1),
      unlist(tail(individual_data[, c("lY1", "lY2", "lY3")], 1))
    ))
    baseline_hazard_sDt <- base_hazard$hazard[base_hazard$time <= (s + Dt)]
    cum_hazard_sDt[i] <- tail(baseline_hazard_sDt, 1) * exp(linear_predictor)
  }
}

# Estimate survival probability at time s + Dt
surv_prob_sDt <- exp(-cum_hazard_sDt)

# Calculate dynamic prediction
DP <- 1 - surv_prob_sDt / surv_prob_s

# Compute the cumulative risk index (CRI) using the dynamic prediction
CRI_TSJM <- Criteria(
  s = s, t = Dt, Survt = dataSurv_v$survtime,
  CR = dataSurv_v$death, P = DP, cause = 1
)$Cri[, 1]

# Print the cumulative risk index
print(CRI_TSJM)
}
