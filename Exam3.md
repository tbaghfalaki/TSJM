Getting Started
---------------

```
library(TSJM)
```
Loading the data from the package includes both longitudinal data in long format and survival data. It's essential to ensure that the same subject (ID) is present in both datasets.

```
data(dataLong)
data(dataSurv)
```

Dividing data to 80% training data and 20% validation set


```
set.seed(2)
INDTRAIN <- sample(dataSurv$id, 0.8 * (dim(dataSurv)[1]))
INDVALID <- dataSurv$id[-INDTRAIN]
dataLong_t <- subset(
  dataLong,
  dataLong$id %in% INDTRAIN
)
dataSurv_t <- subset(
  dataSurv,
  dataSurv$id %in% INDTRAIN
)

dataLong_v <- subset(
  dataLong,
  dataLong$id %in% INDVALID
)
dataSurv_v <- subset(
  dataSurv,
  dataSurv$id %in% INDVALID
)
```

We are considering three markers; therefore, we require three lists as follows: one for the fixed effects model, one for the random effects model, and another for the survival model.

```
formFixed <- list(Y1 ~ obstime, Y2 ~ obstime, Y3 ~ obstime)
formRandom <- list(~obstime, ~obstime, ~obstime)
formGroup <- list(~id, ~id, ~id)
formSurv <- survival::Surv(survtime, death) ~ x1 + x2
```

We need to choose the model for the marker trend among "intercept," "linear," and "quadratic." For instance, if we consider a covariate $x_1$, the options are as follows:

intercept:
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}x_{1i}+b_{0ki}+\varepsilon_{ikt}$


linear:
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}x_{1i}+b_{0ki}+b_{1ki} t+\varepsilon_{ikt}$


quadratic:
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}t^2+\beta_{3k}x_{1i}+b_{0ki}+b_{1ki} t+b_{1ki} t^2+\varepsilon_{ikt}$


This has been done by considering the following command:

```
model <- list("intercept", "linear", "quadratic")
```
Finally, we have to use the TS function with the following arguments:


- formFixed a list of formulas for fixed part of longitudinal model
- formRandom a list of formulas for random part of longitudinal model
- formGroup a list of formulas specifying the cluster variable for Y (e.g. = list (~ subject, ~ subject,...))
- formSurv formula for survival model
- dataLong data set of observed longitudinal variables.
- dataSurv data set of observed survival variables.
- nmark the number of longitudinal markers
- K1 Number of nodes and weights for calculating Gaussian quadrature in the first stage.
- K2 Number of nodes and weights for calculating Gaussian quadrature in the second stage.
- model a list of the models for the longitudinal part which includes "linear" or "quadratic".
- Obstime the observed time in longitudinal data
- Limp the number of multiple imputation; default is 10.
- ncl the number of nodes to be forked for parallel computing
- n.chains1 the number of parallel chains for the model in the first stage; default is 1.
- n.iter1 integer specifying the total number of iterations in the first stage; default is 1000.
- n.burnin1 integer specifying how many of n.iter to discard as burn-in in the first stage; default is 5000.
- n.thin1 integer specifying the thinning of the chains in the first stage; default is 1.
- simplify Logical; the option for simplifying the use of CS and DS; default is TRUE.
- DIC Logical; if TRUE (default), compute deviance, pD, and DIC. The rule pD=var(deviance) / 2 is used.
- quiet Logical, whether to suppress stdout in jags.model().
-----------------

As an example, consider the following command, where this implementation has been performed on training data:


```
TSC0 <- TSC(formFixed, formRandom, formGroup, formSurv,
            nmark = 3, K1 = 15, K2 = 15,
            model = model, n.chains1 = 1, n.iter1 = 500, n.burnin1 = 200, n.thin1 = 1,
            Obstime = "obstime", ncl = 2, Limp = 50,
            DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
)
```

The outputs of this function is as follows: 

```
> TSC0$Longitudinal
[[1]]
[[1]]$Longitudinal_model
                    Est          SD        L_CI      U_CI
(Intercept) -0.01649726 0.044819289 -0.08197063 0.1010317
obstime      0.65105117 0.022232601  0.60667648 0.6888478
sigma2_e     0.30375675 0.006596923  0.29182300 0.3175708
sigma2_b     1.31253403 0.087325790  1.14812838 1.4978071


[[2]]
[[2]]$Longitudinal_model
                   Est          SD       L_CI      U_CI
(Intercept) 0.09573039 0.020999828 0.04830476 0.1289396
obstime     0.59711720 0.022313466 0.56030509 0.6466764
sigma2      0.19875196 0.004595108 0.19054850 0.2080349

[[2]]$Sigma
          Intercept      Time
Intercept 1.0671703 0.4713698
Time      0.4713698 1.0562758


[[3]]
[[3]]$Longitudinal_model
                    Est         SD        L_CI        U_CI
(Intercept)  0.01766145 0.02254862 -0.02003771  0.06819337
obstime      0.52638949 0.06454787  0.38411543  0.65078173
obstime2    -0.68826213 0.03193343 -0.73951093 -0.60657123
sigma2       0.20031963 0.00463674  0.19174606  0.20931106

[[3]]$Sigma
           Intercept       Time      Time2
Intercept  1.0500927  0.5486612 -0.2229174
Time       0.5486612  1.4515064 -0.4888390
Time2     -0.2229174 -0.4888390  0.6345375
```


If we consider n.chains1 or n.chains2 > 1, the values of the Gelman-Rubin criteria are also provided, which helps in checking the convergence of the MCMC.

So far, only the first stage has been completed. The second stage of our approach will be conducted using the following code, leveraging the outputs from TSC0 by running a proportional hazards model with time-dependent covariates:

```
# Merge survival data for modeling
surv_data <- survival::tmerge(dataSurv_t, dataSurv_t, id = id, endpt = event(survtime, death))

# Merge longitudinal data with survival data
long.data1 <- survival::tmerge(surv_data, dataLong_t,
                               id = id, lY1 = tdc(obstime, TSC0$lPredY[, 1]),
                               lY2 = tdc(obstime, TSC0$lPredY[, 2]), lY3 = tdc(obstime, TSC0$lPredY[, 3])
)

# Fit a Cox proportional hazards model using the joint model's longitudinal predictions
cox_model_tsjm <- coxph(Surv(time = tstart, time2 = tstop, endpt) ~ lY1 + lY2 + lY3,
                        data = long.data1, id = id
)

```
