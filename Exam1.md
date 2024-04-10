Getting Started
---------------

```
library(TSJM)
```
Generating summary statistics


Loading data of the package:

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
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}x_1+b_{0ki}+\varepsilon_{ikt}$


linear:
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}x_1+b_{0ki}+b_{1ki} t+\varepsilon_{ikt}$


quadratic:
$Y_{ik}(t)= \beta_{0k}+\beta_{1k}t+\beta_{2k}t^2+\beta_{3k}x_1+b_{0ki}+b_{1ki} t+b_{1ki} t^2+\varepsilon_{ikt}$


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
- ncl the number of nodes to be forked for parallel computing
- n.chains1 the number of parallel chains for the model in the first stage; default is 1.
- n.iter1 integer specifying the total number of iterations in the first stage; default is 1000.
- n.burnin1 integer specifying how many of n.iter to discard as burn-in in the first stage; default is 5000.
- n.thin1 integer specifying the thinning of the chains in the first stage; default is 1.
- n.chains2 the number of parallel chains for the model in the second stage; default is 1.
- n.iter2 integer specifying the total number of iterations in the second stage; default is 1000.
- n.burnin2 integer specifying how many of n.iter to discard as burn-in in the second stage; default is 5000.
- n.thin2 integer specifying the thinning of the chains in the second stage; default is 1.
- simplify Logical; the option for simplifying the use of CS and DS; default is TRUE.
- DIC Logical; if TRUE (default), compute deviance, pD, and DIC. The rule pD=var(deviance) / 2 is used.
- quiet Logical, whether to suppress stdout in jags.model().
-----------------

As an example, consider the following command, where this implementation has been performed on training data:


```
TS0 <- TS(formFixed, formRandom, formGroup, formSurv,
         nmark = 3, K1 = 15, K2 = 15,
         model = model, n.chains1 = 1, n.iter1 = 2000, n.burnin1 = 1000,
         n.thin1 = 1,  n.chains2 = 1, n.iter2 = 3000, n.burnin2 = 1000,
         n.thin2 = 1, Obstime = "obstime", ncl = 3,
         DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
)
```

The outputs of this function is as follows: 

