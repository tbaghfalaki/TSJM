For computing risk prediction, it is necessary to first estimate the parameters of the model and then utilize its object in the functions for dynamic prediction.


# Parameter estimation


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

We are considering four markers; therefore, we require four lists as follows: one for the fixed effects model, one for the random effects model, and another for the survival model.

```
formFixed <- list(Y1 ~ obstime, Y2 ~ obstime, Y3 ~ obstime, Y4 ~ obstime)
formRandom <- list(~obstime, ~obstime, ~obstime, ~obstime)
formGroup <- list(~id, ~id, ~id, ~id)
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
model <- list("linear", "linear", "linear", "linear")
```
Finally, we utilize the TS function, which has been applied to the training data:

```
TS1 <- TS(formFixed, formRandom, formGroup, formSurv,
         nmark = 4, K1 = 15, K2 = 15,
         model = model, n.chains1 = 1, n.iter1 = 2000, n.burnin1 = 1000,
         n.thin1 = 1,  n.chains2 = 1, n.iter2 = 3000, n.burnin2 = 1000,
         n.thin2 = 1, Obstime = "obstime", ncl = 2,
         DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
)
```

The outputs of this function is as follows: 

```
> TS1$Longitudinal
[[1]]
[[1]]$Longitudinal_model
                  Est          SD       L_CI      U_CI
(Intercept) 0.1186460 0.027971557 0.05358023 0.1792465
obstime     0.6404244 0.036004289 0.56167250 0.7022207
sigma2      0.2009000 0.004351306 0.19202376 0.2093111

[[1]]$Sigma
          Intercept      Time
Intercept  1.073042 0.4470210
Time       0.447021 0.8869411


[[2]]
[[2]]$Longitudinal_model
                   Est          SD       L_CI      U_CI
(Intercept) 0.06243818 0.020493833 0.02663963 0.1049948
obstime     0.51864278 0.044125250 0.42139482 0.5939914
sigma2      0.19970782 0.004462027 0.19116195 0.2087919

[[2]]$Sigma
          Intercept      Time
Intercept 1.0535532 0.5454659
Time      0.5454659 1.0804504


[[3]]
[[3]]$Longitudinal_model
                    Est          SD         L_CI        U_CI
(Intercept) -0.08047871 0.024246694 -0.120527757 -0.02565325
obstime      0.08538702 0.049140291  0.007382119  0.18494356
sigma2       0.19739551 0.004319167  0.189187016  0.20599431

[[3]]$Sigma
          Intercept      Time
Intercept 1.0911398 0.4644015
Time      0.4644015 1.1014568


[[4]]
[[4]]$Longitudinal_model
                    Est          SD       L_CI         U_CI
(Intercept) -0.06499082 0.029161950 -0.1171062 -0.005371892
obstime      0.29497339 0.025756866  0.2449191  0.354918011
sigma2       0.19357406 0.004277406  0.1857173  0.202200815

[[4]]$Sigma
          Intercept      Time
Intercept 1.1100642 0.6453169
Time      0.6453169 1.1460680


> TS1$TDsurvival
$S_model
                Est         SD        L_CI         U_CI
x1       0.29541449 0.14038668  0.03089244  0.578761990
x2       0.02109871 0.12670912 -0.23382368  0.264229926
Marker1 -0.11169877 0.05780525 -0.22983375 -0.001480224
Marker2 -0.31113639 0.05081028 -0.41252658 -0.205722528
Marker3  0.18669627 0.05282661  0.08747509  0.294427578
Marker4  0.14789465 0.05159023  0.04769651  0.251546871
h1       1.02437491 0.18043331  0.71654485  1.432232923
h2       0.87988795 0.15838168  0.59846500  1.219036324
h3       0.97621759 0.17490050  0.66560180  1.341800893
h4       0.93385031 0.16009532  0.65304219  1.278264073
h5       1.00604702 0.18025980  0.69027367  1.389592704
```

If we consider n.chains1 or n.chains2 > 1, the values of the Gelman-Rubin criteria are also provided, which helps in checking the convergence of the MCMC.



# Dynamic prediction 

For dynamic prediction, we offer two functions: DP and DP_CI. The distinction between these two functions lies in their computing approach, as described in the main paper. The first function employs a first-order approximation of dynamic prediction, while the second utilizes Monte Carlo approximation, allowing for the computation of credible intervals.

#### DP function
The main arguments of this functions are as follows:
- object an object inheriting from class TS
-  s the landmark time for prediction
-  t the window of prediction for prediction
-  dataLong data set of observed longitudinal variables (validation set).
-  dataSurv data set of observed survival variables (validation set).

```
DP <- DP(TS1,
         s = 0.5, t = 0.5, n.chains = 1, n.iter = 2000, n.burnin = 1000,
         n.thin = 1,
         DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
)
```

The outputs of this function is as follows: 

```
> DP
$DP
     id       est
1     2 0.2976367
2     3 0.3391111
3     5 0.5425427
4    10 0.4484265
5    11 0.3715304
6    18 0.4328885
7    19 0.1572549
8    22 0.2864214
9    24 0.2283791
10   26 0.2647648
11   28 0.4001724
.
.
.
114 540 0.3755988
115 546 0.4698084
116 573 0.5565222
117 574 0.3375089
118 580 0.1553600
119 597 0.6882316
120 600 0.6806330

$s
[1] 0.5

$t
[1] 0.5
```
# Computing AUC and BS for the predictions
For this purpose, we use DPCri package <https://github.com/tbaghfalaki/DPCri>.

Computing the criteria using this package is straightforward, as demonstrated by the following commands:

- s the landmark time for prediction
- t the window of prediction for prediction
- Survt the survival time
- CR the indicator for competing risks or censoring
- P the risk predictions
- cause the main cause for prediction


Consider the following command: 

```
Criteria(
  s = 0.5, t = 0.5, Survt = dataSurv_v$survtime,
  CR = dataSurv_v$death, P = DP$DP$est, cause = 1
)$Cri
```
with the following outputs:

```
       est         sd
AUC 0.6488860 0.10095338
BS  0.2177984 0.02648602
```


