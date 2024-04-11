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
model <- list("linear", "linear", "linear", "linear")
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
> TS0$Longitudinal
[[1]]
[[1]]$Longitudinal_model
                     Est          SD       L_CI      U_CI
(Intercept) -0.003167629 0.051649161 -0.1016439 0.0865190
obstime      0.651195083 0.024117669  0.6057510 0.6996882
sigma2_e     0.304434795 0.006494876  0.2912420 0.3178109
sigma2_b     1.318007024 0.087255212  1.1564923 1.5036663


[[2]]
[[2]]$Longitudinal_model
                  Est         SD       L_CI      U_CI
(Intercept) 0.1079177 0.02551374 0.05766939 0.1556610
obstime     0.7318347 0.02821918 0.68070271 0.7879759
sigma2      0.1988282 0.00442706 0.19091613 0.2079220

[[2]]$Sigma
          Intercept      Time
Intercept 1.0611350 0.5255265
Time      0.5255265 1.0764285


[[3]]
[[3]]$Longitudinal_model
                    Est          SD        L_CI        U_CI
(Intercept) -0.01248891 0.020887352 -0.05146358  0.02923079
obstime      0.53419759 0.073747348  0.39215866  0.64446072
obstime2    -0.63622137 0.052401701 -0.70784160 -0.53403394
sigma2       0.19828364 0.004728854  0.18956590  0.20868304

[[3]]$Sigma
           Intercept       Time      Time2
Intercept  1.0512538  0.5376876 -0.1235974
Time       0.5376876  1.0998180 -0.2271203
Time2     -0.1235974 -0.2271203  0.4990175


> TS0$TDsurvival
$S_model
                Est         SD        L_CI        U_CI
x1       0.27251737 0.14241424 -0.01891624  0.54445881
x2       0.02010240 0.13661497 -0.25578120  0.27991697
Marker1 -0.03986755 0.05585649 -0.14893049  0.07103536
Marker2 -0.27965352 0.05428684 -0.38340941 -0.17559974
Marker3  0.21850350 0.05160108  0.11619111  0.32225479
h1       1.05664630 0.19742373  0.71824690  1.49755694
h2       0.90545060 0.16341794  0.61461051  1.24709989
h3       1.00960526 0.18004706  0.69618933  1.39670388
h4       0.97249676 0.17158193  0.67793851  1.33887908
h5       0.95898814 0.18192923  0.64538832  1.35723983
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
DP <- DP(object = TS0,
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
1     2 0.4172157
2     3 0.4085329
3     5 0.5036053
4    10 0.3905509
5    11 0.4036330
6    18 0.5323239
7    19 0.1724203
8    22 0.3148257
9    24 0.1963160
10   26 0.2776932
11   28 0.4662175
12   32 0.3599169
13   33 0.3803012
14   46 0.4364036
.
.
.
116 573 0.4714233
117 574 0.3506668
118 580 0.1595123
119 597 0.4564295
120 600 0.5671340

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
> Criteria(
+   s = 0.5, t = 0.5, Survt = dataSurv_v$survtime,
+   CR = dataSurv_v$death, P = DP$DP$est, cause = 1
+ )$Cri
```
with the following outputs:

```
          est        sd
AUC 0.5852190 0.1047990
BS  0.2308166 0.0254811
```


