rm(list=ls())
data(dataLong)
data(dataSurv)
set.seed(2)
INDTRAIN <- sample(dataSurv$id, 0.5 * (dim(dataSurv)[1]))
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



formFixed <- list(Y1 ~ obstime, Y2 ~ obstime, Y3 ~ obstime)
formRandom <- list(~obstime, ~obstime, ~obstime)
formGroup <- list(~id, ~id, ~id)
formSurv <- survival::Surv(survtime, death) ~ x1 + x2
model <- list("intercept", "linear", "quadratic")



TS <- TS(formFixed, formRandom, formGroup, formSurv,
         nmark = 3, K1 = 15, K2 = 15,
         model = model, n.chains1 = 1, n.iter1 = 2000, n.burnin1 = 1000,
         n.thin1 = 1,  n.chains2 = 1, n.iter2 = 3000, n.burnin2 = 1000,
         n.thin2 = 1, simplify=TRUE, Obstime = "obstime", ncl = 3,
         DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
)








