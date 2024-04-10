\dontrun{
rm(list=ls())
data(dataLong)
data(dataSurv)
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



formFixed <- list(Y1 ~ obstime, Y2 ~ obstime, Y3 ~ obstime)
formRandom <- list(~obstime, ~obstime, ~obstime)
formGroup <- list(~id, ~id, ~id)
formSurv <- survival::Surv(survtime, death) ~ x1 + x2
model <- list("intercept", "linear", "quadratic")

TS0 <- TS(formFixed, formRandom, formGroup, formSurv,
         nmark = 3, K1 = 15, K2 = 15,
         model = model, n.chains1 = 1, n.iter1 = 20, n.burnin1 = 10,
         n.thin1 = 1,  n.chains2 = 1, n.iter2 = 30, n.burnin2 = 10,
         n.thin2 = 1, simplify=TRUE, Obstime = "obstime", ncl = 2,
         DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
)

model <- list("linear", "linear", "linear")

TS1 <- TS(formFixed, formRandom, formGroup, formSurv,
         nmark = 3, K1 = 15, K2 = 15,
         model = model, n.chains1 = 1, n.iter1 = 20, n.burnin1 = 10,
         n.thin1 = 1,  n.chains2 = 1, n.iter2 = 30, n.burnin2 = 10,
         n.thin2 = 1, simplify=TRUE, Obstime = "obstime", ncl = 2,
         DIC = TRUE, quiet = FALSE, dataLong_t, dataSurv_t
)



DP <- DP(TS1,
         s = 0.5, t = 0.5, n.chains = 1, n.iter = 20, n.burnin = 10,
         n.thin = 1,
         DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
)



DP <- DP_CI(TS0,
         s = 0.5, t = 0.5, n.chains = 1, n.iter = 20, n.burnin = 10,
         n.thin = 1, mi=3,
         DIC = TRUE, quiet = FALSE, dataLong = dataLong_v, dataSurv = dataSurv_v
)


}




