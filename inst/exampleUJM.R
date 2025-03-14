\dontrun{
  library(survival)
data(dataLong)
data(dataSurv)

A0 <- UJM(
  formFixed = Y1 ~ obstime + x1 + x2, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, death) ~ x1+x2,
  dataLong, dataSurv, K = 15, model = "intercept", Obstime = "obstime",
  n.chains = 2, n.iter = 10, n.burnin = 5,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)


A1 <- UJM(
  formFixed = Y1 ~ obstime, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, death) ~ x1+x2,
  dataLong, dataSurv, K = 15, model = "intercept", Obstime = "obstime",
  n.chains = 2, n.iter = 10, n.burnin = 5,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)


A2 <- UJM(
  formFixed = Y1 ~ obstime + x1 + x2, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, death) ~ x1+x2,
  dataLong, dataSurv, K = 15, model = "linear", Obstime = "obstime",
  n.chains = 2, n.iter = 10, n.burnin = 5,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)


A3 <- UJM(
  formFixed = Y1 ~ obstime, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, death) ~ x1+x2,
  dataLong, dataSurv, K = 15, model = "linear", Obstime = "obstime",
  n.chains = 2, n.iter = 10, n.burnin = 5,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)




A4 <- UJM(
  formFixed = Y1 ~ obstime + x1 + x2, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, death) ~ x1+x2,
  dataLong, dataSurv, K = 15, model = "quadratic", Obstime = "obstime",
  n.chains = 2, n.iter = 10, n.burnin = 5,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)


A5 <- UJM(
  formFixed = Y1 ~ obstime, formRandom = ~obstime,
  formGroup = ~id, formSurv = Surv(survtime, death) ~ x1+x2,
  dataLong, dataSurv, K = 15, model = "quadratic", Obstime = "obstime",
  n.chains = 2, n.iter = 10, n.burnin = 5,
  n.thin = 1,
  DIC = TRUE, quiet = FALSE
)


}
