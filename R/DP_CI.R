#'  Dynamic prediction with credible interval
#'
#' @description
#' Dynamic prediction for TSJM with CI
#'
#'
#' @details
#' Estimate DP for joint modeling based on TS approach.
#'
#' @param object an object inheriting from class TS
#' @param dataLong data set of observed longitudinal variables.
#' @param dataSurv data set of observed survival variables.
#' @param s the landmark time for prediction
#' @param t the window of prediction for prediction
#' @param mi the number of multiple imputation for Monte-Carlo approximation; default is 10.
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 1000.
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
#' @param n.thin integer specifying the thinning of the chains; default is 1.
#' @param DIC Logical; if TRUE (default), compute deviance, pD, and DIC. The rule pD=var(deviance) / 2 is used.
#' @param quiet Logical, whether to suppress stdout in jags.model().
#'
#'
#'
#' @importFrom stats quantile rnorm model.frame model.matrix
#'
#' @return
#' - mu.vect list of posterior mean for each parameter
#' - sd.vect list of standard error for each parameter
#' - 2.5% list of posterior mode for each parameter
#' - 97.5% list of posterior median for each parameter
#' - Rhat Gelman and Rubin diagnostic for all parameter
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}
#'
#' @example inst/exampleTS.R
#'
#' @md
#' @export

DP_CI <- function(object, s = s, t = t, mi = mi, n.chains = n.chains, n.iter = n.iter, n.burnin = floor(n.iter / 2),
                  n.thin = max(1, floor((n.iter - n.burnin) / 1000)),
                  DIC = TRUE, quiet = FALSE, dataLong, dataSurv) {
  Dt <- t
  KK <- 1000000

  formFixed <- object$formFixed
  formRandom <- object$formRandom
  formGroup <- object$formGroup
  formSurv <- object$formSurv
  model <- object$model
  Obstime <- object$Obstime
  nmark <- object$nmark
  mu1 <- object$mu1
  peice <- object$peice
  #######


  ########### univariate_jm_random_effect_estimation

  model_fileI1b <- "model{
  for(i in 1:n){
    #Longitudinal observations
    Y1[i]~dnorm(mu1[i],tau1)
    mu1[i]<-inprod(betaL[],X[i,])+b[id[i]]
  }
  #Survival and censoring times
  #Hazard function

  for(k in 1:n2){
  for(j in 1:K){
    # Scaling Gauss-Kronrod/Legendre quadrature
    xk11[k,j]<-(xk[j]+1)/2*Time[k]
    wk11[k,j]<- wk[j]*Time[k]/2
  }
    linearpred[k]<-betaL[1]+b[k]+(betaL[indtime])*Time[k]

    Alpha0[k]<- inprod(betaS[1:NbetasS],XS[k,1:NbetasS])+
                 gamma1*(betaL[1]+b[k])


    Alpha1[k]<- gamma1*(betaL[indtime])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
                (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
                (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
                (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
                (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature

      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
                      (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
                      (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
                      (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
                      (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])
  phi[k]<- KK-Death[k]*log(haz[k])-logSurv[k]
    #Definition of the survival log-likelihood using zeros trick
    zeros[k]~dpois(phi[k])
    #Random effects
    b[k]~dnorm(0,Omega)

  }
Omega<-1/Sigma

  tau1<-1/sigma1
}"

  model_fileIb <- "model{
  for(i in 1:n){
    #Longitudinal observations
    Y1[i]~dnorm(mu1[i],tau1)
    mu1[i]<-inprod(betaL[],X[i,])+b[id[i]]
  }
  #Survival and censoring times
  #Hazard function
  for(k in 1:n2){
  for(j in 1:K){
    # Scaling Gauss-Kronrod/Legendre quadrature
    xk11[k,j]<-(xk[j]+1)/2*Time[k]
    wk11[k,j]<- wk[j]*Time[k]/2
  }
    linearpred[k]<-inprod(betaL[nindtime],Xv[k,])+b[k]+(betaL[indtime])*Time[k]

    Alpha0[k]<- inprod(betaS[1:NbetasS],XS[k,1:NbetasS])+
                 gamma1*(inprod(betaL[nindtime],Xv[k,])+b[k])


    Alpha1[k]<- gamma1*(betaL[indtime])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
                (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
                (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
                (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
                (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature

      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
                      (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
                      (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
                      (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
                      (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])
    phi[k]<- KK-Death[k]*log(haz[k])-logSurv[k]


    #Definition of the survival log-likelihood using zeros trick
    zeros[k]~dpois(phi[k])
    #Random effects
    b[k]~dnorm(0,Omega)

  }

Omega<-1/Sigma

  tau1<-1/sigma1

}"




  model_fileL1b <- "model{
  for(i in 1:n){
    #Longitudinal observations
    Y1[i]~dnorm(mu1[i],tau1)
    mu1[i]<-inprod(betaL[],X[i,])+inprod(b[id[i],],Z[i,])
  }
  #Survival and censoring times
  #Hazard function
  for(k in 1:n2){
  for(j in 1:K){
    # Scaling Gauss-Kronrod/Legendre quadrature
    xk11[k,j]<-(xk[j]+1)/2*Time[k]
    wk11[k,j]<- wk[j]*Time[k]/2
  }}
  for(k in 1:n2){
    linearpred[k]<-betaL[1]+b[k,1]+(betaL[indtime]+b[k,2])*Time[k]

    Alpha0[k]<- inprod(betaS[1:NbetasS],XS[k,1:NbetasS])+
                 gamma1*(betaL[1]+b[k,1])


    Alpha1[k]<- gamma1*(betaL[indtime]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
                (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
                (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
                (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
                (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature

      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
                      (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
                      (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
                      (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
                      (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])
    phi[k]<- KK-Death[k]*log(haz[k])-logSurv[k]


    #Definition of the survival log-likelihood using zeros trick
    zeros[k]~dpois(phi[k])


    #Random effects
    b[k,1:Nb]~dmnorm(mub[],Omega[,])

  }
 Omega[1:Nb,1:Nb]<-inverse(Sigma[,])
  tau1<-1/sigma1
}"

  model_fileLb <- "model{
  for(i in 1:n){
    #Longitudinal observations
    Y1[i]~dnorm(mu1[i],tau1)
    mu1[i]<-inprod(betaL[],X[i,])+inprod(b[id[i],],Z[i,])
  }
  #Survival and censoring times
  #Hazard function
  for(k in 1:n2){
  for(j in 1:K){
    # Scaling Gauss-Kronrod/Legendre quadrature
    xk11[k,j]<-(xk[j]+1)/2*Time[k]
    wk11[k,j]<- wk[j]*Time[k]/2
  }}
  for(k in 1:n2){
    linearpred[k]<-inprod(betaL[nindtime],Xv[k,])+b[k,1]+(betaL[indtime]+b[k,2])*Time[k]
    Alpha0[k]<- inprod(betaS[1:NbetasS],XS[k,1:NbetasS])+
                 gamma1*(inprod(betaL[nindtime],Xv[k,])+b[k,1])


    Alpha1[k]<- gamma1*(betaL[indtime]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
                (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
                (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
                (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
                (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature

      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
                      (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
                      (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
                      (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
                      (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])
   phi[k]<- KK-Death[k]*log(haz[k])-logSurv[k]


    #Definition of the survival log-likelihood using zeros trick
    zeros[k]~dpois(phi[k])

    #Random effects
    b[k,1:Nb]~dmnorm(mub[],Omega[,])

  }
 Omega[1:Nb,1:Nb]<-inverse(Sigma[,])
  tau1<-1/sigma1

}"




  model_fileQ1b <- "model{
  for(i in 1:n){
    #Longitudinal observations
    Y1[i]~dnorm(mu1[i],tau1)
    mu1[i]<-inprod(betaL[],X[i,])+inprod(b[id[i],],Z[i,])
  }
  for(k in 1:n2){
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
    }}
  #Survival and censoring times
  #Hazard function
  for(k in 1:n2){
    linearpred[k]<-betaL[1]+b[k,1]+(betaL[indtime[1]]+b[k,2])*Time[k]+
      (betaL[indtime[2]]+b[k,3])*pow(Time[k],2)


    Alpha0[k]<-  inprod(betaS[1:NbetasS],XS[k,1:NbetasS])+
                 gamma1*(betaL[1]+b[k,1])
    Alpha1[k]<- gamma1*(betaL[indtime[1]]+b[k,2])
    Alpha2[k]<- gamma1*(betaL[indtime[2]]+b[k,3])

    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
                  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
                  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
                  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
                  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k]+Alpha2[k]*pow(Time[k],2))


    for(j in 1:K){
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
                        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
                        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
                        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
                        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j]+Alpha2[k]*pow(xk11[k,j],2))

        }

    logSurv[k]<- -inprod(wk11[k,],chaz[k,])
    phi[k]<- KK-Death[k]*log(haz[k])-logSurv[k]


    #Definition of the survival log-likelihood using zeros trick
    zeros[k]~dpois(phi[k])
    #Random effects
    b[k,1:Nb]~dmnorm(mub[],Omega[,])

  }
   Omega[1:Nb,1:Nb]<-inverse(Sigma[,])
  tau1<-1/sigma1

}"

  model_fileQb <- "model{
  for(i in 1:n){
    #Longitudinal observations
    Y1[i]~dnorm(mu1[i],tau1)
    mu1[i]<-inprod(betaL[],X[i,])+inprod(b[id[i],],Z[i,])
  }
  for(k in 1:n2){
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
    }}
  #Survival and censoring times
  #Hazard function
  for(k in 1:n2){
    linearpred[k]<-inprod(betaL[nindtime],Xv[k,])+b[k,1]+(betaL[indtime[1]]+b[k,2])*Time[k]+
      (betaL[indtime[2]]+b[k,3])*pow(Time[k],2)


    Alpha0[k]<-  inprod(betaS[1:NbetasS],XS[k,1:NbetasS])+
                 gamma1*(inprod(betaL[nindtime],Xv[k,])+b[k,1])
    Alpha1[k]<- gamma1*(betaL[indtime[1]]+b[k,2])
    Alpha2[k]<- gamma1*(betaL[indtime[2]]+b[k,3])

    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
                  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
                  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
                  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
                  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k]+Alpha2[k]*pow(Time[k],2))


    for(j in 1:K){
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
                        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
                        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
                        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
                        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j]+Alpha2[k]*pow(xk11[k,j],2))

        }

    logSurv[k]<- -inprod(wk11[k,],chaz[k,])
   phi[k]<- KK-Death[k]*log(haz[k])-logSurv[k]


    #Definition of the survival log-likelihood using zeros trick
    zeros[k]~dpois(phi[k])
    #Random effects
    b[k,1:Nb]~dmnorm(mub[],Omega[,])

  }
    Omega[1:Nb,1:Nb]<-inverse(Sigma[,])
  tau1<-1/sigma1
}"

  ########
  tmp <- dataSurv[all.vars(formSurv)]
  Time <- tmp[all.vars(formSurv)][, 1] # matrix of observed time such as Time=min(Tevent,Tcens)
  Death <- tmp[all.vars(formSurv)][, 2] # vector of event indicator (delta)
  nTime <- length(Time) # number of subject having Time
  # design matrice
  suppressWarnings({
    mfZ <- stats::model.frame(formSurv, data = tmp, na.action = NULL)
  })
  XS <- stats::model.matrix(formSurv, mfZ, na.action = NULL)[, -1]

  glq <- statmod::gauss.quad(15, kind = "legendre")
  xk <- glq$nodes # Nodes
  wk <- glq$weights # Weights
  K <- length(xk) # K-points
  ################




  data_Long_s <- dataLong[dataLong$obstime <= s, ]



  DP_tot <- NULL
  for (ttt in 1:mi) {
    X <- Z <- Xv <- Zv <- Nb <- list()
    indB <- indtime <- list()

    bhat_mean <- bhat_chain <- list()
    for (j in c(1:nmark)) {
      if (model[[j]] == "intercept") {
        data_long <- data_Long_s[unique(c(all.vars(formGroup[[j]]), all.vars(formFixed[[j]]), all.vars(formRandom[[j]])))]
        y <- data_long[all.vars(formFixed[[j]])][, 1]
        # data_long <- data_long[is.na(y) == FALSE, ]
        # y <- data_long[all.vars(formFixed[[j]])][, 1]

        mfX <- stats::model.frame(formFixed[[j]], data = data_long, na.action = NULL)
        X[[j]] <- stats::model.matrix(formFixed[[j]], mfX, na.action = NULL)
        mfU <- stats::model.frame(formRandom[[j]], data = data_long, na.action = NULL)
        id <- as.integer(data_long[all.vars(formGroup[[j]])][, 1])
        n2 <- length(unique(id))
        n <- length(id)

        M <- table(id)
        id2 <- rep(1:length(M), M)

        Obstime <- Obstime
        Xvtime <- cbind(id, X[[j]][, colnames(X[[j]]) %in% setdiff(colnames(X[[j]]), Obstime)])
        Xv[[j]] <- Xvtime[!duplicated(Xvtime), -1] ### X of data without time and id replications

        indB[[j]] <- 1:dim(X[[j]])[2]
        indtime[[j]] <- indB[[j]][colnames(X[[j]]) %in% Obstime] # index of time
      }
      ####
      if (model[[j]] == "linear") {
        data_long <- data_Long_s[unique(c(all.vars(formGroup[[j]]), all.vars(formFixed[[j]]), all.vars(formRandom[[j]])))]
        y <- data_long[all.vars(formFixed[[j]])][, 1]
        # data_long <- data_long[is.na(y) == FALSE, ]
        # y <- data_long[all.vars(formFixed[[j]])][, 1]

        mfX <- stats::model.frame(formFixed[[j]], data = data_long, na.action = NULL)
        X[[j]] <- stats::model.matrix(formFixed[[j]], mfX, na.action = NULL)
        mfU <- stats::model.frame(formRandom[[j]], data = data_long, na.action = NULL)
        Z[[j]] <- stats::model.matrix(formRandom[[j]], mfU, na.action = NULL)
        id <- as.integer(data_long[all.vars(formGroup[[j]])][, 1])
        n <- length(id)
        M <- table(id)
        id2 <- rep(1:length(M), M)
        n2 <- length(unique(id))
        Nb[[j]] <- dim(Z[[j]])[2]
        Obstime <- Obstime
        Xvtime <- cbind(id, X[[j]][, colnames(X[[j]]) %in% setdiff(colnames(X[[j]]), Obstime)])
        Xv[[j]] <- Xvtime[!duplicated(Xvtime), -1]

        indB[[j]] <- 1:dim(X[[j]])[2]
        indtime[[j]] <- indB[[j]][colnames(X[[j]]) %in% Obstime] # index of time
      }
      #############
      if (model[[j]] == "quadratic") {
        data_long <- data_Long_s[unique(c(all.vars(formGroup[[j]]), all.vars(formFixed[[j]]), all.vars(formRandom[[j]])))]
        y <- data_long[all.vars(formFixed[[j]])][, 1]
        # data_long <- data_long[is.na(y) == FALSE, ]
        # y <- data_long[all.vars(formFixed[[j]])][, 1]


        mfX <- stats::model.frame(formFixed[[j]], data = data_long, na.action = NULL)
        X[[j]] <- stats::model.matrix(formFixed[[j]], mfX, na.action = NULL)
        mfU <- stats::model.frame(formRandom[[j]], data = data_long, na.action = NULL)
        Z[[j]] <- stats::model.matrix(formRandom[[j]], mfU, na.action = NULL)


        colnamesmfX <- colnames(X[[j]])
        colnamesmfU <- colnames(Z[[j]])
        Obstime <- Obstime

        Xtime2 <- mfX[, Obstime]^2

        X[[j]] <- cbind(X[[j]], Xtime2)
        colnames(X[[j]]) <- c(colnamesmfX, "obstime2")
        Z[[j]] <- cbind(Z[[j]], Xtime2)
        colnames(Z[[j]]) <- c(colnamesmfU, "obstime2")
        Obstime2 <- "obstime2"

        id <- as.integer(data_long[all.vars(formGroup[[j]])][, 1])
        M <- table(id)
        id2 <- rep(1:length(M), M)

        n <- length(id)
        n2 <- length(unique(id))
        Nb[[j]] <- dim(Z[[j]])[2]

        Obstime2n <- c(Obstime, Obstime2)
        Xvtime <- cbind(id, X[[j]][, colnames(X[[j]]) %in% setdiff(colnames(X[[j]]), Obstime2n)])
        Xv[[j]] <- Xvtime[!duplicated(Xvtime), -1] ### X of data without time and id replications


        indB[[j]] <- 1:dim(X[[j]])[2]
        indtime[[j]] <- indB[[j]][colnames(X[[j]]) %in% Obstime2n] # index of time
      }
      # }






      ########  BUGS code  ########
      M1 <- table(id)
      NbetasS <- dim(XS)[2]


      if (model[[j]] == "intercept") {
        i.jags <- function() {
          list(
            b = rep(0, n2)
          )
        }
      } else {
        i.jags <- function() {
          list(
            b = matrix(0, n2, Nb[[j]])
          )
        }
      }
      parameters <- c("b")
      ############################### ????
      betaLmc <- object$sim_step1[[j]]$sim$beta
      betaSmc <- object$sim_step1[[j]]$sim$gamma
      gamma1mc <- object$sim_step1[[j]]$sim$alpha
      sigma1mc <- object$sim_step1[[j]]$sim$sigma
      Sigmamc <- object$sim_step1[[j]]$sim$Sigma
      hmc <- object$sim_step1[[j]]$sim$h
      #### Data
      Sample <- sample(1:length(sigma1mc), mi)

      ######

      betaLmc <- object$sim_step1[[j]]$sim$beta
      betaSmc <- object$sim_step1[[j]]$sim$gamma
      gamma1mc <- object$sim_step1[[j]]$sim$alpha
      sigma1mc <- object$sim_step1[[j]]$sim$sigma
      Sigmamc <- object$sim_step1[[j]]$sim$Sigma
      hmc <- object$sim_step1[[j]]$sim$h
      #### Data
      set.seed(3)
      Sample <- sample(1:length(sigma1mc), mi)


      if (is.matrix(betaLmc) == TRUE) {
        betaLmc <- betaLmc[Sample, ]
      } else {
        betaLmc <- betaLmc[Sample]
      }
      gamma1mc <- gamma1mc[Sample]
      sigma1mc <- sigma1mc[Sample]

      betaSmc <- betaSmc[Sample, ]


      if (is.array(Sigmamc) == TRUE) {
        Sigmamc <- Sigmamc[Sample, , ]
      } else {
        Sigmamc <- Sigmamc[Sample]
      }

      hmc <- hmc[Sample, ]



      if (is.matrix(betaLmc) == TRUE) {
        betaL <- betaLmc[ttt, ]
      } else {
        betaL <- betaLmc[ttt]
      }
      gamma1 <- gamma1mc[ttt]
      sigma1 <- sigma1mc[ttt]

      betaS <- betaSmc[ttt, ]


      if (is.array(Sigmamc) == TRUE) {
        Sigma <- Sigmamc[ttt, , ]
      } else {
        Sigma <- Sigmamc[ttt]
      }

      h <- hmc[ttt, ]












      NbetasS <- dim(XS)[2]

      if (is.matrix(Xv[[j]]) == FALSE) {
        if (model[[j]] == "intercept") {
          model_fileLb_last <- textConnection(model_fileI1b)
          d.jags <- list(
            betaL = betaL, betaS = betaS,
            gamma1 = gamma1, sigma1 = sigma1,
            Sigma = Sigma, h = h,
            n = n, Time = rep(s, n2), Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2],
            X = X[[j]], id = id2, indtime = indtime[[j]],
            Death = Death, zeros = rep(0, n2),
            s = peice, xk = xk, wk = wk, K = K, KK = KK
          )
        }

        if (model[[j]] == "linear") {
          model_fileLb_last <- textConnection(model_fileL1b)
          d.jags <- list(
            betaL = betaL, betaS = betaS,
            gamma1 = gamma1, sigma1 = sigma1,
            Sigma = Sigma, h = h,
            n = n, Time = rep(s, n2), Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2],
            X = X[[j]], Z = Z[[j]], id = id2, indtime = indtime[[j]],
            Death = Death, mub = rep(0, Nb[[j]]), Nb = Nb[[j]], zeros = rep(0, n2),
            s = peice, xk = xk, wk = wk, K = K, KK = KK
          )
        }
        if (model[[j]] == "quadratic") {
          model_fileLb_last <- textConnection(model_fileQ1b)
          d.jags <- list(
            betaL = betaL, betaS = betaS,
            gamma1 = gamma1, sigma1 = sigma1,
            Sigma = Sigma, h = h,
            n = n, Time = rep(s, n2), Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2],
            X = X[[j]], Z = Z[[j]], id = id2, indtime = indtime[[j]],
            Death = Death, mub = rep(0, Nb[[j]]), Nb = Nb[[j]], zeros = rep(0, n2),
            s = peice, xk = xk, wk = wk, K = K, KK = KK
          )
        }
      } else {
        if (model[[j]] == "intercept") {
          model_fileLb_last <- textConnection(model_fileIb)
          d.jags <- list(
            betaL = betaL, betaS = betaS,
            gamma1 = gamma1, sigma1 = sigma1,
            Sigma = Sigma, h = h,
            n = n, Time = rep(s, n2), Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2],
            X = X[[j]], id = id2, Xv = Xv[[j]], indtime = indtime[[j]], nindtime = c(1:dim(X[[j]])[2])[-indtime[[j]]],
            Death = Death, zeros = rep(0, n2),
            s = peice, xk = xk, wk = wk, K = K, KK = KK
          )
        }
        if (model[[j]] == "linear") {
          model_fileLb_last <- textConnection(model_fileLb)
          d.jags <- list(
            betaL = betaL, betaS = betaS,
            gamma1 = gamma1, sigma1 = sigma1,
            Sigma = Sigma, h = h,
            n = n, Time = rep(s, n2), Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2],
            X = X[[j]], Z = Z[[j]], id = id2, Xv = Xv[[j]], indtime = indtime[[j]], nindtime = c(1:dim(X[[j]])[2])[-indtime[[j]]],
            Death = Death, mub = rep(0, Nb[[j]]), Nb = Nb[[j]], zeros = rep(0, n2),
            s = peice, xk = xk, wk = wk, K = K, KK = KK
          )
        }
        if (model[[j]] == "quadratic") {
          model_fileLb_last <- textConnection(model_fileQb)
          d.jags <- list(
            betaL = betaL, betaS = betaS,
            gamma1 = gamma1, sigma1 = sigma1,
            Sigma = Sigma, h = h,
            n = n, Time = rep(s, n2), Y1 = y, n2 = n2, XS = XS, NbetasS = dim(XS)[2],
            X = X[[j]], Z = Z[[j]], id = id2, Xv = Xv[[j]], indtime = indtime[[j]],
            nindtime = c(1:dim(X[[j]])[2])[-indtime[[j]]],
            Death = Death, mub = rep(0, Nb[[j]]), Nb = Nb[[j]], zeros = rep(0, n2),
            s = peice, xk = xk, wk = wk, K = K, KK = KK
          )
        }
      }
      sim1 <- jagsUI::jags(
        data = d.jags,
        inits = i.jags,
        parameters.to.save = parameters,
        model.file = model_fileLb_last,
        n.chains = n.chains,
        parallel = FALSE,
        n.adapt = FALSE,
        n.iter = n.iter,
        n.burnin = n.burnin,
        n.thin = n.thin,
        DIC = TRUE
      )

      bhat_mean[[j]] <- sim1$mean$b
      bhat_chain[[j]] <- sim1$sims.list$b
    }





    ################################### ???
   n2 <- length(unique(id))
    sigma <- c()
    # mu1 <- matrix(0, n2, nmark)
    betaL <- b <- list()

    for (j in c(1:nmark)) {
      betaL[[j]] <- object$sim_step1[[j]]$PMean$beta
      sigma <- append(sigma, object$sim_step1[[j]]$PMean$sigma)
      # mu1[,j]<- object$sim_step1[[j]]$PMean$linearpred
      b[[j]] <- bhat_mean[[j]]
    }
    indtime <- nindtime <- list()
    for (j in c(1:nmark)) {
      indB <- 1:dim(X[[j]])[2]
      if (model[[j]] == "quadratic") {
        indtime[[j]] <- indB[colnames(X[[j]]) %in% c(Obstime, Obstime2)]
      } else {
        indtime[[j]] <- indB[colnames(X[[j]]) %in% Obstime]
      }
      nindtime[[j]] <- c(1:dim(X[[j]])[2])[-indtime[[j]]]
    }
    LP1 <- LP2 <- LP3 <- matrix(0, n2, nmark)
    for (i in 1:n2) {
      for (j in c(1:nmark)) {
        if (is.matrix(Xv[[j]]) == TRUE) {
          if (model[[j]] == "intercept") {
            LP1[i, j] <- betaL[[j]][nindtime[[j]]] %*% Xv[[j]][i, ] + b[[j]][i]
          } else {
            LP1[i, j] <- betaL[[j]][nindtime[[j]]] %*% Xv[[j]][i, ] + b[[j]][i, 1]
          }
        } else {
          if (model[[j]] == "intercept") {
            LP1[i, j] <- betaL[[j]][nindtime[[j]]] + b[[j]][i]
          } else {
            LP1[i, j] <- betaL[[j]][nindtime[[j]]] + b[[j]][i, 1]
          }
        }
        if (model[[j]] == "intercept") {
          LP2[i, j] <- betaL[[j]][indtime[[j]][1]]
        } else {
          LP2[i, j] <- betaL[[j]][indtime[[j]][1]] + b[[j]][i, 2]
        }
        LP3[i, j] <- 0


        if (model[[j]] == "quadratic") (LP3[i, j] <- betaL[[j]][indtime[[j]][2]] + b[[j]][i, 3])
      }
    }

    ########################### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # betaS <- object$Estimation$Survival_model$gamma$Est
    # alpha <- object2$Estimation$Survival_model$alpha$Est
    # h <- object2$Estimation$Survival_model$lambda$Est

    betaS <- object$sim_step2$mean$betaS
    alpha <- object$sim_step2$mean$alpha
    h <- object$sim_step2$mean$h
    ########################### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    step <- function(x) {
      z <- 0
      if (x >= 0) (z <- 1)
      z
    }
    inprod <- function(a, b) {
      z <- a %*% b
      z
    }
    Surv_d <- c()
    chaz <- matrix(0, n2, K)
    Alpha0 <- Alpha1 <- Alpha2 <- c()
    DP <- c()

    for (k in 1:n2) {
      Alpha0[k] <- betaS %*% XS[k, ] + alpha %*% LP1[k, ]
      Alpha1[k] <- alpha %*% LP2[k, ]
      Alpha2[k] <- alpha %*% LP3[k, ]




      xk11 <- wk11 <- c()

      for (j in 1:K) {
        # Scaling Gauss-Kronrod/Legendre quadrature
        xk11[j] <- (xk[j] + 1) / 2 * s
        wk11[j] <- wk[j] * s / 2
        #  Hazard function at Gauss-Kronrod/Legendre nodes
        chaz[k, j] <- ((h[1] * step(peice[1] - xk11[j])) +
          (h[2] * step(xk11[j] - peice[1]) * step(peice[2] - xk11[j])) +
          (h[3] * step(xk11[j] - peice[2]) * step(peice[3] - xk11[j])) +
          (h[4] * step(xk11[j] - peice[3]) * step(peice[4] - xk11[j])) +
          (h[5] * step(xk11[j] - peice[4]))) * exp(Alpha0[k] + Alpha1[k] * xk11[j] + Alpha2[k] * xk11[j]^2)
      }


      Surv_d[k] <- exp(-inprod(wk11, chaz[k, ]))
    }

    Surv_n <- c()
    chaz <- matrix(0, n2, K)
    for (k in 1:n2) {
      xk11 <- wk11 <- c()

      for (j in 1:K) {
        # Scaling Gauss-Kronrod/Legendre quadrature
        xk11[j] <- (xk[j] + 1) / 2 * (s + Dt)
        wk11[j] <- wk[j] * (s + Dt) / 2
        #  Hazard function at Gauss-Kronrod/Legendre nodes
        chaz[k, j] <- ((h[1] * step(peice[1] - xk11[j])) +
          (h[2] * step(xk11[j] - peice[1]) * step(peice[2] - xk11[j])) +
          (h[3] * step(xk11[j] - peice[2]) * step(peice[3] - xk11[j])) +
          (h[4] * step(xk11[j] - peice[3]) * step(peice[4] - xk11[j])) +
          (h[5] * step(xk11[j] - peice[4]))) * exp(Alpha0[k] + Alpha1[k] * xk11[j] + +Alpha2[k] * xk11[j]^2)
      }

      Surv_n[k] <- exp(-inprod(wk11, chaz[k, ]))
    }
    Surv_d[Surv_d == 0] <- 0.000001
    DP <- 1 - Surv_n / Surv_d


    DP_tot <- rbind(DP_tot, DP)
  }
  #####################


  DP <- apply(DP_tot, 2, mean)
  DPQ <- t(apply(DP_tot, 2, quantile, c(0.025, 0.975)))

  DP_last <- cbind(unique(id), DP, DPQ)
  colnames(DP_last) <- c("id", "est", "lower", "upper")
  DP_last <- data.frame(DP_last)
  list(DP = DP_last, s = s, t = Dt)
}
