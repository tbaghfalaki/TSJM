#' Two stage Joint Modeling
#'
#' @description
#' Employ a two-stage approach to efficiently select variables for joint modeling of multiple longitudinal markers and time-to-event outcomes within a Bayesian framework.
#'
#' @details
#' A two-stage approach within a Bayesian framework is designed to handle multiple longitudinal measurements and survival outcomes. In the first stage, we address the estimation of a one-marker joint model for the event along with each longitudinal marker. Leveraging these estimates, we derive predictions for the expected values or slopes of individual marker trajectories. In the second stage, we utilize a proportional hazard model that incorporates the expected current values and/or slopes of all markers as time-dependent covariates.
#'
#' @param formFixed a list of formulas for fixed part of longitudinal model
#' @param formRandom a list of formulas for random part of longitudinal model
#' @param formGroup a list of formulas specifying the cluster variable for Y (e.g. = list (~ subject, ~ subject,...))
#' @param formSurv formula for survival model
#' @param dataLong data set of observed longitudinal variables.
#' @param dataSurv data set of observed survival variables.
#' @param nmark the number of longitudinal markers
#' @param K1 Number of nodes and weights for calculating Gaussian quadrature in the first stage.
#' @param K2 Number of nodes and weights for calculating Gaussian quadrature in the second stage.
#' @param model a list of the models for the longitudinal part which includes "linear" or "quadratic".
#' @param Obstime the observed time in longitudinal data
#' @param ncl the number of nodes to be forked for parallel computing
#' @param Limp the number of multiple imputation; default is 10.
#' @param n.chains1 the number of parallel chains for the model in the first stage; default is 1.
#' @param n.iter1 integer specifying the total number of iterations in the first stage; default is 1000.
#' @param n.burnin1 integer specifying how many of n.iter to discard as burn-in in the first stage; default is 5000.
#' @param n.thin1 integer specifying the thinning of the chains in the first stage; default is 1.
#' @param DIC Logical; if TRUE (default), compute deviance, pD, and DIC. The rule pD=var(deviance) / 2 is used.
#' @param quiet Logical, whether to suppress stdout in jags.model().
#'
#' @importFrom stats quantile rnorm model.frame model.matrix
#'
#'
#' @return
#' - MCMC chains for the unknown parameters
#' - mu.vect list of posterior mean for each parameter
#' - sd.vect list of standard error for each parameter
#' - 2.5% list of posterior mode for each parameter
#' - 97.5% list of posterior median for each parameter
#' - Rhat Gelman and Rubin diagnostic for all parameter
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}
#'
#'
#' @example inst/exampleTSC.R
#'
#' @md
#' @export


TSC <- function(formFixed, formRandom, formGroup, formSurv, nmark, K1 = K1, K2 = K2,
                model = model, n.chains1 = n.chains1, n.iter1 = n.iter1, n.burnin1 = floor(n.iter1 / 2),
                n.thin1 = max(1, floor((n.iter1 - n.burnin1) / 1000)),
                Obstime = "obstime", ncl = ncl, Limp = 10,
                DIC = TRUE, quiet = FALSE, dataLong, dataSurv) {
  j <- 1:nmark
  boot_fx <- function(j) {
    A1 <- UJM(
      formFixed = formFixed[[j]], formRandom = formRandom[[j]],
      formGroup = formGroup[[j]], formSurv = formSurv, dataLong = dataLong,
      dataSurv = dataSurv, K = K1, model = model[[j]], Obstime = Obstime,
      n.chains = n.chains1, n.iter = n.iter1, n.burnin = n.burnin1,
      n.thin = n.thin1,
      DIC = DIC, quiet = quiet
    )
    list(sim = A1$MCMC, PMean = A1$PMean, Long = A1$Estimation)
  }

  results <- parallelsugar::mclapply(j, boot_fx, mc.cores = ncl)


  #############
  gamma <- sigma <- c()
  K <- K2
  X <- Z <- Xv <- Zv <- Nb <- list()
  indB <- indtime <- list()
  for (j in 1:nmark) {
    if (model[[j]] == "intercept") {
      data_long <- dataLong[unique(c(all.vars(formGroup[[j]]), all.vars(formFixed[[j]]), all.vars(formRandom[[j]])))]
      y <- data_long[all.vars(formFixed[[j]])][, 1]
      data_long <- data_long[is.na(y) == FALSE, ]
      y <- data_long[all.vars(formFixed[[j]])][, 1]
      mfX <- stats::model.frame(formFixed[[j]], data = data_long) # , na.action = NULL)
      X[[j]] <- stats::model.matrix(formFixed[[j]], mfX) # , na.action = NULL)
      mfU <- stats::model.frame(formRandom[[j]], data = data_long) # , na.action = NULL)
      id0 <- as.integer(data_long[all.vars(formGroup[[j]])][, 1])
      n2 <- length(unique(id0))

      M <- table(id0)
      id02 <- rep(1:length(M), M)

      Obstime <- Obstime
      Xvtime <- cbind(id0, X[[j]][, colnames(X[[j]]) %in% setdiff(colnames(X[[j]]), Obstime)])
      Xv[[j]] <- Xvtime[!duplicated(Xvtime), -1] ### X of data without time and id0 replications

      indB[[j]] <- 1:dim(X[[j]])[2]
      indtime[[j]] <- indB[[j]][colnames(X[[j]]) %in% Obstime] # index of time
    }
    ####
    if (model[[j]] == "linear") {
      data_long <- dataLong[unique(c(all.vars(formGroup[[j]]), all.vars(formFixed[[j]]), all.vars(formRandom[[j]])))]
      y <- data_long[all.vars(formFixed[[j]])][, 1]
      data_long <- data_long[is.na(y) == FALSE, ]
      y <- data_long[all.vars(formFixed[[j]])][, 1]

      mfX <- stats::model.frame(formFixed[[j]], data = data_long, na.action = NULL)
      X[[j]] <- stats::model.matrix(formFixed[[j]], mfX, na.action = NULL)
      mfU <- stats::model.frame(formRandom[[j]], data = data_long, na.action = NULL)
      Z[[j]] <- stats::model.matrix(formRandom[[j]], mfU, na.action = NULL)
      id0 <- as.integer(data_long[all.vars(formGroup[[j]])][, 1])

      n2 <- length(unique(id0))
      Nb[[j]] <- dim(Z[[j]])[2]
      Obstime <- Obstime
      Xvtime <- cbind(id0, X[[j]][, colnames(X[[j]]) %in% setdiff(colnames(X[[j]]), Obstime)])
      Xv[[j]] <- Xvtime[!duplicated(Xvtime), -1]

      indB[[j]] <- 1:dim(X[[j]])[2]
      indtime[[j]] <- indB[[j]][colnames(X[[j]]) %in% Obstime] # index of time
    }
    #############
    if (model[[j]] == "quadratic") {
      data_long <- dataLong[unique(c(all.vars(formGroup[[j]]), all.vars(formFixed[[j]]), all.vars(formRandom[[j]])))]
      y <- data_long[all.vars(formFixed[[j]])][, 1]
      data_long <- data_long[is.na(y) == FALSE, ]
      y <- data_long[all.vars(formFixed[[j]])][, 1]


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

      id0 <- as.integer(data_long[all.vars(formGroup[[j]])][, 1])


      n2 <- length(unique(id0))
      Nb[[j]] <- dim(Z[[j]])[2]

      Obstime2n <- c(Obstime, Obstime2)
      Xvtime <- cbind(id0, X[[j]][, colnames(X[[j]]) %in% setdiff(colnames(X[[j]]), Obstime2n)])
      Xv[[j]] <- Xvtime[!duplicated(Xvtime), -1] ### X of data without time and id0 replications


      indB[[j]] <- 1:dim(X[[j]])[2]
      indtime[[j]] <- indB[[j]][colnames(X[[j]]) %in% Obstime2n] # index of time
    }
  }
  ##########
  n <- length(id0)

  tmp <- dataSurv[all.vars(formSurv)]
  Time <- tmp[all.vars(formSurv)][, 1] # matrix of observed time such as Time=min(Tevent,Tcens)
  Death <- tmp[all.vars(formSurv)][, 2] # vector of event indicator (delta)
  nTime <- length(Time) # number of subject having Time
  peice <- quantile(Time, seq(.2, 0.8, length = 4))

   # design matrice
  suppressWarnings({
    mfZ <- stats::model.frame(formSurv, data = tmp, na.action = NULL)
  })
  XS <- stats::model.matrix(formSurv, mfZ, na.action = NULL)[, -1]
  #########
  n2 <- dim(dataSurv)[1]
  gamma <- sigma <- c()
  mu1 <- matrix(0, n2, nmark)
  betaL <- b <- list()
  for (j in 1:nmark) {
    betaL[[j]] <- results[[j]]$PMean$beta
    gamma <- append(gamma, results[[j]]$PMean$alpha)
    sigma <- append(sigma, results[[j]]$PMean$sigma)
    mu1[, j] <- results[[j]]$PMean$linearpred
    b[[j]] <- results[[j]]$PMean$b
  }

  indtime <- nindtime <- list()
  for (j in 1:nmark) {
    indB <- 1:dim(X[[j]])[2]
    if (model[[j]] == "intercept") {
      indtime[[j]] <- indB[colnames(X[[j]]) %in% Obstime]
    }
    if (model[[j]] == "linear") {
      indtime[[j]] <- indB[colnames(X[[j]]) %in% Obstime]
    }
    if (model[[j]] == "quadratic") {
      indtime[[j]] <- indB[colnames(X[[j]]) %in% c(Obstime, Obstime2)]
    }
    nindtime[[j]] <- c(1:dim(X[[j]])[2])[-indtime[[j]]]
  }


  Lp1 <- Lp2 <- Lp3 <- matrix(0, n2, nmark)
  for (i in 1:n2) {
    for (j in 1:nmark) {
      if (is.matrix(Xv[[j]]) == TRUE) {
        if (model[[j]] != "intercept") {
          Lp1[i, j] <- betaL[[j]][nindtime[[j]]] %*% Xv[[j]][i, ] + b[[j]][i, 1]
        } else {
          Lp1[i, j] <- betaL[[j]][nindtime[[j]]] %*% Xv[[j]][i, ] + b[[j]][i]
        }
      } else {
        if (model[[j]] != "intercept") {
          Lp1[i, j] <- betaL[[j]][nindtime[[j]]] + b[[j]][i, 1]
        } else {
          Lp1[i, j] <- betaL[[j]][nindtime[[j]]] + b[[j]][i]
        }
      }
      if (model[[j]] != "intercept") {
        Lp2[i, j] <- betaL[[j]][indtime[[j]][1]] + b[[j]][i, 2]
      } else {
        Lp2[i, j] <- betaL[[j]][indtime[[j]][1]]
      }
      Lp3[i, j] <- 0


      if (model[[j]] == "quadratic") (Lp3[i, j] <- betaL[[j]][indtime[[j]][2]] + b[[j]][i, 3])
    }
  }

  LP1 <- matrix(rep(Lp1[1, ], M[1]), ncol = ncol(Lp1), byrow = TRUE)
  LP2 <- matrix(rep(Lp2[1, ], M[1]), ncol = ncol(Lp2), byrow = TRUE)
  LP3 <- matrix(rep(Lp3[1, ], M[1]), ncol = ncol(Lp3), byrow = TRUE)
  for (i in 2:n2) {
    LP1 <- rbind(LP1, matrix(rep(Lp1[i, ], M[i]), ncol = ncol(Lp1), byrow = TRUE))
    LP2 <- rbind(LP2, matrix(rep(Lp2[i, ], M[i]), ncol = ncol(Lp2), byrow = TRUE))
    LP3 <- rbind(LP3, matrix(rep(Lp3[i, ], M[i]), ncol = ncol(Lp3), byrow = TRUE))
  }

  lPredY <- LP1 + LP2 * data_long[, Obstime] + LP3 * data_long[, Obstime]^2

  ###################################
  indexchain <- 1:length(results[[1]]$sim$sigma)

  samplec <- sample(indexchain, Limp)

  n2 <- dim(dataSurv)[1]
  LPredY <- list()
  rrr <- 0
  for (samplec1 in samplec) {
    rrr <- rrr + 1
    gamma <- sigma <- c()
    mu1 <- matrix(0, n2, nmark)
    betaL <- b <- list()

    for (j in 1:nmark) {
      betaL[[j]] <- results[[j]]$sim$beta[samplec1, ]
      gamma <- append(gamma, results[[j]]$sim$alpha[samplec1])
      sigma <- append(sigma, results[[j]]$sim$sigma[samplec1])

      if (is.matrix(results[[j]]$sim$b) == TRUE) {
        b[[j]] <- results[[j]]$sim$b[samplec1, ]
      } else {
        b[[j]] <- results[[j]]$sim$b[samplec1, , ]
      }
    }

    indtime <- nindtime <- list()
    for (j in 1:nmark) {
      indB <- 1:dim(X[[j]])[2]
      if (model[[j]] == "intercept") {
        indtime[[j]] <- indB[colnames(X[[j]]) %in% Obstime]
      }
      if (model[[j]] == "linear") {
        indtime[[j]] <- indB[colnames(X[[j]]) %in% Obstime]
      }
      if (model[[j]] == "quadratic") {
        indtime[[j]] <- indB[colnames(X[[j]]) %in% c(Obstime, Obstime2)]
      }
      nindtime[[j]] <- c(1:dim(X[[j]])[2])[-indtime[[j]]]
    }


    Lp1 <- Lp2 <- Lp3 <- matrix(0, n2, nmark)
    for (i in 1:n2) {
      for (j in 1:nmark) {
        if (is.matrix(Xv[[j]]) == TRUE) {
          if (model[[j]] != "intercept") {
            Lp1[i, j] <- betaL[[j]][nindtime[[j]]] %*% Xv[[j]][i, ] + b[[j]][i, 1]
          } else {
            Lp1[i, j] <- betaL[[j]][nindtime[[j]]] %*% Xv[[j]][i, ] + b[[j]][i]
          }
        } else {
          if (model[[j]] != "intercept") {
            Lp1[i, j] <- betaL[[j]][nindtime[[j]]] + b[[j]][i, 1]
          } else {
            Lp1[i, j] <- betaL[[j]][nindtime[[j]]] + b[[j]][i]
          }
        }
        if (model[[j]] != "intercept") {
          Lp2[i, j] <- betaL[[j]][indtime[[j]][1]] + b[[j]][i, 2]
        } else {
          Lp2[i, j] <- betaL[[j]][indtime[[j]][1]]
        }
        Lp3[i, j] <- 0


        if (model[[j]] == "quadratic") (Lp3[i, j] <- betaL[[j]][indtime[[j]][2]] + b[[j]][i, 3])
      }
    }

    LP1 <- matrix(rep(Lp1[1, ], M[1]), ncol = ncol(Lp1), byrow = TRUE)
    LP2 <- matrix(rep(Lp2[1, ], M[1]), ncol = ncol(Lp2), byrow = TRUE)
    LP3 <- matrix(rep(Lp3[1, ], M[1]), ncol = ncol(Lp3), byrow = TRUE)
    for (i in 2:n2) {
      LP1 <- rbind(LP1, matrix(rep(Lp1[i, ], M[i]), ncol = ncol(Lp1), byrow = TRUE))
      LP2 <- rbind(LP2, matrix(rep(Lp2[i, ], M[i]), ncol = ncol(Lp2), byrow = TRUE))
      LP3 <- rbind(LP3, matrix(rep(Lp3[i, ], M[i]), ncol = ncol(Lp3), byrow = TRUE))
    }

    lPredY1 <- LP1 + LP2 * data_long[, Obstime] + LP3 * data_long[, Obstime]^2

    LPredY[[rrr]] <- lPredY1
  }














  Longitudinal <- list()
  for (j in 1:nmark) {
    Longitudinal[[j]] <- results[[j]]$Long
  }


  list(
    formFixed = formFixed, formRandom = formRandom, formGroup = formGroup, formSurv = formSurv,
    model = model, Obstime = Obstime,
    mu1 = mu1,peice=peice,
    nmark = nmark, XS = XS, lPredY = lPredY,
    sim_step1 = results, Longitudinal = Longitudinal, LPredY = LPredY
  )
}
