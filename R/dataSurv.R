#' Simulated competing risks data
#'
#' Simulated survival data were generated in which "survtime" follows a survival time with a event indicator called "CR" and some covariates.
#'
#' @docType data
#'
#'
#' @name dataSurv
#' @format A data frame which contains id, survtime, CR, w1, x1, w2, x2.
#' \describe{
#'   \item{id}{patients identifier}
#'   \item{survtime}{survival time (the response variable)}
#'   \item{death}{event indicator, 1=death, 0=alive (censored)}
#'   \item{x1}{a binary covariate}
#'   \item{x2}{a continuous covariate}
#' }
#' @seealso \code{\link{UJM}}
"dataSurv"

