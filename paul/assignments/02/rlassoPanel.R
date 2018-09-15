#' Rigorous Lasso for panel data
#'
#' Estimate the Lasso for panel data with fixed effects using the
#' approach of Belloni, Chernozhukov, Hansen, & Kozbur (2016). Based
#' on rlasso from the hdm package.
#'
#'
#' @param x matrix of regressors
#' @param y dependent variable
#' @param f group identifier for fixed effects
#' @param post Whether to use post lasso
#' @param penalty list with options for the calculation of the penalty.
#' \itemize{
#' \item{\code{c} and \code{gamma}}{ constants for the penalty with default \code{c=1.1} and \code{gamma=0.1}}
#' \item{\code{homoscedastic}}{ logical, if homoscedastic errors are considered (default \code{FALSE}). Option \code{none} is described below.}
#' \item{\code{X.dependent.lambda}}{ logical,  \code{TRUE}, if the penalization parameter depends on the the design of the matrix \code{x}. \code{FALSE}, if independent of the design matrix  (default).}
#' \item{\code{numSim}}{ number of simulations for the dependent methods, default=5000}
#' \item{\code{lambda.start}}{ initial penalization value, compulsory for method "none"}
#' }
#' @param control list with control values.
#' \code{numIter} number of iterations for the algorithm for
#' the estimation of the variance and data-driven penalty, ie. loadings,
#' \code{tol} tolerance for improvement of the estimated variances.
#'\code{threshold} is applied to the final estimated lasso
#' coefficients. Absolute values below the threshold are set to zero.
#'
#' @return
rlassoPanel <-
  function(x, y, f, post=TRUE,
           penalty = list(homoscedastic=FALSE,
                          X.dependent.lambda=FALSE,
                          lambda.start = NULL, c = 1.1,
                          gamma=.1/log(n)),
           control = list(numIter = 15, tol = 10^-5, threshold =
                                                       NULL))
{
  x <- as.matrix(x)
  y <- as.matrix(y)

    n <- dim(x)[1]
  p <- dim(x)[2]

  if (is.null(colnames(x))) colnames(x) <- paste("V", 1:p, sep = "")

  ind.names <- 1:p
  # set options to default values if missing
  if (!exists("homoscedastic", where = penalty))  penalty$homoscedastic = "FALSE"
  if (!exists("X.dependent.lambda", where = penalty))  penalty$X.dependent.lambda = "FALSE"
  if (!exists("gamma", where = penalty))  penalty$gamma = 0.1/log(n)
  if (penalty$homoscedastic=="none" & !exists("lambda.start", where=penalty)) stop("lambda.start must be provided!")
  # checking input numIter, tol
  if (!exists("numIter", where = control)) {
    control$numIter = 15
  }
  if (!exists("tol", where = control)) {
    control$tol = 10^-5
  }
  if (post==FALSE & (!exists("c", where = penalty))) {
    penalty$c = 0.5
  }


  # remove na's
  drop <- is.na(rowSums(cbind(x,y)))
  x.in <- x
  y.in <- y
  x <- x[!drop,]
  y <- y[!drop,]

  # Partial out fixed effects
  xr <- lm(x ~ as.factor(f))$residuals
  yr <- lm(y ~ as.factor(f))$residuals

  XX <- crossprod(xr)
  Xy <- crossprod(xr,y)

  pen <- lambdaCalculationPanel(penalty = penalty, y = yr, x = xr, f=f)
  lambda <- pen$lambda
  Ups0 <- Ups1 <- pen$Ups0
  lambda0 <- pen$lambda0

  mm <- 1
  s0 <- sqrt(var(y))

  while (mm <= control$numIter) {
    # calculation parameters
    if (mm==1 && post) {
      coefTemp <- hdm:::LassoShooting.fit(xr, yr, lambda/2, XX = XX, Xy = Xy)$coefficients
    } else {
      coefTemp <- hdm:::LassoShooting.fit(xr, yr, lambda, XX = XX, Xy = Xy)$coefficients
    }
    coefTemp[is.na(coefTemp)] <- 0
    ind1 <- (abs(coefTemp) > 0)
    x1 <- as.matrix(xr[, ind1, drop = FALSE])
    if (dim(x1)[2] == 0) {
      intercept.value <- mean(yr)
      coef <- rep(0,p)
      names(coef) <- colnames(x) #names(coefTemp)

      est <- list(coefficients = coef, beta=rep(0,p), intercept=intercept.value, index = rep(FALSE, p),
                  lambda = lambda, lambda0 = lambda0, loadings = Ups0, residuals = yr -
                    mean(yr), sigma = var(yr), iter = mm, call = match.call(),
                  options = list(post = post,
                                 control = control))
      #if (model) {
      #  est$model <- x
      #} else {
      #  est$model <- NULL
      #}
      est$tss <- est$rss <- sum((yr - mean(yr))^2)
      est$dev <- yr - mean(yr)
      class(est) <- "rlasso"
      return(est)
    }

    # refinement variance estimation
    if (post) {
      reg <- lm(yr ~ -1 + x1)
      coefT <- coef(reg)
      coefT[is.na(coefT)] <- 0
      e1 <- yr - x1 %*% coefT
      coefTemp[ind1] <- coefT
    }
    if (!post) {
      e1 <- yr - x1 %*% coefTemp[ind1]
    }
    s1 <- sqrt(var(e1))

    lc <- lambdaCalculationPanel(penalty, y=e1, x=x,f=f)
    Ups1 <- lc$Ups0
    lambda <- lc$lambda

    mm <- mm + 1
    if (abs(s0 - s1) < control$tol) {
      break
    }
    s0 <- s1
  }

  if (dim(x1)[2] == 0) {
    coefTemp = NULL
    ind1 <- rep(0, p)
  }
  coefTemp <- as.vector(coefTemp)
  coefTemp[abs(coefTemp) < control$threshold] <- 0
  ind1 <- as.vector(ind1)
  coefTemp <- as.vector(as.vector(coefTemp))
  names(coefTemp) <- names(ind1) <- colnames(x)
  intercept.value <- NA

  beta <- coefTemp

  s1 <- sqrt(var(e1))
  est <- list(coefficients = beta, beta=coefTemp, intercept=intercept.value, index = ind1, lambda = lambda,
              lambda0 = lambda0, loadings = Ups1, residuals = as.vector(e1), sigma = s1,
              iter = mm, call = match.call(), options = list(post = post,
                                                             control = control, penalty = penalty
                                                             ))
  est$tss <- sum((yr - mean(yr))^2)
  est$rss <- sum(est$residuals^2)
  est$dev <- yr - mean(yr)
  class(est) <- "rlassoPanel"
  return(est)
}

#' Function for Calculation of the penalty parameter
#'
#' This function implements different methods for calculation of the penalization parameter \eqn{\lambda}. Further details can be found under \link{rlassoPanel}.
#'
#' @param penalty list with options for the calculation of the penalty.
#' \itemize{
#' \item{\code{c} and \code{gamma}}{ constants for the penalty with default \code{c=1.1} and \code{gamma=0.1}}
#' }
#' @param x matrix of regressor variables
#' @param y residual which is used for calculation of the variance or
#'   the data-dependent loadings
#' @param f Grouping factor, we "cluster" on this factor. That is, we
#'   allow for dependence within levels of the factor, and assume
#'   independence across levels
#' @return The functions returns a list with the penalty \code{lambda} which is the product of \code{lambda0} and \code{Ups0}. \code{Ups0}
#' denotes either the variance (\code{independent} case) or the data-dependent loadings for the regressors. \code{method} gives the selected method for the calculation.
#' @export
lambdaCalculationPanel <- function(penalty = list(c = 1.1, gamma = 0.1),
                              y = NULL, x = NULL, f=NULL) {
  if (!exists("c", where = penalty)) {
    penalty$c = 1.1
  }
  if (!exists("gamma", where = penalty)) {
    penalty$gamma = 0.1
  }

  p <- dim(x)[2]
  n <- dim(x)[1]
  lambda0 <- 2 * penalty$c * sqrt(n) * qnorm(1 - penalty$gamma/(2 * p))
  #Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
  df <- data.frame(y=y,f=f,x)
  xxee <- by(df, f, FUN=function(d) {
    x <- as.matrix(d[,3:ncol(d)])
    y <- as.matrix(d$y)
    diag(t(x) %*% (y %*% t(y)) %*% x)
  }, simplify=TRUE)
  xxee <- matrix(unlist(xxee),ncol=ncol(x),byrow=TRUE)
  Ups0 <- 1/sqrt(n)*sqrt(colSums(xxee))
  lambda <- lambda0 * Ups0

  return(list(lambda0 = lambda0, lambda = lambda, Ups0 = Ups0, method = penalty))
}
