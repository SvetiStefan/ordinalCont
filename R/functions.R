#' Generalized logistic g function
#'
#' This function compute a parametric version of the g function following Richards (1959): g(v) = M + 1/B * log(Tv^T/(1-v^T))
#' @param v Vector of standarized scores from the continuous ordinal scale
#' @param M The offset of the curve
#' @param B The slope of the curve
#' @param T The symmetry of the curve
#' @keywords Richards, generalized logistic function
#' @export
#' @examples
#' g_glf()

g_glf <- function(v, M, B, T){
  return(M+log(T*v^T/(1-v^T))/B)
}


#' Derivative of generalized logistic g function
#'
#' This function compute a parametric version of the g function following Richards (1959): dg(v)/dv = T/B * 1/[v*(1-v^T)]
#' @param v Vector of standarized scores from the continuous ordinal scale
#' @param B The slope of the curve
#' @param T The symmetry of the curve
#' @keywords Richards, derivative, generalized logistic function
#' @export
#' @examples
#' dg_glf()

dg_glf <- function(v, B, T){
  return(T/B/v/(1-v^T))
}


#' Log-likelihood function for single subject using the generalized logistic function as g function and without random effects
#'
#' This function compute the Log-likelihood function for single subject using the generalized logistic function as g function and without random effects.
#' @param x Design matrix (fixed effects)
#' @param beta Regression coefficients
#' @param v Vector of standarized scores from the continuous ordinal scale
#' @param M The offset of the curve
#' @param B The slope of the curve
#' @param T The symmetry of the curve
#' @keywords likelihood, log-likelihood
#' @export
#' @examples
#' loglik_glf()

negloglik_glf <- function(pars, v, x){
  beta <- pars[1:3]
  M <- pars[4]
  B <- pars[5]
  T <- pars[6]
  
  g <- g_glf(v, M, B, T)
  dg <- dg_glf(v, B, T)
  if (any(dg<=0)) return(Inf)
  xb <- x %*% beta
  return(sum(log(dg) + g + xb -2*log(1+exp(g+xb))))
}

contOrdEst <- function(pars, v, x){
  fit <- optim(pars,negloglik_glf, v=v, x=x)
  
  ## compute QR-decomposition of x
  qx <- qr(x)
  
  ## compute (x’x)^(-1) x’y
  coef <- fit$par
  
  ## degrees of freedom and standard deviation of residuals
  df <- nrow(x)-ncol(x)
  fitted.values <- g_glf(v, coef[4], coef[5], coef[6]) + x%*%coef[1:3]
  sigma2 <- sum((v - fitted.values)^2)/df
  
  ## compute sigma^2 * (x’x)^-1
  vcov <- sigma2 * chol2inv(qx$qr)
  colnames(vcov) <- rownames(vcov) <- colnames(x)
  list(coefficients = coef,
       vcov = vcov,
       sigma = sqrt(sigma2),
       df = df)
}



contOrd <- function(x, ...) UseMethod("contOrd")

contOrd.default <- function(x, v, par, ...)
{
    x <- as.matrix(x)
    v <- as.numeric(v)
    est <- contOrdEst(pars, v, x)
    coef <- est$coefficients
    est$fitted.values <- as.vector(g_glf(v, coef[4], coef[5], coef[6]) + x%*%coef[1:3])
    est$residuals <- v - est$fitted.values
    est$call <- match.call()
    class(est) <- "contOrd"
    est 
}

print.contOrd <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}
