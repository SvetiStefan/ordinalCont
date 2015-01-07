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

negloglik_glf <- function(par, v, x){
  beta <- par[1:3]
  M <- par[4]
  B <- par[5]
  T <- par[6]
  
  g <- g_glf(v, M, B, T)
  dg <- dg_glf(v, B, T)
  if (any(dg<=0)) return(Inf)
  xb <- x %*% beta
  return(-sum(log(dg) + g + xb -2*log(1+exp(g+xb))))
}

contOrdEst <- function(par, v, x){
  fit <- optim(par,negloglik_glf, v=v, x=x)
  ## compute QR-decomposition of x
  qx <- qr(x)
  
  ## compute (x’x)^(-1) x’y
  coef <- fit$par
  names(coef) <- names(par)
  
  ## degrees of freedom and standard deviation of residuals
  df <- nrow(x)-ncol(x)
  require(boot) ##### put this somewhere else
  fitted.values <- inv.logit(g_glf(v, coef[4], coef[5], coef[6]) + x%*%coef[1:3])
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
    est <- contOrdEst(par, v, x)
    coef <- est$coefficients
    est$fitted.values <- as.vector(inv.logit(g_glf(v, coef[4], coef[5], coef[6]) + x%*%coef[1:3]))
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

summary.contOrd <- function(object, ...)
{
  se <- sqrt(diag(object$vcov))
  tval <- coef(object) / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df=object$df))
  res <- list(call=object$call,
              coefficients=TAB)
  class(res) <- "summary.contOrd"
  res
}

print.summary.contOrd <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}


contOrd.formula <- function(formula, data=list(), ...)
{
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  v <- model.response(mf)
  est <- contOrd.default(x, v, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}

predict.contOrd <- function(object, newdata=NULL, ...)
{
  if(is.null(newdata))
    y <- fitted(object)
  else{
    if(!is.null(object$formula)){
      ## model has been fitted using formula interface
      x <- model.matrix(object$formula, newdata)
    }
    else{
      x <- newdata
    }
    y <- as.vector(x %*% coef(object))
  }
  y 
}

plot.contOrd <- function(object, ...)
{
  plot(resid(object), main='Residuals',ylab='')
  lines(c(-10,length(resid(object))+20),c(0,0))
}
