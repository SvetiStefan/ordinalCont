#' Generalized logistic g function
#'
#' This function compute a parametric version of the g function following Richards (1959): g(v) = M + 1/B * log(Tv^T/(1-v^T))
#' @param v Vector of standarized scores from the continuous ordinal scale.
#' @param par Vector of M, the offset of the curve; B, the slope of the curve, and T, the symmetry of the curve.
#' @keywords Richards, generalized logistic function.
#' @export
#' @examples
#' g_glf()

g_glf <- function(v, par){
  #par = c(B, T) #par = c(M, B, T)
  #just  comment
  #return(par[1]+log(par[3]*v^par[3]/(1-v^par[3]))/par[2])
  return(log(par[2]*v^par[2]/(1-v^par[2]))/par[1])
}


#' Derivative of generalized logistic g function
#'
#' This function compute a parametric version of the g function following Richards (1959): dg(v)/dv = T/B * 1/[v*(1-v^T)]
#' @param v Vector of standarized scores from the continuous ordinal scale.
#' @param par Vector of B, the slope of the curve, and T, the symmetry of the curve.
#' @keywords Richards, derivative, generalized logistic function.
#' @export
#' @examples
#' dg_glf()

dg_glf <- function(v, par){
  #par = c(B, T)
  return(par[2]/par[1]/v/(1-v^par[2]))
}


#' Log-likelihood function for single subject using the generalized logistic function as g function and without random effects
#'
#' This function compute the Log-likelihood function for single subject using the generalized logistic function as g function and without random effects.
#' @param x Design matrix (fixed effects).
#' @param beta Regression coefficients.
#' @param v Vector of standarized scores from the continuous ordinal scale.
#' @param M The offset of the curve.
#' @param B The slope of the curve.
#' @param T The symmetry of the curve.
#' @keywords likelihood, log-likelihood.
#' @export
#' @examples
#' negloglik_glf()

negloglik_glf <- function(par, v, d.matrix, len_beta){
  x <- d.matrix
  beta <- par[1:len_beta]
  #par_g <- par[(len_beta+1):(len_beta+3)]
  par_g <- par[(len_beta+1):(len_beta+2)]
  par_dg <- par[(len_beta+1):(len_beta+2)]
  g <- g_glf(v, par_g)
  dg <- dg_glf(v, par_dg)
  if (any(dg<=0)) return(Inf)
  xb <- x %*% beta
  return(-sum(log(dg) + g + xb -2*log(1+exp(g+xb))))
}

contOrdEst <- function(start, v, x){
  require(numDeriv)
  require(boot)
  len_beta <- ncol(x)
  fit <- optim(par=start,negloglik_glf, v=v, d.matrix=x, len_beta=len_beta, method="BFGS", hessian = F)
  ## compute QR-decomposition of x
  #qx <- qr(x)
  #Hessian
  #H=fit$hessian
  H=hessian(negloglik_glf,fit$par,v=v, d.matrix=x,len_beta=len_beta)
  qrH <- qr(H)
  if(qrH$rank < nrow(H))
    stop("Cannot compute vcov: \nHessian is numerically singular")
  vcov <- solve.qr(qrH)
  
  ## compute (x’x)^(-1) x’y
  coef <- fit$par
  names(coef) <- names(start)
  len_beta = ncol(x)
  beta <- coef[1:len_beta]
  par_g <- coef[(len_beta+1):(len_beta+2)]
  
  ## degrees of freedom and standard deviation of residuals
  df <- nrow(x)-ncol(x)
  fitted.values <- inv.logit(g_glf(v, par_g) + x%*%beta)
  sigma2 <- sum((v - fitted.values)^2)/df
  
  ## compute sigma^2 * (x’x)^-1
  #vcov <- sigma2 * chol2inv(qx$qr)
  colnames(vcov) <- rownames(vcov) <- c(colnames(x),"B","T")
  list(coefficients = coef,
       vcov = vcov,
       sigma = sqrt(sigma2),
       df = df)
}

#' Continuous ordinal regression
#'
#' This function performs the comtinuous ordinal regression with logt link using the generalized logistic function as g function and without random effects.
#' @param formula A formula object (fixed effects).
#' @keywords likelihood, log-likelihood.
#' @export
#' @examples
#' contOrd()


contOrd <- function(x, ...) UseMethod("contOrd")

#' Continuous ordinal regression
#'
#' This function performs the comtinuous ordinal regression with logt link using the generalized logistic function as g function and without random effects.
#' @param formula A formula object (fixed effects).
#' @keywords likelihood, log-likelihood.
#' @export
#' @examples
#' contOrd.default()

contOrd.default <- function(x, v, start=NULL, ...)
{
    x <- as.matrix(x)
    v <- as.numeric(v)
    if (is.null(start)) beta_start <- set.beta_start(x,v)
    len_beta = length(beta_start)
    glf_start <- set.glf_start(x,v)
    start <- c(beta_start, glf_start)
    est <- contOrdEst(start, v, x)
    coef <- est$coefficients
    beta <- coef[1:len_beta]
    par_g <- coef[(len_beta+1):(len_beta+2)]
    est$fitted.values <- as.vector(inv.logit(g_glf(v, par_g) + x%*%beta))
    est$residuals <- v - est$fitted.values
    est$call <- match.call()
    class(est) <- "contOrd"
    est 
}

set.beta_start <- function(x,v){
  as.numeric(-coef(glm(v~0+x,family=gaussian(link="logit"))))
}

set.glf_start <- function(x,v){
  c(1,1)
}

#' Continuous ordinal regression
#'
#' This function performs the comtinuous ordinal regression with logt link using the generalized logistic function as g function and without random effects.
#' @param formula A formula object (fixed effects).
#' @keywords likelihood, log-likelihood.
#' @export
#' @examples
#' print()

print.contOrd <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients, ...)
}

#' @title Summarizing Continuous Ordinal Fits
#' @description summary method for class "contOrd"
#'
#' @usage ## S3 method for class 'contOrd'
#' summary(object, correlation = FALSE, symbolic.cor = FALSE, ...)
#' ## S3 method for class 'summary.contOrd'
#' print(x, digits = max(3, getOption("digits") - 3),
#'      symbolic.cor = x$symbolic.cor,
#'      signif.stars = getOption("show.signif.stars"), ...)

#' @param object An object of class "contOrd", usually, a result of a call to contOrd.
#' @param x An object of class "summary.contOrd", usually, a result of a call to summary.contOrd.
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#' @keywords summary
#' @export
#' @examples
#' summary.contOrd()

summary.contOrd <- function(object, ...)
{
  se <- sqrt(diag(object$vcov))
  tval <- coef(object)[1:length(se)] / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df=object$df))
  res <- list(call=object$call,
              coefficients=TAB)
  class(res) <- "summary.contOrd"
  res
}

#' @title Summarizing Continuous Ordinal Fits
#' @description summary method for class "contOrd"
#'
#' @usage ## S3 method for class 'contOrd'
#' summary(object, correlation = FALSE, symbolic.cor = FALSE, ...)
#' ## S3 method for class 'summary.contOrd'
#' print(x, digits = max(3, getOption("digits") - 3),
#'      symbolic.cor = x$symbolic.cor,
#'      signif.stars = getOption("show.signif.stars"), ...)

#' @param object An object of class "contOrd", usually, a result of a call to contOrd.
#' @param x An object of class "summary.contOrd", usually, a result of a call to summary.contOrd.
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#' @keywords summary
#' @export
#' @examples
#' print.summary.contOrd()

print.summary.contOrd <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE, ...)
}

#' Continuous ordinal regression
#'
#' This function performs the comtinuous ordinal regression with logt link using the generalized logistic function as g function and without random effects.
#' @param formula A formula object (fixed effects).
#' @keywords likelihood, log-likelihood.
#' @export
#' @examples
#' contOrd.formula()


contOrd.formula <- function(formula, data=list(), ...)
{
  if (any(sapply(attributes(terms(formula))$term.labels,function(x)grepl("|", x, fixed=T)))) stop("Random effects not yet supported.")
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  v <- model.response(mf)
  est <- contOrd.default(x, v, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' @title Predict method for Continuous Ordinal Fits
#' 
#' @description Predicted values based on contOrd object.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()

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
