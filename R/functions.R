#' Continuous ordinal regression
#'
#' This function performs the continuous ordinal regression with logit link using the generalized logistic function as g function and without random effects.
#' @param formula a formula expression as for regression models, of the form response ~ predictors. Only fixed effects are supported. The model must have an intercept: attempts to remove one will lead to a warning and will be ignored (TODO).
#' @param data  an optional data frame in which to interpret the variables occurring in the formulas.
#' @param start initial values for the parameters in the format c(alpha, beta, zeta), where alpha are the threshold parameters (adjusted for potential nominal effects), beta are the regression parameters and zeta are the scale parameters. (CHANGETHIS)
#' @param control a list of control parameters passed on to clm.control.
#' @param link link function, i.e., the type of location-scale distribution assumed for the latent distribution. The default "logit" link gives the proportional odds model.
#' @param ... additional arguments are passed on to clm.control.
#' @keywords likelihood, log-likelihood, ordinal regression.
#' @export
#' @examples
#' ocm() #write something here using a data set included with the package


ocm <- function(formula, data, start=NULL, control=list(), link = c("logit"), ...)
{
    if (any(sapply(attributes(terms(formula))$term.labels,function(x)grepl("|", x, fixed=T)))) stop("Random effects not yet supported.")
    mf <- model.frame(formula=formula, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    v <- model.response(mf)
    
    x <- as.matrix(x)
    v <- as.numeric(v)
    if (is.null(start)) beta_start <- set.beta_start(x,v)
    len_beta = length(beta_start)
    glf_start <- set.glf_start(x,v)
    start <- c(beta_start, glf_start)
    est <- ocmEst(start, v, x)
    coef <- est$coefficients
    beta <- coef[1:len_beta]
    par_g <- coef[(len_beta+1):(len_beta+2)]
    est$fitted.values <- as.vector(inv.logit(g_glf(v, par_g) + x%*%beta))
    est$residuals <- v - est$fitted.values
    est$call <- match.call()
    class(est) <- "ocm"
    est$formula <- formula
    est
}  
  
  

#' Continuous ordinal regression
#'
#' This function performs the comtinuous ordinal regression with logt link using the generalized logistic function as g function and without random effects.
#' @param formula A formula object (fixed effects).
#' @keywords likelihood, log-likelihood.
#' @export
#' @examples
#' print()

print.ocm <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients, ...)
}

#' @title Summarizing Continuous Ordinal Fits
#' @description summary method for class "ocm"
#'
#' @usage ## S3 method for class 'ocm'
#' summary(object, correlation = FALSE, symbolic.cor = FALSE, ...)
#' ## S3 method for class 'summary.ocm'
#' print(x, digits = max(3, getOption("digits") - 3),
#'      symbolic.cor = x$symbolic.cor,
#'      signif.stars = getOption("show.signif.stars"), ...)

#' @param object An object of class "ocm", usually, a result of a call to ocm.
#' @param x An object of class "summary.ocm", usually, a result of a call to summary.ocm.
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#' @keywords summary
#' @export
#' @examples
#' summary.ocm()

summary.ocm <- function(object, ...)
{
  se <- sqrt(diag(object$vcov))
  tval <- coef(object)[1:length(se)] / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df=object$df))
  res <- list(call=object$call,
              coefficients=TAB)
  class(res) <- "summary.ocm"
  res
}

#' @title Summarizing Continuous Ordinal Fits
#' @description summary method for class "ocm"
#'
#' @usage ## S3 method for class 'ocm'
#' summary(object, correlation = FALSE, symbolic.cor = FALSE, ...)
#' ## S3 method for class 'summary.ocm'
#' print(x, digits = max(3, getOption("digits") - 3),
#'      symbolic.cor = x$symbolic.cor,
#'      signif.stars = getOption("show.signif.stars"), ...)

#' @param object An object of class "ocm", usually, a result of a call to ocm.
#' @param x An object of class "summary.ocm", usually, a result of a call to summary.ocm.
#' @param digits The number of significant digits to use when printing.
#' @param ... Further arguments passed to or from other methods.
#' @keywords summary
#' @export
#' @examples
#' print.summary.ocm()

print.summary.ocm <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE, ...)
}




#' @title Predict method for Continuous Ordinal Fits
#' 
#' @description Predicted values based on ocm object.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()

predict.ocm <- function(object, newdata=NULL, ...)
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

plot.ocm <- function(object, ...)
{
  plot(resid(object), main='Residuals',ylab='')
  lines(c(-10,length(resid(object))+20),c(0,0))
}
