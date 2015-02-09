#' Continuous ordinal regression
#'
#' This function performs the continuous ordinal regression with logit link using the generalized logistic function as g function and without random effects.
#' @param formula a formula expression as for regression models, of the form response ~ predictors. Only fixed effects are supported. The model must have an intercept: attempts to remove one will lead to a warning and will be ignored (TODO).
#' @param data  an optional data frame in which to interpret the variables occurring in the formulas.
#' @param start initial values for the parameters in the format c(alpha, beta, zeta), where alpha are the threshold parameters (adjusted for potential nominal effects), beta are the regression parameters and zeta are the scale parameters. (CHANGETHIS)
#' @param control a list of control parameters passed on to clm.control.
#' @param link link function, i.e., the type of location-scale distribution assumed for the latent distribution. The default "logit" link gives the proportional odds model.
#' @param gfun A smooth monotonic function capable of capturing the non-linear nature of the ordinal measure. It defaults to the generalized logistic function, which is currently the only possibility.
#' @param ... additional arguments are passed on to clm.control.
#' @keywords likelihood, log-likelihood, ordinal regression.
#' @export
#' @examples
#' # Change this with something that uses the data set included with the package (Gillian?)
#' fit = ocm(vas ~ lasert1+lasert2+lasert3, data=pain) 


ocm <- function(formula, data, start=NULL, control=list(), link = c("logit"), gfun = c("glf"), ...)
{
    if (any(sapply(attributes(terms(formula))$term.labels,function(x)grepl("|", x, fixed=T)))) 
      stop("Random effects not yet supported.")
    if (missing(formula)) 
      stop("Model needs a formula")
    link <- match.arg(link)
    gfun <- match.arg(gfun)
    
    mf <- model.frame(formula=formula, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    v <- model.response(mf)
    
    x <- as.matrix(x)
    v <- as.numeric(v)
    if (is.null(start)) beta_start <- set.beta_start(x,v)
    len_beta = length(beta_start)
    glf_start <- set.glf_start(x,v)
    start <- c(beta_start, glf_start)
    est <- ocmEst(start, v, x, link, gfun)
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
#' @param x An object of class "ocm", usually, a result of a call to ocm.
#' @param ... Further arguments passed to or from other methods.
#' @keywords likelihood, log-likelihood.
#' @method print ocm
#' @export

print.ocm <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients, ...)
}

#' @title Summarizing Continuous Ordinal Fits
#' @description summary method for class "ocm"
#' @param object An object of class "ocm", usually, a result of a call to ocm.
#' @param ... Further arguments passed to or from other methods.
#' @method summary ocm
#' @keywords summary
#' @export

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
  print(res, ...)
}

#' @title Summarizing Continuous Ordinal Fits
#' @description summary method for class "summary.ocm"
#' @param x An object of class "summary.ocm", usually, a result of a call to summary.ocm.
#' @param ... Further arguments passed to or from other methods.
#' @keywords summary
#' @export

print.summary.ocm <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE, ...)
}




#' @title Predict method for Continuous Ordinal Fits
#' 
#' @description Predicted values based on ocm object.
#' @param object An ocm object.
#' @param newdata optionally, a data frame in which to look for variables with which to predict. Note that all predictor variables should be present having the same names as the variables used to fit the model.
#' @param ... Further arguments passed to or from other methods.
#' @keywords predict
#' @export

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

#' @title Plot method for Continuous Ordinal Fits
#' 
#' @description Plot based on ocm object.
#' @param x An ocm object.
#' @param ... Further arguments passed to or from other methods.
#' @keywords plot
#' @export

plot.ocm <- function(x, ...)
{
  plot(resid(x), main='Residuals',ylab='')
  lines(c(-10,length(resid(x))+20),c(0,0))
}
