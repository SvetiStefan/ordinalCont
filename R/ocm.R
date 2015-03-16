#' Ordinal regression for continuous scales
#'
#' This function performs  continuous ordinal regression with logit link using the 
#' generalized logistic function as g function and without random effects.
#' @param formula a formula expression as for regression models, of the form 
#' response ~ predictors. Only fixed effects are supported. 
#' The model must have an intercept: attempts to remove one will lead to a warning and will be ignored (TODO).
#' @param data  an optional data frame in which to interpret the variables occurring in the formulas.
#' @param start initial values for the parameters in the format c(alpha, beta, zeta), where alpha are the threshold parameters (adjusted for potential nominal effects), beta are the regression parameters and zeta are the scale parameters. (CHANGETHIS)
#' @param control a list of control parameters passed on to clm.control.
#' @param link link function, i.e. the type of location-scale distribution assumed for the latent distribution. The default "logit" link gives the proportional odds model.
#' @param gfun A smooth monotonic function capable of capturing the non-linear nature of the ordinal measure. It defaults to the generalized logistic function, which is currently the only possibility.
#' @param ... additional arguments are passed on to clm.control.
#' @keywords likelihood, log-likelihood, ordinal regression.
#' @return an object of type 'ocm'
#' @export
#' @examples
#' # Change data set
#' fit = ocm(vas ~ lasert1+lasert2+lasert3, data=pain)


ocm <- function(formula, data, start=NULL, control=list(), link = c("logit"), gfun = c("glf"), ...)
{
  if (any(sapply(attributes(terms(formula))$term.labels,function(x)grepl("|", x, fixed=T)))) 
    stop("Random effects specified. Please call ocmm.")
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
  est$len_beta <- len_beta
  #fitted.values are the cumulative probabilities gamma(v|x)
  #est$fitted.values <- as.vector(inv.logit(g_glf(v, par_g) + x%*%beta))
  #The intercept is the parameter M of the glf. Fitted values v* = g^{-1}(W*) = g^{-1}(-x'B), where W=W*+epsilon
  #On the ordinal scale [0,1] (v* and v-v*):
  #est$fitted.values <- as.vector(g_glf_inv(-x%*%beta, par_g))
  #est$residuals <- v - est$fitted.values
  #On the latent scale (W* and epsilon=W-W*):
  est$fitted.values <- as.vector(-x%*%beta+coef[1])
  est$residuals <- coef[1] + g_glf(v, par_g) - est$fitted.values
  est$v <- v
  est$x <- x
  est$sample.size <- nrow(x)
  est$call <- match.call()
  est$no.pars <- length(coef)
  est$data <- data
  est$link <- link
  est$gfun <- gfun
  est$formula <- formula
  class(est) <- "ocm"
  est
}  

#' Log-likelihood function for the fixed-effects model, using the generalized logistic function as 
#' g function and the logit link function.
#'
#' This function compute the log-likelihood function for a fixed-effects model using the 
#' generalized logistic function as g function and the logit link function.
#' @param par Vector of B, the slope of the curve, and T, the symmetry of the curve.
#' @param v Vector of standarized scores from the continuous ordinal scale.
#' @param d.matrix Design matrix (fixed effects).
#' @param len_beta Length of the regression coefficients vector.
#' @keywords likelihood, log-likelihood.

negloglik_glf <- function(par, v, d.matrix, len_beta){
  return(-sum(logdensity_glf(par, v, d.matrix, len_beta)))
}

logdensity_glf <- function(par, v, d.matrix, len_beta){
  x <- d.matrix
  beta <- par[1:len_beta]
  par_g <- par[(len_beta+1):(len_beta+2)]
  par_dg <- par[(len_beta+1):(len_beta+2)]
  g <- g_glf(v, par_g)
  dg <- dg_glf(v, par_dg)
  if (any(dg<=0)) return(Inf)
  xb <- x %*% beta
  return(log(dg) + g + xb -2*log(1+exp(g+xb)))
}

ocmEst <- function(start, v, x, link, gfun){
  len_beta <- ncol(x)
  if (gfun == "glf") {
    if (link == "logit"){
      fit <- optim(par=start,negloglik_glf, v=v, d.matrix=x, len_beta=len_beta, method="BFGS", hessian = T)
    } else {
      stop("link function not implemented.")
    }
  } else {
    stop("g function not implemented.")
  }
  ## compute QR-decomposition of x
  #qx <- qr(x)
  #Hessian
  H=fit$hessian
  #require(numDeriv)
  #H=hessian(negloglik_glf,fit$par,v=v, d.matrix=x,len_beta=len_beta)
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
  df <- nrow(x)-ncol(x)-length(par_g)
  fitted.values <- inv.logit(g_glf(v, par_g) + x%*%beta)
  sigma2 <- sum((v - fitted.values)^2)/df
  
  ## compute sigma^2 * (x’x)^-1
  #vcov <- sigma2 * chol2inv(qx$qr)
  colnames(vcov) <- rownames(vcov) <- c(colnames(x),"B","T")
  list(coefficients = coef,
       vcov = vcov,
       sigma = sqrt(sigma2),
       df = df,
       logLik = -fit$value)
}
