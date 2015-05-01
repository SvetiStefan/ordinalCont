#' Ordinal regression for continuous scales
#'
#' This function performs continuous ordinal regression with logit link using the 
#' generalized logistic function as g function. Random effects are not supported.
#' @param formula a formula expression as for regression models, of the form 
#' response ~ predictors. Only fixed effects are supported. 
#' The model must have an intercept: attempts to remove one will lead to a warning and will be 
#' ignored (TODO).
#' @param data  an optional data frame in which to interpret the variables occurring in the 
#' formulas.
#' @param start initial values for the parameters in the format c(alpha, beta, zeta), where 
#' alpha are the threshold parameters (adjusted for potential nominal effects), beta are the 
#' regression parameters and zeta are the scale parameters. (CHANGETHIS)
#' @param control a list of control parameters passed on to clm.control.
#' @param link link function, i.e. the type of location-scale distribution assumed for the latent 
#' distribution. The default "logit" link gives the proportional odds model.
#' @param gfun A smooth monotonic function capable of capturing the non-linear nature of the 
#' ordinal measure. It defaults to the generalized logistic function, which is currently the only 
#' possibility.
#' @param ... additional arguments are passed on to clm.control.
#' @keywords likelihood, log-likelihood, ordinal regression.
#' @details Ordinal regression analysis is a convenient tool for analyzing ordinal response variables 
#' in the presence of covariates. We extend this methodology to the case of continuous self-rating 
#' scales such as the Visual Analog Scale (VAS) used in pain assessment, or the Linear Analog 
#' Self-Assessment (LASA) scales in quality of life studies. Subjects are
#' typically given a linear scale of 100 mm and asked to put a mark where they perceive
#' themselves. These scales  measure subjects' 
#' perception of an intangible quantity, and cannot be handled as ratio variables because of their 
#' inherent nonlinearity.  We express  the likelihood in terms of a function (the "g function")
#'  connecting the  
#' scale with an underlying continuous latent  variable. In the current version the g function 
#' is taken as 
#' the generalized logistic function (Richards 1959). This has 3 parameters: 
#'  \code{M}, the offset, \code{B}, the slope, and \code{T}, the symmetry of the curve.
#' The link function is the inverse of the CDF of the assumed underlying distribution of the 
#' latent variable. Currently 
#' the logit link, which corresponds to a standard logistic distribution, is implemented. 
#' (This implies a proportional odds model.)
#' A regression framework supporting fixed effects
#'  is implemented. The likelihood is maximized using XXXXX MAURIZIO which function?.
#' 
#' @seealso \code{\link{ocmm}}
#' @return an object of type \code{ocm}. Parameter estimates are in \code{coefficients}. 
#' The last 3 elements of \code{coefficients} are the parameters of the g function: 
#' \code{M},  \code{B},  and \code{T}.
#'  @references Manuguerra M, Heller GZ (2010). Ordinal Regression Models for Continuous 
#'  Scales, \emph{The International Journal of Biostatistics}: 6(1), Article 14.
#'@references Richards, F. (1959). A flexible growth function for empirical use, 
#' \emph{Journal of Experimental Botany}, 10, 290-301.
#' @author Maurizio Manuguerra, Gillian Heller


#' @export
#' @examples
#' ANZ0001.ocm <- ANZ0001[ANZ0001$cycleno==0 | ANZ0001$cycleno==5,]
#' ANZ0001.ocm$cycleno[ANZ0001.ocm$cycleno==5] <- 1
#' fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' fit.phys 	  <- ocm(phys 	  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' fit.pain 	  <- ocm(pain 	  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' fit.mood 	  <- ocm(mood 	  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' fit.nausvom  <- ocm(nausvom  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' fit.appetite <- ocm(appetite ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' summary(fit.overall)
#' summary(fit.phys)
#' summary(fit.pain)
#' summary(fit.mood)
#' summary(fit.nausvom)
#' summary(fit.appetite)
#' par(mfrow=c(2,3))
#' plot(fit.overall, CIs='vcov', R=100)
#' plot(fit.phys, CIs='vcov', R=100)
#' plot(fit.pain, CIs='vcov', R=100)
#' plot(fit.mood, CIs='vcov', R=100)
#' plot(fit.nausvom, CIs='vcov', R=100)
#' plot(fit.appetite, CIs='vcov', R=100)
#' par(mfrow=c(1,1))


ocm <- function(formula, data, weights, start=NULL, control=list(), link = c("logit"), 
                gfun = c("glf"), ...)
{
  #FIXME check for the intercept in formula.
  if (any(sapply(attributes(terms(formula))$term.labels,function(x)grepl("|", x, fixed=T)))) 
    stop("Random effects specified. Please call ocmm.")
  if (missing(formula)) 
    stop("Model needs a formula")
  #formula = update(formula,.~.-1) ##no_intercept --> problems with contrasts, 
  #better to drop the intercept later
  link <- match.arg(link)
  gfun <- match.arg(gfun)
  if(missing(weights)) weights <- rep(1, nrow(data))
  keep <- weights > 0
  data <- data[keep,]
  weights <- weights[keep]
    
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  lenx <- ncol(x)
  v <- model.response(mf)
  xnames <- dimnames(x)[[2]][2:lenx]
  x <- as.matrix(x)[,2:lenx]
  v <- as.numeric(v)
  if (is.null(start)) {
    beta_start <- set.beta_start(x,v)
    len_beta = length(beta_start)
    names(beta_start) <- xnames[1:len_beta]
    if (gfun == 'glf') {
      gfun_start <- set.glf_start(x,v)
      names(gfun_start) <- c("M", "B", "T")
    }
    start <- c(beta_start, gfun_start)
    len_gfun <- length(gfun_start)
  }
  est <- ocmEst(start, v, x, weights, link, gfun)
  coef <- est$coefficients
  beta <- coef[1:len_beta]
  par_g <- coef[(len_beta+1):(len_beta+len_gfun)]
  est$len_beta <- len_beta
  est$len_gfun <- len_gfun
  est$fitted.values <- as.vector(-x%*%beta)
  est$residuals <- g_glf(v, par_g) - est$fitted.values
  est$v <- v
  est$x <- x
  est$sample.size <- nrow(x)
  est$nobs <- sum(weights)
  est$call <- match.call()
  est$no.pars <- length(coef)
  est$data <- data
  est$link <- link
  est$gfun <- gfun
  est$formula <- formula
  class(est) <- "ocm"
  est
}  

#' @title Log-likelihood function for the fixed-effects model, using the generalized logistic 
#' function as g function and the logit link function
#'
#' @details This function computes minus the log-likelihood function for a fixed-effects model using the 
#' generalized logistic function as g function and the logit link function.
#' @param par vector of regression coefficients (first \code{len_beta} elements), 
#' and \code{M},  \code{B}, \code{T}, (offset, slope and symmetry of the g function - 
#' last 3 elements)
#' @param v vector of standardized scores from the continuous ordinal scale
#' @param d.matrix design matrix (fixed effects)
#' @param wts ????
#' @param len_beta length of the regression coefficients vector
#' @keywords likelihood, log-likelihood.
#' @return minus the log-likelihood at parameter values \code{par} 

negloglik_glf <- function(par, v, d.matrix, wts, len_beta){
  return(-sum(wts * logdensity_glf(par, v, d.matrix, len_beta)))
}

logdensity_glf <- function(par, v, d.matrix, len_beta){
  x <- d.matrix
  beta <- par[1:len_beta]
  par_g <- par[(len_beta+1):(len_beta+3)]
  par_dg <- par[(len_beta+2):(len_beta+3)]
  g <- g_glf(v, par_g)
  dg <- dg_glf(v, par_dg)
  if (any(dg<=0)) return(Inf)
  xb <- x %*% beta
  return(log(dg) + g + xb -2*log(1+exp(g+xb)))
}

ocmEst <- function(start, v, x, weights, link, gfun){
  len_beta <- ncol(x)
  if (gfun == "glf") {
    if (link == "logit"){
      fit <- optim(par=start,negloglik_glf, v=v, d.matrix=x, wts=weights, len_beta=len_beta, 
                   method="BFGS", hessian = T)
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
  par_g <- coef[(len_beta+1):(len_beta+3)]
  
  ## degrees of freedom and standard deviation of residuals
  df <- nrow(x)-ncol(x)-length(par_g)
  fitted.values <- inv.logit(g_glf(v, par_g) + x%*%beta)
  sigma2 <- sum((v - fitted.values)^2)/df
  
  ## compute sigma^2 * (x’x)^-1
  #vcov <- sigma2 * chol2inv(qx$qr)
  colnames(vcov) <- rownames(vcov) <- c(colnames(x),"M", "B", "T")
  list(coefficients = coef,
       vcov = vcov,
       sigma = sqrt(sigma2),
       df = df,
       logLik = -fit$value)
}
