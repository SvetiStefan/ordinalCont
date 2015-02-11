#' Generalized logistic g function
#'
#' This function computes a parametric version of the g function following Richards (1959): g(v) = M + 1/B * log(Tv^T/(1-v^T)). M is omitted as an intercept is always fitted.
#' @param v Vector of standarized scores from the continuous ordinal scale.
#' @param par Vector of M, the offset of the curve; B, the slope of the curve, and T, the symmetry of the curve.
#' @keywords Richards, generalized logistic function.

g_glf <- function(v, par){
  return(log(par[2]*v^par[2]/(1-v^par[2]))/par[1])
}


#' Derivative of generalized logistic g function
#'
#' This function compute a parametric version of the g function following Richards (1959): \deqn{\frac{dg(v)}{dv} = T/B  \frac{1}{v(1-v^{T})}}
#' @param v Vector of standarized scores from the continuous ordinal scale.
#' @param par Vector of B, the slope of the curve, and T, the symmetry of the curve.
#' @keywords Richards, derivative, generalized logistic function.

dg_glf <- function(v, par){
  #par = c(B, T)
  return(par[2]/par[1]/v/(1-v^par[2]))
}


#' Log-likelihood function for the fixed-effects model, using the generalized logistic function as g function and the logit link function.
#'
#' This function compute the log-likelihood function for a fixed-effects model using the generalized logistic function as g function and the logit link function.
#' @param par Vector of B, the slope of the curve, and T, the symmetry of the curve.
#' @param v Vector of standarized scores from the continuous ordinal scale.
#' @param d.matrix Design matrix (fixed effects).
#' @param len_beta Length of the regression coefficients vector.
#' @keywords likelihood, log-likelihood.

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


set.beta_start <- function(x,v){
  vv=ifelse(v<median(v),0,1)
  #as.numeric(-coef(glm(vv~0+x,family=gaussian(link="logit"))))
  as.numeric(-coef(glm(vv~0+x,family=binomial(link="logit"))))
}

set.glf_start <- function(x,v){
  c(1,1)
}


inv.logit <- function(x){
  ifelse(is.finite(x),exp(x)/(1+exp(x)),sign(x)*Inf)
}