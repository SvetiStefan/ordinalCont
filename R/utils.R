#' Generalized logistic g function
#'
#' This function computes a parametric version of the g function following Richards (1959): 
#' \deqn{g(v) = M + \frac{1}{B} \log\left(\frac{Tv^T}{1-v^T}\right)}
#'  M is omitted as an intercept is always fitted.
#' @param v vector of standardized scores from the continuous ordinal scale, 0<v<1.
#' @param par vector of 2 elements: B, the slope of the curve, and T, the symmetry of the curve. 
#' @keywords Richards, generalized logistic function.
#' @details The generalized logistic functions maps from (0,1) to \eqn{(-\infty,\infty)}. 
#' B is the slope of the curve,  T is the symmetry and M is the offset. 
#' M is absorbed into the intercept of the model and so is set to zero in this function.
#' @return A vector of length equal to the length of v, with values g(v).
#' 
#' @references Richards, F. (1959). A flexible growth function for empirical use, 
#' \emph{Journal of Experimental Botany}, 10, 290–301.
#


g_glf <- function(v, par){
  return(log(par[2]*v^par[2]/(1-v^par[2]))/par[1])
}


#' Derivative of generalized logistic g function
#'
#' This function computes the derivative of the generalized logistic function as in Richards (1959): 
#' \deqn{\frac{dg(v)}{dv} = \frac{T}{B}  \frac{1}{v(1-v^{T})}}
#' @param v vector of standardized scores from the continuous ordinal scale, 0<v<1.
#' @param par vector of 2 elements: B, the slope of the curve, and T, the symmetry of the curve.
#' @keywords Richards, derivative, generalized logistic function.
#' @return A vector of length equal to the length of v, with values g'(v).
#' @references Richards, F. (1959). A flexible growth function for empirical use, 
#' \emph{Journal of Experimental Botany}, 10, 290–301.
#

dg_glf <- function(v, par){
  #par = c(B, T)
  return(par[2]/par[1]/v/(1-v^par[2]))
}


#' Inverse of generalized logistic g function
#'
#' This function computes the inverse of a parametric version of the g function 
#' following Richards (1959): 
#' \deqn{g^{-1}(W) = \left( \frac{e^{BW}}{T+e^{BW}}  \right)^{\frac{1}{T}}} M is omitted 
#' as it is absorbed into the 
#' model intercept.
#' @param W vector of scores on the latent scale \eqn{(-\infty,\infty)}. M, the 
#' offset of the g function, is 
#' estimated as the intercept of the model. 
#' \deqn{W=-x'\beta [without intercept] = -x'\beta - M} [as the model is 
#' estimated with an intercept].
#' @param par vector of 2 elements: B, the slope of the curve, and T, the symmetry 
#' of the curve. 
#' @return A vector of length equal to the length of W, with values \eqn{g^{-1}(W)}
#' @references Richards, F. (1959). A flexible growth function for empirical use, 
#' \emph{Journal of Experimental Botany}, 10, 290–301.
#' @keywords Richards, generalized logistic function.

g_glf_inv <- function(W, par){
  exp.part <- exp(par[1]*W)
  return((exp.part/(par[2]+exp.part))^(1/par[2]))
}

####################


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

####################

negloglik_glf_rnd <- function(par, v, d.matrix, rnd.matrix, len_beta, rnd, n_nodes, quad, iclusters)
  { #b is the random effect
  densities_by_cluster <- sapply(iclusters, negloglik_glf_rnd2, par=par, b=b, v=v, d.matrix=d.matrix, 
                                 rnd.matrix=rnd.matrix, len_beta=len_beta, n_nodes=n_nodes, quad=quad)
  loglik = -sum(log(densities_by_cluster))
  return(loglik)
  #return(-prod(densities_by_cluster))
}

negloglik_glf_rnd2 <- function(indices, par, b, v, d.matrix, rnd.matrix, len_beta, n_nodes, quad){
  require(fastGHQuad)
  rule <- gaussHermiteData(n_nodes)
  x <- d.matrix[indices,]
  y <- v[indices]
  beta <- par[1:len_beta]
  par_g <- par[(len_beta+1):(len_beta+2)]
  par_dg <- par[(len_beta+1):(len_beta+2)]
  sigma_rnd <- par[len_beta+3]
  g <- g_glf(y, par_g)
  dg <- dg_glf(y, par_dg)
  if (any(dg<=0)) return(Inf)
  xb <- as.numeric(x %*% beta)
  if (quad=="Laplace")
    return(aghQuad(g=density_glf_Laplace, muHat=0, sigmaHat=1, rule=rule, all_pars=cbind(g,dg,xb), 
                   sigma_rnd=sigma_rnd))
  else if (quad=="GH")
    return(ghQuad(f=density_glf_GH, rule=rule, all_pars=cbind(g,dg,xb), sigma_rnd=sigma_rnd))
}

density_glf_Laplace <- function(b, all_pars, sigma_rnd){
  lik_cluster <- apply(all_pars, 1, function(x, b, sigma_rnd){x[2]*exp(x[1]+x[3]+sqrt(2)*sigma_rnd*b)/(1+exp(x[1]+x[3]+sqrt(2)*sigma_rnd*b))^2/sqrt(pi)*exp(-b^2)},b=b,sigma_rnd=sigma_rnd)
  lik_cluster <- apply(lik_cluster,1,prod)
  return(lik_cluster)
}

density_glf_GH <- function(b, all_pars, sigma_rnd){
  lik_cluster <- apply(all_pars, 1, function(x, b, sigma_rnd){x[2]*exp(x[1]+x[3]+sqrt(2)*sigma_rnd*b)/(1+exp(x[1]+x[3]+sqrt(2)*sigma_rnd*b))^2/sqrt(pi)},b=b,sigma_rnd=sigma_rnd)
  lik_cluster <- apply(lik_cluster,1,prod)
  return(lik_cluster)
}

ocmmEst <- function(start, v, x, z, link, gfun, rnd=NULL, n_nodes, quad, iclusters){
  len_beta <- ncol(x)
#  if (!is.null(z)){
#    cat("Random effects over",paste(rnd,collapse=', '),'\n') #ready for multiple rnd effects
#    stop("Work in progress...")
#  }
  if (gfun == "glf") {
    if (link == "logit"){
      fit <- optim(par=start,negloglik_glf_rnd, v=v, d.matrix=x, rnd.matrix=z, len_beta=len_beta, rnd=rnd, n_nodes=n_nodes, quad=quad, iclusters=iclusters, method="BFGS", hessian = T)
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
  sigma_rnd <- coef[len_beta+3]
  
  ## degrees of freedom and standard deviation of residuals
  df <- nrow(x)-ncol(x)-length(par_g)
  fitted.values <- inv.logit(g_glf(v, par_g) + x%*%beta)
  sigma2 <- sum((v - fitted.values)^2)/df
  
  ## compute sigma^2 * (x’x)^-1
  #vcov <- sigma2 * chol2inv(qx$qr)
  colnames(vcov) <- rownames(vcov) <- c(colnames(x),"B","T","sigma_rnd")
  list(coefficients = coef,
       vcov = vcov,
       sigma = sqrt(sigma2),
       sigma_rnd = sigma_rnd,
       df = df,
       logLik = -fit$value)
}

####################



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

#Functions used in plot.ocm to bootstrapping data (random-x or fixed-x resampling) and find CIs.
rnd.x.bootstrap <- function(data, indices, fit){
  data <- data[indices,]
  mod <- update(fit, .~., data = data)
  coefficients(mod)
}
fix.x.bootstrap <- function(data, indices, fit){
  WminusM = as.numeric(-fit$x %*% fit$coefficients[1:fit$len_beta])
  data$new_v <- g_glf_inv(WminusM + residuals(fit)[indices], tail(fit$coefficients,2))
  mod <- update(fit, new_v ~., data = data)
  coefficients(mod)
}
param.bootstrap <- function(data, indices, fit){
  WminusM = as.numeric(-fit$x %*% fit$coefficients[1:fit$len_beta])
  data$new_v <- g_glf_inv(WminusM + rlogis(fit$sample.size), tail(fit$coefficients,2))
  mod <- update(fit, new_v ~., data = data)
  coefficients(mod)
}


# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)