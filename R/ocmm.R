#' Ordinal regression for continuous scales - with random effects 
#' Gillian: Please do not work on this as the code will be merged with ocm when fully working. 
#'
#' This function performs the continuous ordinal regression with logit link using the generalized logistic function as g function and without random effects.
#' @param formula a formula expression as for regression models, of the form response ~ predictors. Only fixed effects are supported. The model must have an intercept: attempts to remove one will lead to a warning and will be ignored (TODO).
#' @param data  an optional data frame in which to interpret the variables occurring in the formulas.
#' @param start initial values for the parameters in the format c(alpha, beta, zeta), where alpha are the threshold parameters (adjusted for potential nominal effects), beta are the regression parameters and zeta are the scale parameters. (CHANGETHIS)
#' @param control a list of control parameters passed on to clm.control.
#' @param link link function, i.e., the type of location-scale distribution assumed for the latent distribution. The default "logit" link gives the proportional odds model.
#' @param gfun A smooth monotonic function capable of capturing the non-linear nature of the ordinal measure. It defaults to the generalized logistic function, which is currently the only possibility.
#' @param quad A string indicating the type of quadrature used to integrate over the random effects. Can take values "Laplace" (Adaptive Gauss-Hermite quadrature using Laplace approximation; the default) or "GH" (Gauss-Hermite quadrature).
#' @param n_nodes Order of Gauss-Hermite rule used (number of nodes). 
#' @param ... additional arguments are passed on to clmm.control.
#' @keywords likelihood, log-likelihood, ordinal regression.
#' @export
#' @examples
#' # Change data set
#' #fitLaplace <- ocmm(vas ~ lasert1+lasert2+lasert3+ (1|ID), data=pain, quad="Laplace")
#' #fitGH <- ocmm(vas ~ lasert1+lasert2+lasert3+ (1|ID), data=pain, quad="GH") 
#' fit.overall.rnd  <- ocmm(overall  ~ cycleno + age + bsa + treatment + (1|randno), data=ANZ0001)
#' fit.phys.rnd     <- ocmm(phys 	   ~ cycleno + age + bsa + treatment + (1|randno), data=ANZ0001)
#' fit.pain.rnd 	  <- ocmm(pain 	   ~ cycleno + age + bsa + treatment + (1|randno), data=ANZ0001)
#' fit.mood.rnd 	  <- ocmm(mood 	   ~ cycleno + age + bsa + treatment + (1|randno), data=ANZ0001)
#' fit.nausvom.rnd  <- ocmm(nausvom	 ~ cycleno + age + bsa + treatment + (1|randno), data=ANZ0001)
#' fit.appetite.rnd <- ocmm(appetite ~ cycleno + age + bsa + treatment + (1|randno), data=ANZ0001)
#' summary(fit.overall.rnd)
#' summary(fit.phys.rnd)
#' summary(fit.pain.rnd)
#' summary(fit.mood.rnd)
#' summary(fit.nausvom.rnd)
#' summary(fit.appetite.rnd)


ocmm <- function(formula, data, weights, start=NULL, control=list(), link = c("logit"), gfun = c("glf"), quad=c("Laplace","GH"), n_nodes=10, ...)
{
  if (missing(formula)) 
    stop("Model needs a formula")
  terms <- attributes(terms(formula))$term.labels
  #print(attributes(terms(formula)))
  i_rnd <- which(sapply(attributes(terms(formula))$term.labels,function(x)grepl("|", x, fixed=T)))
  if (length(i_rnd) > 1) stop("Only one random effect supported. Please respecify your model.")
  if (length(i_rnd) == 0) stop("No random effects specified. Please call ocm.")
  #for (rndterms in terms[i_rnd]){ #not needed as only one rnd effect is supported.
  rndterms <- terms[i_rnd]
  fixterms <- terms[-i_rnd]
  ibar <- regexpr("|",rndterms,fixed=T)
  left <- trim(substr(rndterms,1,ibar-1))
  right <- trim(substr(rndterms,ibar+1, nchar(rndterms)))
  #cat("\nLeft:",as.numeric(left)," - Right:",right,"\n")
  if (as.numeric(left)!=1) stop("Only random effects on the intercept are supported in this version of ordinalCont.")
  if (any(sapply(c(":","*","|"), function(x)grepl(x, right,fixed=T)))) stop("Syntax incorrect or feature not implemented.")
  #cat("\nGoing to stratify by",right,"..\n")
  
  link <- match.arg(link)
  gfun <- match.arg(gfun)
  quad = match.arg(quad)
  if(missing(weights)) weights <- rep(1, nrow(data))
  keep <- weights > 0
  data <- data[keep,]
  weights <- weights[keep]
  
  to_stratify <- as.factor(data[,right])
  iclusters <- lapply(levels(to_stratify),function(x)which(to_stratify==x))
  
  mf <- model.frame(formula=formula, data=data)
  x.complete <- model.matrix(attr(mf, "terms"), data=mf)
  v <- model.response(mf)
  
  x <- as.matrix(x.complete)[,-c(1,i_rnd+1)] # 1 for the intercept, +1 is to consider intercept
  z <- as.matrix(x.complete)[,(i_rnd+1)] # +1 is to consider intercept
  v <- as.numeric(v)
  if (is.null(start)) {
    beta_start <- set.beta_start(x,v)
    len_beta = length(beta_start)
    if (gfun == 'glf') {
      gfun_start <- set.glf_start(x,v)
    }
    start <- c(beta_start, gfun_start, 1) #1 is the variance of the single rnd effect
    len_gfun <- length(gfun_start)
  }
  est <- ocmmEst(start, v, x, z, weights, link, gfun, rnd=right, n_nodes=n_nodes, quad=quad, iclusters)
  coef <- est$coefficients
  beta <- coef[1:len_beta]
  par_g <- coef[(len_beta+1):(len_beta+len_gfun)]
  #sigma_rnd <- coef[len_beta+len_gfun+1]
  est$len_beta <- len_beta
  est$len_gfun <- len_gfun
  est$len_rnd <- 1
  est$rnd <- right
  #fitted.values are the cumulative probabilities gamma(v|x)
  #est$fitted.values <- as.vector(inv.logit(g_glf(v, par_g) + x%*%beta))
  #The intercept is the parameter M of the glf. Fitted values v* = g^{-1}(W*) = g^{-1}(-x'B), where W=W*+epsilon
  #On the ordinal scale [0,1] (v* and v-v*):
  #est$fitted.values <- as.vector(g_glf_inv(-x%*%beta, par_g))
  #est$residuals <- v - est$fitted.values
  #On the latent scale (W* and epsilon=W-W*):
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
  class(est) <- "ocmm"
  est
  
}




negloglik_glf_rnd <- function(par, v, d.matrix, rnd.matrix, wts, len_beta, rnd, n_nodes, quad, iclusters)
  { #b is the random effect
  #densities_by_cluster <- NULL
  #for (ind in iclusters) densities_by_cluster <- c(densities_by_cluster, negloglik_glf_rnd2(ind, par=par, b=b, v=v, d.matrix=d.matrix, rnd.matrix=rnd.matrix, len_beta=len_beta, n_nodes=n_nodes, quad=quad))
  densities_by_cluster <- sapply(iclusters, negloglik_glf_rnd2, par=par, b=b, v=v, d.matrix=d.matrix, rnd.matrix=rnd.matrix, wts=wts, len_beta=len_beta, n_nodes=n_nodes, quad=quad)
  #if (length(densities_by_cluster) != length(wts)) {
  #  print(paste(length(densities_by_cluster), length(wts)))
  #  print( summary(densities_by_cluster))
  #  print(par)
  #}
  #loglik = -sum(wts * log(densities_by_cluster))
  loglik = -sum(log(densities_by_cluster))
  return(loglik)
}

negloglik_glf_rnd2 <- function(indices, par, b, v, d.matrix, rnd.matrix, wts, len_beta, n_nodes, quad){
  require(fastGHQuad)
  rule <- gaussHermiteData(n_nodes)
  x <- d.matrix[indices,]
  y <- v[indices]
  w <- wts[indices]
  beta <- par[1:len_beta]
  par_g <- par[(len_beta+1):(len_beta+3)]
  par_dg <- par[(len_beta+2):(len_beta+3)]
  sigma_rnd <- par[len_beta+4]
  g <- g_glf(y, par_g)
  dg <- dg_glf(y, par_dg)
  if (any(dg<=0)) return(Inf)
  xb <- as.numeric(x %*% beta)
  if (quad=="Laplace")
    return(aghQuad(g=density_glf_Laplace, muHat=0, sigmaHat=1, rule=rule, all_pars=cbind(g,dg,xb), sigma_rnd=sigma_rnd, w=w))
  else if (quad=="GH")
    return(ghQuad(f=density_glf_GH, rule=rule, all_pars=cbind(g,dg,xb), sigma_rnd=sigma_rnd, w=w))
}

density_glf_Laplace <- function(b, all_pars, sigma_rnd, w){
  lik_cluster <- apply(all_pars, 1, function(x, b, sigma_rnd){x[2]*exp(x[1]+x[3]+sqrt(2)*sigma_rnd*b)/(1+exp(x[1]+x[3]+sqrt(2)*sigma_rnd*b))^2/sqrt(pi)*exp(-b^2)},b=b,sigma_rnd=sigma_rnd)
  lik_cluster <- apply(lik_cluster,1,function(x,w)return(prod(x^w)),w=w)
  #lik_cluster <- apply(lik_cluster,1,prod)
  return(lik_cluster)
}

density_glf_GH <- function(b, all_pars, sigma_rnd, w){
  lik_cluster <- apply(all_pars, 1, function(x, b, sigma_rnd){x[2]*exp(x[1]+x[3]+sqrt(2)*sigma_rnd*b)/(1+exp(x[1]+x[3]+sqrt(2)*sigma_rnd*b))^2/sqrt(pi)},b=b,sigma_rnd=sigma_rnd)
  lik_cluster <- apply(lik_cluster,1,function(x,w)return(prod(x^w)),w=w)
  #lik_cluster <- apply(lik_cluster,1,prod)
  return(lik_cluster)
}

ocmmEst <- function(start, v, x, z, weights, link, gfun, rnd=NULL, n_nodes, quad, iclusters){
  len_beta <- ncol(x)
#  if (!is.null(z)){
#    cat("Random effects over",paste(rnd,collapse=', '),'\n') #ready for multiple rnd effects
#    stop("Work in progress...")
#  }
  if (gfun == "glf") {
    if (link == "logit"){
      fit <- optim(par=start,negloglik_glf_rnd, v=v, d.matrix=x, rnd.matrix=z, wts=weights, len_beta=len_beta, rnd=rnd, n_nodes=n_nodes, quad=quad, iclusters=iclusters, method="BFGS", hessian = T)
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
  sigma_rnd <- coef[len_beta+4]
  
  ## degrees of freedom and standard deviation of residuals
  df <- nrow(x)-ncol(x)-length(par_g)
  fitted.values <- inv.logit(g_glf(v, par_g) + x%*%beta)
  sigma2 <- sum((v - fitted.values)^2)/df
  
  ## compute sigma^2 * (x’x)^-1
  #vcov <- sigma2 * chol2inv(qx$qr)
  colnames(vcov) <- rownames(vcov) <- c(colnames(x),"M", "B", "T", "sigma_rnd")
  list(coefficients = coef,
       vcov = vcov,
       sigma = sqrt(sigma2),
       sigma_rnd = sigma_rnd,
       df = df,
       logLik = -fit$value)
}
