#' Ordinal regression for continuous scales
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
#' # Change data set
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
    est$len_beta <- len_beta
    #fitted.values are the cumulative probabilities gamma(v|x)
    #est$fitted.values <- as.vector(inv.logit(g_glf(v, par_g) + x%*%beta))
    #The intercept is the parameter M of the glf. Fitted values v* = g^{-1}(W*) = g^{-1}(-x'B), where W=W*+epsilon
    #On the ordinal scale [0,1] (v* and v-v*):
    #est$fitted.values <- as.vector(g_glf_inv(-x%*%beta, par_g))
    #est$residuals <- v - est$fitted.values
    #On the latent scale (W* and epsilon=W-W*):
    est$fitted.values <- as.vector(-x%*%beta+coef[1])
    est$residuals <- g_glf(v, par_g) - est$fitted.values
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
#' @param ... additional arguments are passed on to clm.control.
#' @keywords likelihood, log-likelihood, ordinal regression.
#' @export
#' @examples
#' # Change data set
#' fitLaplace = ocmm(vas ~ lasert1+lasert2+lasert3+ (1|ID), data=pain, quad="Laplace")
#' fitGH = ocmm(vas ~ lasert1+lasert2+lasert3+ (1|ID), data=pain, quad="GH") 

ocmm <- function(formula, data, start=NULL, control=list(), link = c("logit"), gfun = c("glf"), quad=c("Laplace","GH"), n_nodes=10, ...)
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
  cat("\nGoing to stratify by",right,"..\n")
  
  link <- match.arg(link)
  gfun <- match.arg(gfun)
  quad = match.arg(quad)
  
  to_stratify <- as.factor(data[,right])
  iclusters <- lapply(levels(to_stratify),function(x)which(to_stratify==x))
    
  mf <- model.frame(formula=formula, data=data)
  x.complete <- model.matrix(attr(mf, "terms"), data=mf)
  v <- model.response(mf)
  
  x <- as.matrix(x.complete)[,-(i_rnd+1)] #+1 for the intercept
  z <- as.matrix(x.complete)[,(i_rnd+1)] #+1 for the intercept
  v <- as.numeric(v)
  if (is.null(start)) beta_start <- set.beta_start(x,v)
  len_beta = length(beta_start)
  glf_start <- set.glf_start(x,v)
  start <- c(beta_start, glf_start, 1) #1 is the variance of the single rnd effect
  ### As ocm + z --- end ###
  est <- ocmmEst(start, v, x, z, link, gfun, rnd=right, n_nodes=n_nodes, quad=quad, iclusters)
  coef <- est$coefficients
  beta <- coef[1:len_beta]
  par_g <- coef[(len_beta+1):(len_beta+2)]
  sigma_rnd <- coef[len_beta+3]
  est$len_beta <- len_beta
  #fitted.values are the cumulative probabilities gamma(v|x)
  #est$fitted.values <- as.vector(inv.logit(g_glf(v, par_g) + x%*%beta))
  #The intercept is the parameter M of the glf. Fitted values v* = g^{-1}(W*) = g^{-1}(-x'B), where W=W*+epsilon
  #On the ordinal scale [0,1] (v* and v-v*):
  #est$fitted.values <- as.vector(g_glf_inv(-x%*%beta, par_g))
  #est$residuals <- v - est$fitted.values
  #On the latent scale (W* and epsilon=W-W*):
  est$fitted.values <- as.vector(-x%*%beta+coef[1])
  est$residuals <- g_glf(v, par_g) - est$fitted.values
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
  formula <- object$formula
  params <- coef(object)
  if(is.null(newdata)){
    x <- object$x 
  }else{
    x <- model.matrix(object$formula, newdata)
  }
  len_beta <- ncol(x)
  ndens <- 100
  v <- seq(0.01, 0.99, length.out = ndens)
  modes <- NULL
  densities <- NULL
  #FIXME: rewrite efficiently
  for (subject in 1:nrow(x)){
    d.matrix <- matrix(rep(x[subject,], ndens), nrow = ndens, dimnames = list(as.character(1:ndens), colnames(x)), byrow = TRUE)
    densities <- rbind(densities, t(logdensity_glf(par = params, v = v, d.matrix = d.matrix, len_beta = len_beta)))
    modes <- c(modes, v[which.max(logdensity_glf(par = params, v = v, d.matrix = d.matrix, len_beta = len_beta))])
  }
  #y = logdensity_glf(par = params, v = v, d.matrix = x, len_beta = len_beta)
  #plot(v,y)
  pred <- list(mode = modes, density = densities, x = v, formula = formula, newdata = newdata)
  class(pred) <- "predict.ocm"
  return(pred)
}

#' @title Print the output of predict method
#' @description print method for class "predict.ocm"
#' @param x An object of class "predict.ocm".
#' @param ... Further arguments passed to or from other methods.
#' @keywords predict
#' @export

print.predict.ocm <- function(x, ...)
{
  cat("\nThe data set use by the predict method contains",length(x$mode),"records.\n")
  cat("Call:\n")
  print(x$formula)
  cat("\nSummary of modes:\n")
  print(summary(x$mode), ...)
}

#' @title Plot the probability densities as from the output of the predict method
#' @description plot method for class "predict.ocm"
#' @param x An object of class "predict.ocm".
#' @param ... Further arguments passed to or from other methods.
#' @keywords predict, plot
#' @export

plot.predict.ocm <- function(x, ...)
{
  cat("Call:\n")
  print(x$formula)
  cat("The data set used in the predict methos contains ",nrow(x$density)," records.\n")
  #cat("Please press 'enter' to start/advance plotting and q to quit.\n")
  for (i in 1:nrow(x$density)){
    input <- readline(paste("Press 'enter' to plot the probability density of record ",i,", 'q' to quit: ",sep=''))
    if (input == "q") break()
    plot(x$x, exp(x$density[i,]), ylab="Probability Density", main=paste("Record", i), xlab=paste("mode =", round(x$mode[i],3)), t='l')
    lines(rep(x$mode[i],2), c(0, max(exp(x$density[i,]))), lty=21)
  }
}


#' @title Plot method for Continuous Ordinal Fits
#' 
#' @description This function plots the g function as fitted in an ocm call.
#' @param x An ocm object.
#' @param CIs Indicates if confidence bands for the g function should be computed based on the Wald 95\% CIs or by bootstrapping. In  the latter case, bootstrapping can be performed using a random-x or a fixed-x resampling. 95\% CIs computed with either of the bootstrapping options are obtained with simple percentiles. 
#' @param R The number of bootstrap replicates. 
#' @param ... Further arguments passed to or from other methods.
#' @keywords plot
#' @export

plot.ocm <- function(x, CIs = c('simple','rnd.x.bootstrap','fix.x.bootstrap','param.bootstrap'), R = 1000, ...)
{
  #FIXME: this works for glf only: make general?
  #FIXME: with bootstrapping, when a variable is a factor, it can go out of observation for some level making optim fail.
  CIs <- match.arg(CIs)
  R <- as.integer(R)
  M <- x$coefficients[1]
  params <- tail(coef(x), 2)
  v <- seq(0.01, 0.99, by=0.01)
  gfun <- M + g_glf(v, params)
  xlim <- c(0,1)
  ylim <- c(min(gfun), max(gfun))
  if (CIs=='simple') {
    #FIXME this is a very simple version, not the bootstrap one.
    sds <- sqrt(diag(x$vcov))
    sdM <- sds[1]
    sM <- rnorm(R, M, sdM)
    sdparams <- tail(sds, 2)
    sparams <- matrix(rnorm(2*R, params, sdparams), ncol = 2, byrow = T)
    all_gfuns <- NULL
    for (i in 1:R){
        all_gfuns <- rbind(all_gfuns, sM[i] + g_glf(v, sparams[i,]))
    }
    ci_low  <- apply(all_gfuns, 2, function(x)quantile(x, 0.025))
    ci_high <- apply(all_gfuns, 2, function(x)quantile(x, 0.975)) 
    ylim <- c(min(ci_low), max(ci_high))
  } else if (CIs=='rnd.x.bootstrap' | CIs=='fix.x.bootstrap'| CIs=='param.bootstrap'){
    require(boot)
    bs <- boot(x$data, eval(parse(text=CIs)), R, fit = x)
    all_gfuns <- NULL
    for (i in 1:R){
      all_gfuns <- rbind(all_gfuns, bs$t[i,1] + g_glf(v, tail(bs$t[i,],2)))
    }
    ci_low  <- apply(all_gfuns, 2, function(x)quantile(x, 0.025))
    ci_median <- apply(all_gfuns, 2, function(x)quantile(x, 0.5))
    ci_high <- apply(all_gfuns, 2, function(x)quantile(x, 0.975)) 
    ylim <- c(min(ci_low), max(ci_high))
  }
  plot(v, gfun, main='g function', xlim = xlim, ylim = ylim, xlab = 'Continuous ordinal scale', ylab = '', t='l')
  lines(c(.5,.5), ylim, col='grey')
  lines(xlim, c(0, 0), col='grey')
  #CIs
  lines(v, ci_low, lty = 2)
  lines(v, ci_high, lty = 2)
  if (CIs=='rnd.x.bootstrap' | CIs=='fix.x.bootstrap') lines(v, ci_median, lty = 2)
}

#' @title Anova method for Continuous Ordinal Fits
#' 
#' @description Comparison of continuous ordinal models in likelihood ratio tests.
#' @param object An ocm object.
#' @param ... one or more additional ocm objects.
#' @keywords anova
#' @export


anova.ocm <- function(object, ...)
  ### requires that ocm objects have components:
  ###  no.pars: no. parameters used
  ###  call$formula
  ###  link (character)
  ###  gfun (character)
  ###  logLik
  ###
{
  mc <- match.call()
  dots <- list(...)
  ## remove 'test' and 'type' arguments from dots-list:
  not.keep <- which(names(dots) %in% c("test", "type"))
  if(length(not.keep)) {
    message("'test' and 'type' arguments ignored in anova.ocm\n")
    dots <- dots[-not.keep]
  }
  if(length(dots) == 0)
    stop('anova is not implemented for a single "ocm" object')
  mlist <- c(list(object), dots)
  if(!all(sapply(mlist, function(model)
    inherits(model, c("ocm", "ocmm")))))
    stop("only 'ocm' and 'ocmm' objects are allowed")
  nfitted <- sapply(mlist, function(x) length(x$fitted.values))
  if(any(nfitted != nfitted[1L]))
    stop("models were not all fitted to the same dataset")
  no.par <- sapply(mlist, function(x) x$no.pars)
  ## order list with increasing no. par:
  ord <- order(no.par, decreasing=FALSE)
  mlist <- mlist[ord]
  no.par <- no.par[ord]
  no.tests <- length(mlist)
  ## extract formulas, links, gfun:
  forms <- sapply(mlist, function(x) deparse(x$call$formula))
  links <- sapply(mlist, function(x) x$link)
  gfun <- sapply(mlist, function(x) x$gfun)
  models <- data.frame(forms)
  models.names <- c('formula', "link", "gfun")
  models <- cbind(models, data.frame(links, gfun))
  ## extract AIC, logLik, statistics, df, p-values:
  AIC <- sapply(mlist, function(x) -2*x$logLik + 2*x$no.pars)
  logLiks <- sapply(mlist, function(x) x$logLik)
  statistic <- c(NA, 2*diff(sapply(mlist, function(x) x$logLik)))
  df <- c(NA, diff(no.par))
  pval <- c(NA, pchisq(statistic[-1], df[-1], lower.tail=FALSE))
  pval[!is.na(df) & df==0] <- NA
  ## collect results in data.frames:
  tab <- data.frame(no.par, AIC, logLiks, statistic, df, pval)
  tab.names <- c("no.par", "AIC", "logLik", "LR.stat", "df",
                 "Pr(>Chisq)")
  colnames(tab) <- tab.names
  #mnames <- sapply(as.list(mc), deparse)[-1]
  #rownames(tab) <- rownames(models) <- mnames[ord]
  rownames(tab) <- rownames(models) <- paste("Model ",1:length(mlist),":",sep='')
  colnames(models) <- models.names
  attr(tab, "models") <- models
  attr(tab, "heading") <-
    "Likelihood ratio tests of ordinal regression models for continuous scales:\n"
  class(tab) <- c("anova.ocm", "data.frame")
  tab
}

#' @export

#' @title Print anova.ocm objects
#' 
#' @description Print the results of the comparison of continuous ordinal models in likelihood ratio tests.
#' @param x An object of class "anova.ocm".
#' @param ... Further arguments passed to or from other methods.
#' @keywords summary, anova
#' @export

print.anova.ocm <-
  function(x, digits=max(getOption("digits") - 2, 3),
           signif.stars=getOption("show.signif.stars"), ...)
  {
    if (!is.null(heading <- attr(x, "heading")))
      cat(heading, "\n")
    models <- attr(x, "models")
    #row.names(models) <- paste("Model ",1:nrow(models),":",sep='')
    print(models, right=FALSE)
    cat("\n")
    printCoefmat(x, digits=digits, signif.stars=signif.stars,
                 tst.ind=4, cs.ind=NULL, # zap.ind=2, #c(1,5),
                 P.values=TRUE, has.Pvalue=TRUE, na.print="", ...)
    return(invisible(x))
  }
