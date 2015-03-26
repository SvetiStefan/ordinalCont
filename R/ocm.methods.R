#' Print continuous ordinal regression objects
#'
#' \code{print.ocm} is the ordinalCont specific method for the generic function \code{print}, 
#' which prints objects of class \code{'ocm'}.
#' @param x An object of class \code{'ocm'}, usually, a result of a call to \code{ocm}.
#' @param ... Further arguments passed to or from other methods.
#' @return Prints an \code{ocm} object
#' @keywords likelihood, log-likelihood.
#' @method print ocm
#' @seealso \code{\link{ocm}}, \code{\link{summary.ocm}}
#' @examples xxx add example
#' @export

print.ocm <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients, ...)
}

#' @title Summarizing Continuous Ordinal Fits
#' @description Summary method for class \code{'ocm'}
#' @param object An object of class \code{"ocm"}, usually, a result of a call to \code{ocm}.
#' @param ... Further arguments passed to or from other methods.
#' @method summary ocm
#' @keywords summary
#' @seealso \code{\link{ocm}}, \code{\link{print.ocm}}
#' @examples xxx add example
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
              coefficients=TAB,
              len_beta=object$len_beta,
              len_gfun=object$len_gfun)
  class(res) <- "summary.ocm"
  print(res, ...)
}

#' @title Summarizing Continuous Ordinal Fits
#' @description Summary method for class \code{"summary.ocm"}
#' @param x An object of class \code{"summary.ocm"}, usually a result of a call to \code{summary.ocm}.
#' @param ... further arguments passed to or from other methods.
#' @details The table of parameter estimates is printed.
#' @examples xxx add example
#' @keywords summary
#' @export

print.summary.ocm <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Coefficients:\n")
  printCoefmat(x$coefficients[1:x$len_beta,], P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE, ...)
  cat("\n")
  cat("g function:\n")
  printCoefmat(x$coefficients[(x$len_beta+1):(x$len_beta+x$len_gfun),], P.values = TRUE, has.Pvalue = TRUE, ...)
}




#' @title Predict method for Continuous Ordinal Fits
#' 
#' @description Predicted values based on \code{ocm} object.
#' @param object an object of class \code{"ocm"}, usually, a result of a call to \code{ocm}
#' @param newdata optionally, a data frame in which to look for variables with 
#' which to predict. 
#' Note that all predictor variables should be present, having the same names as the variables 
#' used to fit the model. If \code{NULL}, predictions are computed for the original dataset.
#' @param ... Further arguments passed to or from other methods.
#' @keywords predict
#' @return MAURIZIO we need to specify this (I'm not sure what you've done)
#' @details MAURIZIO we need to specify this (I'm not sure what you've done)
#' @examples xxx add example
#' @seealso \code{\link{ocm}}
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
#' @description Print method for class \code{"predict.ocm"}
#' @param x An object of class \code{"predict.ocm"}
#' @param ... Further arguments passed to or from other methods
#' @keywords predict
#' @details The table of predictions from \code{predict.ocm} is printed.
#' @seealso \code{\link{predict.ocm}},\code{\link{ocm}}
#' @examples xxx add example
#' @export

print.predict.ocm <- function(x, ...)
{
  cat("\nThe data set used by the predict method contains",length(x$mode),"records.\n")
  cat("Call:\n")
  print(update(x$formula, .~.+1))
  cat("\nSummary of modes:\n")
  print(summary(x$mode), ...)
}

#' @title Plot  probability densities  from  output of  predict method
#' @description Plot method for class \code{"predict.ocm"}
#' @param x An object of class \code{"predict.ocm"}
#' @param records An integer or a vector of integers. The number of the record/s 
#' in the data set for which the density has to be plotted. If not specified, the 
#' function will  plot all of them.
#' @param ... Further arguments passed to or from other methods.
#' @details The probability densities from \code{predict.ocm}  are plotted.
#' @seealso \code{\link{predict.ocm}},\code{\link{ocm}}
#' @examples xxx add example
#' @keywords predict, plot
#' @export

plot.predict.ocm <- function(x, records=NULL, ...)
{
  if (is.null(records)) records=1:nrow(x$density)
  cat("Call:\n")
  print(x$formula)
  cat("The data set used in the predict methos contains ",nrow(x$density)," records.\n")
  #cat("Please press 'enter' to start/advance plotting and q to quit.\n")
  for (i in records){
    input <- readline(paste("Press 'enter' to plot the probability density of record ",i,", 'q' to quit: ",sep=''))
    if (input == "q") break()
    plot(x$x, exp(x$density[i,]), ylab="Probability Density", main=paste("Record", i), 
         xlab=paste("mode =", round(x$mode[i],3)), t='l')
    lines(rep(x$mode[i],2), c(0, max(exp(x$density[i,]))), lty=21)
  }
}


#' @title Plot method for Continuous Ordinal Fits
#' 
#' @description Plots the g function as fitted in an \code{ocm} call.
#' @param x an object of class \code{ocm}
#' @param CIs method used for confidence bands for the g function. \code{"no"} = no CIS; \code{"vcov"} = Wald; 
#' \code{"rnd.x.bootstrap"} = random-x bootstrap; \code{"fix.x.bootstrap"} = bootstrap with fixed-x 
#' resampling; \code{"param.bootstrap"} = parametric bootstrap. 
#' MAURIZIO can you add CIs="NULL" for the option of no confidence bands? M: Done, but used "no" instead of NULL.
#' And make this the default? (In this case R=NULL would have to be the default.)
#' @param R the number of bootstrap replicates [ignored if CIs='no']
#' @param ... further arguments passed to or from other methods
#' @details The fitted g function of an \code{ocm} object is plotted. 
#' If \code{CIs} is not \code{"no"}, 95\% confidence bands are also plotted.
#' Confidence bands computed with any of the bootstrapping options are 
#' obtained with simple percentiles. 
#' @keywords plot
#' @export
#' @seealso \code{\link{ocm}}
#' @examples
#' fit <- ocm(vas ~ lasert1 + lasert2 + lasert3, data = pain)
#' plot(fit, CIs="vcov")

plot.ocm <- function(x, CIs = c('no', 'vcov','rnd.x.bootstrap','fix.x.bootstrap','param.bootstrap'), R = 1000, ...)
{
  #FIXME: this works for glf only: make general?
  #FIXME: with bootstrapping, when a variable is a factor, it can go out of observation for some level making optim fail.
  CIs <- match.arg(CIs)
  R <- as.integer(R)
  len_beta <- x$len_beta
  indices = c(len_beta+1, len_beta+2, len_beta+3)
  params_g <- coef(x)[indices]
  v <- seq(0.01, 0.99, by=0.01)
  gfun <- g_glf(v, params_g)
  xlim <- c(0,1)
  ylim <- c(min(gfun), max(gfun))
  if (CIs=='vcov'){
    #require(MASS)
    vcov_g <- x$vcov[indices, indices]
    #rparams <- mvrnorm(R, params_g, vcov_g, empirical=TRUE)
    rparams <- mvrnormR(R, params_g, vcov_g)
    #FIXME write efficiently
    all_gfuns <- NULL
    for (i in 1:R) all_gfuns <- rbind(all_gfuns, g_glf(v, rparams[i,]))
    ci_low  <- apply(all_gfuns, 2, function(x)quantile(x, 0.025))
    ci_median <- apply(all_gfuns, 2, function(x)quantile(x, 0.5))
    ci_high <- apply(all_gfuns, 2, function(x)quantile(x, 0.975)) 
    ylim <- c(min(ci_low), max(ci_high))
  } else if (CIs=='rnd.x.bootstrap' | CIs=='fix.x.bootstrap'| CIs=='param.bootstrap'){
    require(boot)
    bs <- boot(x$data, eval(parse(text=CIs)), R, fit = x)
    all_gfuns <- NULL
    for (i in 1:R){
      all_gfuns <- rbind(all_gfuns, g_glf(v, bs$t[i,indices]))
    }
    ci_low  <- apply(all_gfuns, 2, function(x)quantile(x, 0.025))
    ci_median <- apply(all_gfuns, 2, function(x)quantile(x, 0.5))
    ci_high <- apply(all_gfuns, 2, function(x)quantile(x, 0.975)) 
    ylim <- c(min(ci_low), max(ci_high))
  }
  plot(v, gfun, main='g function (95% CIs)', xlim = xlim, ylim = ylim, xlab = 'Continuous ordinal scale', ylab = '', t='l')
  lines(c(.5,.5), ylim, col='grey')
  lines(xlim, c(0, 0), col='grey')
  #CIs
  if (CIs != 'no'){
    lines(v, ci_low, lty = 2)
    lines(v, ci_high, lty = 2)
    if (CIs=='vcov' | CIs=='rnd.x.bootstrap' | CIs=='fix.x.bootstrap') lines(v, ci_median, lty = 2)
  }
}

#' @title Anova method for Continuous Ordinal Fits [GH up to here 19/3/15]
#' @description Comparison of continuous ordinal models in likelihood ratio tests. 
#' @param object an object of class \code{ocm}
#' @param ... one or more additional \code{ocm} objects.
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


#' @title Print anova.ocm objects
#' 
#' @description Print the results of the comparison of continuous ordinal models in likelihood ratio tests.
#' @param x An object of class "anova.ocm".
#' @param ... Further arguments passed to or from other methods.
#' @keywords summary, anova
#' @export

print.anova.ocm <- function(x, digits=max(getOption("digits") - 2, 3), signif.stars=getOption("show.signif.stars"), ...){
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

#' @export

coef.ocm <- function(x, ...){
    x$coefficients
  }

#' title Confidence intervals
#' @export
#' @examples
#' fit <- ocm(vas ~ lasert1 + lasert2 + lasert3, data = pain)
#' confint(fit)

confint.ocm <- function(object, parm, level = 0.95){
    stopifnot(is.numeric(level) && length(level) == 1 && level > 0 && level < 1)
    cf <- coef(object)    
    pnames <- names(cf)
    if (missing(parm)) 
        parm <- pnames
    else if (is.numeric(parm)) 
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    fac <- qt(a, object$df)
    pct <- format.perc(a, 3)
    #pct <- as.character(round(a*100,3))
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
        pct))
    ses <- sqrt(diag(object$vcov))[parm]
    ci[] <- cf[parm] + ses %o% fac
    ci
}

#' @export

logLik.ocm <- function(object, ...)
  structure(object$logLik, df = object$df, nobs=object$nobs,
            class = "logLik")

#' @export

extractAIC.ocm <- function(fit, scale = 0, k = 2, ...) {
  edf <- fit$df
  c(edf, -2*fit$logLik + k * edf)
}

#' @export

nobs.ocm <- function(object, ...) object$nobs

#' @export

model.frame.ocm <- function(object, ...) {
  if(is.null(mod <- object$data[,all.vars(object$formula)]))
    stop("Cannot extract model.frame.")
  else
    mod
}

#' @export

model.matrix.ocm <- function(object, ...) object$x

#' @export

terms.ocm <- function(object, ...) terms(object$formula)

#' @export

vcov.ocm <- function(object, ...) object$vcov
