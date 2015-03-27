#' Print continuous ordinal regression objects
#'
#' This function prints an ocmm object 
#' @param x An object of class "ocm", usually, a result of a call to ocm.
#' @param ... Further arguments passed to or from other methods.
#' @keywords likelihood, log-likelihood.
#' @method print ocmm
#' @export

print.ocmm <- function(x, ...)
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
#' @method summary ocmm
#' @keywords summary
#' @export

summary.ocmm <- function(object, ...)
{
  se <- sqrt(diag(object$vcov))
  tval <- coef(object)[1:length(se)] / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df=object$df))
  TABrnd <- data.frame(Groups = object$rnd,
                  Name = "Intercept",
                  #FIXME make general
                  Variance = round(object$sigma_rnd^2,3),
                  Std.Dev. = round(object$sigma_rnd,3))
  res <- list(call=object$call,
              coefficients=TAB,
              coefficients_rnd=TABrnd,
              len_beta=object$len_beta,
              len_gfun=object$len_gfun,
              len_rnd=object$len_rnd,
              rnd=object$rnd)
  class(res) <- "summary.ocmm"
  print(res, ...)
}

#' @title Summarizing Continuous Ordinal Fits
#' @description summary method for class "summary.ocm"
#' @param x An object of class "summary.ocmm", usually, a result of a call to summary.ocmm.
#' @param ... Further arguments passed to or from other methods.
#' @keywords summary
#' @export

print.summary.ocmm <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Random effects:\n")
  #printCoefmat(x$coefficients_rnd, P.values = FALSE, has.Pvalue = FALSE, signif.legend = FALSE, ...)
  #FIXME make general and good looking
  #cat(names(x$coefficients_rnd),"\n")
  print(x$coefficients_rnd)
  cat("\n")
  cat("Coefficients:\n")
  printCoefmat(x$coefficients[1:x$len_beta,], P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE, ...)
  cat("\n")
  cat("g function:\n")
  printCoefmat(x$coefficients[(x$len_beta+1):(x$len_beta+x$len_gfun),], P.values = TRUE, has.Pvalue = TRUE, ...)
}




#' @title Plot method for Continuous Ordinal Fits
#' 
#' @description This function plots the g function as fitted in an ocm call.
#' @param x An ocm object.
#' @param CIs Indicates if confidence bands for the g function should be computed (based on the Wald 95\% CIs).
#' @param R The number of bootstrap replicates. 
#' @param ... Further arguments passed to or from other methods.
#' @keywords plot
#' @export

plot.ocmm <- function(x, CIs = c('no','vcov'), R = 1000, main="g function (95% CIs)", xlab="Continuous ordinal scale", ylab="", ...)
{
  #FIXME: this works for glf only: make general?
  CIs <- match.arg(CIs)
  R <- as.integer(R)
  len_beta <- x$len_beta
  indices = c(len_beta+1, len_beta+2, len_beta+3)
  params_g <- coef(x)[indices]
  v <- seq(0.01, 0.99, by=0.01)
  gfun <- g_glf(v, params_g)
  xlim <- c(0,1)
  ylim <- c(min(gfun), max(gfun))
  if (CIs=='vcov') {
    #require(MASS)
    vcov_g <- x$vcov[indices, indices]
    #rparams <- mvrnorm(R, params_g, vcov_g, empirical=TRUE)
    rparams <- mvrnormR(R, params_g, vcov_g)
    #FIXME write efficiently
    all_gfuns <- NULL
    for (i in 1:R) all_gfuns <- rbind(all_gfuns, g_glf(v, rparams[i,]))
    ci_low  <- apply(all_gfuns, 2, function(x)quantile(x, 0.025))
    ci_high <- apply(all_gfuns, 2, function(x)quantile(x, 0.975)) 
    ylim <- c(min(ci_low), max(ci_high))
  }
  plot(v, gfun, main=main, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, t='l')
  lines(c(.5,.5), ylim, col='grey')
  lines(xlim, c(0, 0), col='grey')
  #CIs
  if (CIs != 'no'){
    lines(v, ci_low, lty = 2)
    lines(v, ci_high, lty = 2)
  }
}

#' @title Anova method for Continuous Ordinal Fits
#' 
#' @description Comparison of continuous ordinal models in likelihood ratio tests.
#' @param object An ocm object.
#' @param ... one or more additional ocm objects.
#' @keywords anova
#' @export
#' @examples
#' \dontrun{
#' fitLaplace = ocmm(vas ~ lasert1+lasert2+lasert3+ (1|ID), data=pain, quad="Laplace")
#' anova(fitLaplace, update(fitLaplace, . ~ . + localisa))
#' }


anova.ocmm <- function(object, ...)
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
  class(tab) <- c("anova.ocmm", "data.frame")
  tab
}


#' @title Print anova.ocm objects
#' 
#' @description Print the results of the comparison of continuous ordinal models in likelihood ratio tests.
#' @param x An object of class "anova.ocm".
#' @param ... Further arguments passed to or from other methods.
#' @keywords summary, anova
#' @export

print.anova.ocmm <-
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

#' @export

vcov.ocmm <- function(object, ...) vcov.ocm(object)

#' @export

nobs.ocmm <- function(object, ...) nobs.ocm(object)

#' @export

coef.ocmm <- function(object, ...) coef.ocm(object)

#' @export

logLik.ocmm <- function(object, ...) logLik.ocm(object)

#' @export

extractAIC.ocmm <- function(object, ...) extractAIC.ocm(object)