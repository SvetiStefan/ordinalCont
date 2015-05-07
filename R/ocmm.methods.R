#' Print continuous ordinal regression objects
#'
#' This function prints an ocmm object 
#' @param x An object of class \code{"ocmm"}, usually, a result of a call to ocm.
#' @param ... Further arguments passed to or from other methods.
#' @examples
#'\dontrun{
#' fit.overall.rnd  <- ocmm(overall  ~ cycleno + age + bsa + treatment + (1|randno), data=ANZ0001)
#' print(fit.overall.rnd)
#' }
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
#' @description summary method for class ``ocmm"
#' @param object An object of class ``ocmm", usually, a result of a call to ocm.
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
#' @description summary method for class ``summary.ocmm"
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
  print(x$coefficients_rnd, row.names=F)
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
#' @param CIs Indicates if confidence bands for the g function should be computed (based on the Wald 95\% CIs). \code{"no"} = no CIS [default]; \code{"vcov"} = Wald
#' @param R The number of bootstrap replicates. 
#' @param main the title of the plot. Defauts to ``g function (95% CIs)"
## #' @param xlab the label of the \code{x} axis. Defaults to ``Continuous ordinal scale" 
#' @param ylab the label of the \code{y} axis. Defaults to an emtpy string
#' @param CIcol the color of the confidence interval bands. Defaults to `lightblue'
#' @param ... Further arguments passed to or from other methods.
#' @details The fitted g function of an \code{ocmm} object is plotted. 
#' @keywords plot
#' @examples
#' \dontrun{
#' fit.overall.rnd  <- ocmm(overall  ~ cycleno + age + bsa + treatment + (1|randno), data=ANZ0001)
#' plot(fit.overall.rnd, CIs="vcov", R=100)
#' }
#' @export

plot.ocmm <- function(x, CIs = c('no','vcov'), R = 1000, main="g function (95% CIs)", xlab="Continuous ordinal scale", ylab="", CIcol='lightblue', ...)
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
  #CIs
  if (CIs != 'no'){
    #lines(v, ci_low, lty = 2)
    #lines(v, ci_high, lty = 2)
    polygon(c(v, rev(v)),c(ci_low,rev(ci_high)), col = CIcol)
    lines(v,gfun) #to superimpose gfun estimate on shaded area
    #if (CIs=='vcov' | CIs=='rnd.x.bootstrap' | CIs=='fix.x.bootstrap') lines(v, ci_median, lty = 2)
  }
  lines(c(.5,.5), ylim, col='grey')
  lines(xlim, c(0, 0), col='grey')
}

#' @title Anova method for Continuous Ordinal Fits
#' 
#' @description Comparison of continuous ordinal models in likelihood ratio tests.
#' @param object An ocm object.
#' @param ... one or more additional ocm objects.
#' @keywords anova
#' @return The method returns an object of class \code{anova.ocmm} and \code{data.frame}, reporting for each model, in hierarchical order:
#' \itemize{
#'   \item no.par the number of parameters
#'   \item AIC the Akaike information criterion
#'   \item loglik the log-likelihood
#'   \item LR.stat the likelihood ratio statistic
#'   \item df the difference in the degrees of freedom in the models being compared
#'   \item Pr(>Chisq) the p-value from the likelihood ratio test 
#' }
#' @export
#' @examples
#' \dontrun{
#' fit.overall.rnd  <- ocmm(overall  ~ cycleno + bsa + treatment + (1|randno), data=ANZ0001)
#' anova(fit.overall.rnd, update(fit.overall, .~. + age))
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
    stop('anova is not implemented for a single "ocmm" object')
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
#' @param x An object of class ``anova.ocmm".
#' @param digits controls the number of digits to print. Defaults to the maximum between the value returned by (getOption("digits") - 2) and 3
#' @param signif.stars a logical. Should the significance stars be printed? Defaults to the value returned by getOption("show.signif.stars")
#' @param ... Further arguments passed to or from other methods.
#' @return Prints \code{anova.ocmm} object
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

#' @title Calculate Variance-Covariance Matrix for a Fitted Model Object
#' @param object An \code{ocmm} object.
#' @param ... Further arguments to be passed to methods
#' @details For the generalized logistic g-function, the variance-covariance matrix of model parameters will be 
#' of dimension (\code{len_beta} +3)x(\code{len_beta} +3), where \code{len_beta}  is the number of 
#' beta coefficients in the model.
#' @export
#' @method vcov ocmm
#'  @return Variance-covariance matrix of model parameters
#' @seealso \code{\link{ocmm}}

vcov.ocmm <- function(object, ...) {
  vcov.ocm(object)
}


#' @title Extract the Log-Likelihood
#' @param object an \code{ocmm} object.
#' @param ... Further arguments to be passed to methods.
#' @usage \\method{logLik}{ocmm}(object, ...)
#' @method logLik ocmm
#' @seealso \code{\link{ocmm}}
#' @return Returns the log-likelihood of an \code{ocmm} object
#' @export


logLik.ocmm <- function(object, ...) {
  structure(object$logLik, df = object$df, nobs=object$nobs, class = "logLik.ocmm")
}


#' @title Extract AIC from a fitted Continuous Ordinal Model
#' @param fit \code{ocmm} object
#' @param scale parameter currently not used. For compatibility with general extractAIC method.
#' @param k  `weight' of the equivalent degrees of freedom (=: edf) 
#'  in the AIC formula. Defaults to 2.
#' @param ... further arguments (currently unused)
#' @details The generalised AIC is computed:
#' \deqn{-2\ell +k\cdot edf}
#' where \eqn{\ell} is the log likelihood, k=2 gives the AIC, and 
#' k=log(n) gives the BIC.
#' @seealso \code{\link{ocmm}}
#' @return Generalised AIC of \code{ocmm} object \code{fit}
#' @references  Akaike, H (1983). 
#' Information measures and model selection, 
#' \emph{Bulletin of the International Statistical Institute}, 50:277-290.
#' @export
#' @method extractAIC ocmm

extractAIC.ocmm <- function(fit, scale = 0, k = 2, ...) {
  edf <- fit$df
  c(edf, -2*fit$logLik + k * edf)
}


