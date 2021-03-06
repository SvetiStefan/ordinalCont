% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/ocm.methods.R
\name{plot.ocm}
\alias{plot.ocm}
\title{Plot method for Continuous Ordinal Fits}
\usage{
\method{plot}{ocm}(x, CIs = c("no", "vcov", "rnd.x.bootstrap",
  "fix.x.bootstrap", "param.bootstrap"), R = 1000,
  main = "g function (95\% CIs)", xlab = "Continuous ordinal scale",
  ylab = "", CIcol = "lightblue", ...)
}
\arguments{
\item{x}{an object of class \code{ocm}}

\item{CIs}{method used for confidence bands for the g function. \code{"no"} = no CIS [default]; \code{"vcov"} = Wald;
\code{"rnd.x.bootstrap"} = random-x bootstrap; \code{"fix.x.bootstrap"} = bootstrap with fixed-x
resampling; \code{"param.bootstrap"} = parametric bootstrap}

\item{R}{the number of bootstrap replicates. Ignored if CIs=\code{"no"}}

\item{main}{title of the plot. Defauts to ``g function (95\% CIs)''}

\item{xlab}{label of the x axis. Defaults to ``Continuous ordinal scale''}

\item{ylab}{label of the \code{y} axis. Defaults to an emtpy string}

\item{CIcol}{color of the confidence interval bands. Defaults to ``lightblue''}

\item{...}{further arguments passed to or from other methods}
}
\description{
Plots the g function as fitted in an \code{ocm} call.
}
\details{
The fitted g function of an \code{ocm} object is plotted.
If \code{CIs} is not \code{"no"}, 95\% confidence bands are also plotted.
Confidence bands computed with any of the bootstrapping options are
obtained with simple percentiles.
}
\examples{
ANZ0001.ocm <- ANZ0001[ANZ0001$cycleno==0 | ANZ0001$cycleno==5,]
ANZ0001.ocm$cycleno[ANZ0001.ocm$cycleno==5] <- 1
fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
plot(fit.overall, CIs="vcov")
\dontrun{
plot(fit.overall, CIs="rnd.x.bootstrap", R=100)
plot(fit.overall, CIs="fix.x.bootstrap", R=100)
plot(fit.overall, CIs="param.bootstrap", R=100)
}
}
\author{
Maurizio Manuguerra, Gillian Heller
}
\seealso{
\code{\link{ocm}}
}
\keyword{plot}

