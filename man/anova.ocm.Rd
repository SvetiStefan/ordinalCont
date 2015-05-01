% Please edit documentation in R/ocm.methods.R
\name{anova.ocm}
\alias{anova.ocm}
\title{Anova method for Continuous Ordinal Fits}
\usage{
\method{anova}{ocm}(object, ...)
}
\arguments{
\item{object}{an object of class \code{ocm}}

\item{...}{one or more additional \code{ocm} objects.}
}
\value{
object of class \code{anova.ocm}, consisting of anova table giving AIC,
 log likelihood, likelihood ratio statistic, df and p-value
 for each model, in hierarchical order
}
\description{
Comparison of continuous ordinal models using likelihood ratio tests.
}
\details{
Likelihood ratio testing of nested models is performed.
}
\author{
Maurizio Manuguerra
}
\seealso{
\code{\link{ocm}}, \code{\link{print.anova.ocm}}
}
\keyword{anova}

