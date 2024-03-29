% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RxODE.R
\name{coef.RxODE}
\alias{coef.RxODE}
\title{Return the RxODE coefficients}
\usage{
\method{coef}{RxODE}(object, ...)
}
\arguments{
\item{object}{is an RxODE object}

\item{...}{ignored arguments}
}
\value{
a rxCoef object with the following
\itemize{
\item \code{params}  is a list of strings for parameters for the RxODE object
\item \code{state} is a list of strings for the names of each state in
the RxODE object.
\item \code{ini} is the model specified default values for the
parameters.
\item \code{RxODE} is the referring RxODE object
}
}
\description{
This returns the parameters , state variables
}
\author{
Matthew L.Fidler
}
\keyword{internal}
