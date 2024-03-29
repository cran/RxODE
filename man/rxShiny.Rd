% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rxShiny.R
\name{rxShiny}
\alias{rxShiny}
\alias{rxShiny.rxSolve}
\alias{rxShiny.default}
\title{Use Shiny to help develop an RxODE model}
\usage{
rxShiny(
  object,
  params = NULL,
  events = NULL,
  inits = NULL,
  ...,
  data = data.frame()
)

\method{rxShiny}{rxSolve}(
  object,
  params = NULL,
  events = NULL,
  inits = NULL,
  ...,
  data = data.frame()
)

\method{rxShiny}{default}(
  object = NULL,
  params = NULL,
  events = NULL,
  inits = NULL,
  ...,
  data = data.frame()
)
}
\arguments{
\item{object}{A RxODE family of objects. If not supplied a
2-compartment indirect effect model is used.  If it is
supplied, use the model associated with the RxODE object for
the model exploration.}

\item{params}{Initial parameters for model}

\item{events}{Event information (currently ignored)}

\item{inits}{Initial estimates for model}

\item{...}{Other arguments passed to rxShiny.  Currently doesn't
do anything.}

\item{data}{Any data that you would like to plot.  If the data has
a \code{time} variable as well as a compartment or calculated
variable that matches the RxODE model, the data will be added
to the plot of a specific compartment or calculated variable.}
}
\value{
Nothing; Starts a shiny server
}
\description{
Use Shiny to help develop an RxODE model
}
\author{
Zufar Mulyukov and Matthew L. Fidler
}
