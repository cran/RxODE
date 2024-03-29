% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/et.R
\name{eventTable}
\alias{eventTable}
\title{Create an event table object}
\usage{
eventTable(amount.units = NA, time.units = NA)
}
\arguments{
\item{amount.units}{string denoting the amount dosing units, e.g.,
\dQuote{mg}, \dQuote{ug}. Default to \code{NA} to denote
unspecified units.  It could also be a solved RxODE object.  In
that case, eventTable(obj) returns the eventTable that was used
to solve the RxODE object.}

\item{time.units}{string denoting the time units, e.g.,
\dQuote{hours}, \dQuote{days}. Default to \code{"hours"}.

An \code{eventTable} is an object that consists of a data.frame
storing ordered time-stamped events of an (unspecified) PK/PD
dynamic system, units (strings) for dosing and time records, plus a
list of functions to add and extract event records.

Currently, events can be of two types: dosing events that represent
inputs to the system and sampling time events that represent
observations of the system with \sQuote{amount.units} and
\sQuote{time.units}, respectively.}
}
\value{
A modified data.frame with the following accessible functions:
\itemize{
\item \code{get.EventTable()} returns the current event table
\item \code{\link[=add.dosing]{add.dosing()}}  adds dosing records to the event table.
\item \code{get.dosing()} returns a data.frame of dosing records.
\item \code{clear.dosing()} clears or deletes all dosing from event table
\item `\code{\link[=add.sampling]{add.sampling()}} adds sampling time observation records to the
event table.
\item \code{get.sampling()}returns a data.frame of sampled observation records.
\item \code{clear.sampling()} removes all sampling from event table.
\item \code{get.obs.rec()} returns a logical vector indicating whether each
event record represents an observation or not.
\item \code{get.nobs()} returns the number of observation (not dosing) records.
\item \code{get.units()} returns a two-element character vector with the
dosing and time units, respectively
\item \code{copy()} makes a copy of the current event table. To create
a copy of an event table object use \code{qd2 <- qd$copy()}
\item \code{expand()} Expands the event table for multi-subject solving.
This is done by \code{qd$expand(400)} for a 400 subject data expansion
}
}
\description{
Initializes an object of class \sQuote{EventTable} with methods for
adding and querying dosing and observation records
}
\examples{
# create dosing and observation (sampling) events
# QD 50mg dosing, 5 days followed by 25mg 5 days
#
qd <- eventTable(amount.units = "mg", time.units = "days")
#
qd$add.dosing(dose = 50, nbr.doses = 5, dosing.interval = 1, do.sampling = FALSE)
#
# sample the system's drug amounts hourly the first day, then every 12 hours
# for the next 4 days
qd$add.sampling(seq(from = 0, to = 1, by = 1 / 24))
qd$add.sampling(seq(from = 1, to = 5, by = 12 / 24))
#
# print(qd$get.dosing())     # table of dosing records
print(qd$get.nobs()) # number of observation (not dosing) records
#
# BID dosing, 5 days
bid <- eventTable("mg", "days") # only dosing
bid$add.dosing(
  dose = 10000, nbr.doses = 2 * 5,
  dosing.interval = 12, do.sampling = FALSE
)
#
# Use the copy() method to create a copy (clone) of an existing
# event table (simple assignments just create a new reference to
# the same event table object (closure)).
#
bid.ext <- bid$copy() # three-day extension for a 2nd cohort
bid.ext$add.dosing(
  dose = 5000, nbr.doses = 2 * 3,
  start.time = 120, dosing.interval = 12, do.sampling = FALSE
)

# You can also use the Piping operator to create a table

qd2 <- eventTable(amount.units = "mg", time.units = "days") \%>\%
  add.dosing(dose = 50, nbr.doses = 5, dosing.interval = 1, do.sampling = FALSE) \%>\%
  add.sampling(seq(from = 0, to = 1, by = 1 / 24)) \%>\%
  add.sampling(seq(from = 1, to = 5, by = 12 / 24))
# print(qd2$get.dosing())     # table of dosing records
print(qd2$get.nobs()) # number of observation (not dosing) records

# Note that piping with \%>\% will update the original table.

qd3 <- qd2 \%>\% add.sampling(seq(from = 5, to = 10, by = 6 / 24))
print(qd2$get.nobs())
print(qd3$get.nobs())
}
\seealso{
\code{\link[=et]{et()}}, \code{\link[=RxODE]{RxODE()}}
}
\author{
Matthew Fidler, Melissa Hallow and Wenping Wang
}
\concept{Nonlinear regression}
\concept{Pharmacodynamics (PD)}
\concept{Pharmacokinetics (PK)}
\concept{ordinary differential equations}
\keyword{data}
\keyword{models}
