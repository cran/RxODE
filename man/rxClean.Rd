% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{rxClean}
\alias{rxClean}
\title{Cleanup anonymous DLLs by unloading them}
\usage{
rxClean(wd)
}
\arguments{
\item{wd}{What directory should be cleaned; (DEPRECIATED), this no
longer does anything.

This unloads all RxODE anonymous dlls.}
}
\value{
TRUE if successful
}
\description{
This cleans up any RxODE loaded DLLs
}
\author{
Matthew L. Fidler
}
