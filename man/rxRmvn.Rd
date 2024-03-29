% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{rxRmvn}
\alias{rxRmvn}
\title{Simulate from a (truncated) multivariate normal}
\usage{
rxRmvn(
  n,
  mu = NULL,
  sigma,
  lower = -Inf,
  upper = Inf,
  ncores = 1,
  isChol = FALSE,
  keepNames = TRUE,
  a = 0.4,
  tol = 2.05,
  nlTol = 1e-10,
  nlMaxiter = 100L
)
}
\arguments{
\item{n}{Number of random row vectors to be simulated OR the
matrix to use for simulation (faster).}

\item{mu}{mean vector}

\item{sigma}{Covariance matrix for multivariate normal or a list
of covariance matrices. If a list of covariance matrix, each
matrix will simulate \code{n} matrices and combine them to a full
matrix}

\item{lower}{is a vector of the lower bound for the truncated
multivariate norm}

\item{upper}{is a vector of the upper bound for the truncated
multivariate norm}

\item{ncores}{Number of cores used in the simulation}

\item{isChol}{A boolean indicating if \code{sigma} is a cholesky
decomposition of the covariance matrix.}

\item{keepNames}{Keep the names from either the mean or covariance
matrix.}

\item{a}{threshold for switching between methods; They can be
tuned for maximum speed;  There are three cases that are considered:

case 1: a < l < u

case 2: l < u < -a

case 3: otherwise

where l=lower and u = upper}

\item{tol}{When case 3 is used from the above possibilities, the
tol value controls the acceptance rejection and
inverse-transformation;

When abs(u-l)>tol, uses accept-reject from randn}

\item{nlTol}{Tolerance for newton line-search}

\item{nlMaxiter}{Maximum iterations for newton line-search}
}
\value{
If \code{n==integer} (default) the output is an (n x d) matrix
where the i-th row is the i-th simulated vector.

If \code{is.matrix(n)} then the random vector are store in \code{n},
which is provided by the user, and the function returns
\code{NULL} invisibly.
}
\description{
This is simulated with the fast, thread-safe threefry simulator
and can use multiple cores to generate the random deviates.
}
\examples{

## From mvnfast
## Unlike mvnfast, uses threefry simulation

d <- 5
mu <- 1:d

# Creating covariance matrix
tmp <- matrix(rnorm(d^2), d, d)
mcov <- tcrossprod(tmp, tmp)


set.seed(414)
rxRmvn(4, 1:d, mcov)

set.seed(414)
rxRmvn(4, 1:d, mcov)

set.seed(414)
rxRmvn(4, 1:d, mcov, ncores = 2) # r.v. generated on the second core are different

###### Here we create the matrix that will hold the simulated
#  random variables upfront.
A <- matrix(NA, 4, d)
class(A) <- "numeric" # This is important. We need the elements of A to be of class "numeric".

set.seed(414)
rxRmvn(A, 1:d, mcov, ncores = 2) # This returns NULL ...
A # ... but the result is here

## You can also simulate from a truncated normal:

rxRmvn(10, 1:d, mcov, lower = 1:d - 1, upper = 1:d + 1)


# You can also simulate from different matrices (if they match
# dimensions) by using a list of matrices.

matL <- lapply(1:4, function(...) {
  tmp <- matrix(rnorm(d^2), d, d)
  tcrossprod(tmp, tmp)
})


rxRmvn(4, setNames(1:d, paste0("a", 1:d)), matL)
}
\references{
John K. Salmon, Mark A. Moraes, Ron O. Dror, and David
E. Shaw (2011). Parallel Random Numbers: As Easy as 1, 2, 3.
D. E. Shaw Research, New York, NY 10036, USA.

The thread safe multivariate normal was inspired from the \code{mvnfast} package by Matteo Fasiolo \url{https://CRAN.R-project.org/package=mvnfast}

The concept of the truncated multivariate normal was
taken from Zdravko Botev Botev (2017) \doi{10.1111/rssb.12162}
and Botev and L'Ecuyer (2015) \doi{10.1109/WSC.2015.7408180} and
converted to thread safe simulation;
}
\author{
Matthew Fidler, Zdravko Botev and some from Matteo Fasiolo
}
