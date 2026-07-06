\name{random.test.data}
\alias{random.test.data}
\title{Random Data Generator}
\description{
Generates a synthetic data set composed of one dependent variable, a
number of independent variables drawn from \code{Uniform(0, 1)}, and
optionally the X and Y coordinates of a regular \code{nrows} x
\code{ncols} grid. Useful for examples, unit testing and small
simulation studies.
}
\usage{
random.test.data(nrows = 10, ncols = 10, vars.no = 3,
                 dep.var.dis = c("normal", "poisson"),
                 xycoords = TRUE)
}
\arguments{
  \item{nrows}{integer; number of rows of the regular grid. Default
    \code{10}.}
  \item{ncols}{integer; number of columns of the regular grid. Default
    \code{10}.}
  \item{vars.no}{integer; total number of model variables (1 dependent
    plus \code{vars.no - 1} independent). Must be at least \code{2}.}
  \item{dep.var.dis}{distribution of the dependent variable. One of
    \code{"normal"} (default, drawn from \code{rnorm()}) or
    \code{"poisson"} (drawn from \code{rpois()} with \code{lambda = 7}).}
  \item{xycoords}{logical; if \code{TRUE} (default) two columns \code{X}
    and \code{Y} containing the grid coordinates are appended.}
}
\details{
The function is mostly intended to provide reproducible toy data for
the package examples. The independent variables are named
\code{X1, X2, \dots, X(vars.no - 1)}.

For reproducible results call \code{set.seed()} before
\code{random.test.data()}.
}
\value{
A data frame with \code{nrows * ncols} rows containing one column
\code{dep}, \code{vars.no - 1} columns named \code{X1, X2, ...} and,
when \code{xycoords = TRUE}, two extra columns \code{X} and \code{Y}.
}
\author{Stamatis Kalogirou \email{stamatis.science@gmail.com}}
\examples{
set.seed(1)
RDF <- random.test.data(12, 12, 3)
head(RDF)
}
\keyword{datagen}
