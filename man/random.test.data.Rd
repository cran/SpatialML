\name{random.test.data}
\alias{random.test.data}
\title{Radmom data generator}
\description{
Generates datasets with random data for modelling including a dependent variable, independent variables and X,Y coordinates.
}
\usage{
random.test.data(nrows = 10, ncols = 10, vars.no = 3, dep.var.dis = "normal",
                xycoords = TRUE)
}
\arguments{
  \item{nrows}{an integer referring to the number of rows for a regular grid}
  \item{ncols}{an integer referring to the number of columns for a regular grid}
  \item{vars.no}{an integer referring to the number of independent variables}
  \item{dep.var.dis}{a character referring to the distribution of the dependent variable. Options are "normal" (default) and "poisson"}
  \item{xycoords}{a logical value indicating whether X,Y coordinates will be created (default) or not.}
}
\details{The creation of a random dataset was necessary here to provide examples to some functions. However, random datasets may be used in simulation studies.}

\value{a dataframe}

\author{Stamatis Kalogirou <stamatis@lctools.science>}

\examples{
RDF <- random.test.data(12,12,3)
}
\keyword{random data}
