\name{rf.mtry.optim}
\alias{rf.mtry.optim}
\title{Optimal mtry}

\description{This function calculates the optimal mtry for a give ramdom forest (RF) model in the range of values 1-7. The optimal mtry value can then be used in the grf model.}

\usage{rf.mtry.optim(formula, dataset)}

\arguments{
  \item{formula}{the local model to be fitted using the same syntax used in the \code{ranger} function of the R package \code{\link{ranger}}. This is a string that is passed to the sub-models' \code{ranger} function. For more details look at the class \code{\link{formula}}.
}
  \item{dataset}{a numeric data frame of at least two suitable variables (one dependent and one independent)}
}

\details{Based on the \code{train} function of the \code{caret} package, this function sets up a grid of tuning parameters for a number of random forest routines, fits each model and calculates a resampling based performance measure to choose the best mtry value.}

\value{A list is returned of class train as in the function \code{train} in the \code{caret} package.}

\references{Kuhn, M. (2008). Building Predictive Models in R Using the caret Package. Journal of Statistical Software, 28(5), 1 - 26. doi:http://dx.doi.org/10.18637/jss.v028.i05}

\author{Stamatis Kalogirou <stamatis@lctools.science>, Stefanos Georganos <sgeorgan@ulb.ac.be>}

\note{This function is under development.}

\examples{

  \donttest{
      data(Income)
      Coords <- Income[ ,1:2]
      results <- rf.mtry.optim(Income01 ~ UnemrT01 + PrSect01, Income)
      }
}
