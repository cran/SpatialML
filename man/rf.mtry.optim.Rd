\name{rf.mtry.optim}
\alias{rf.mtry.optim}
\title{Optimal mtry}

\description{This function calculates the optimal mtry for a given Random Forest (RF) model in a specified range of values. The optimal mtry value can then be used in the grf model.}

\usage{rf.mtry.optim(formula, dataset, min.mtry=NULL, max.mtry=NULL, mtry.step,
       cv.method="repeatedcv", cv.folds=10, ...)}

\arguments{
  \item{formula}{the model to be fitted using the function \code{train} of the R package \code{caret}.}

  \item{dataset}{a numeric data frame of at least two suitable variables (one dependent and one independent)}

  \item{min.mtry}{the minimum mtry value for its optimisation (function \code{expand.grid})}

  \item{max.mtry}{the maximum mtry value for its optimisation (function \code{expand.grid})}

  \item{mtry.step}{the step in the sequence of mtry values for its optimisation (function \code{expand.grid})}

  \item{cv.method}{the resampling method in the function \code{trainControl} of the R package \code{caret}. Default option is "repeatedcv" and alternative option is "cv".}

  \item{cv.folds}{the number of folds (argument "number" in the function \code{trainControl}). Default value is 10)}

  \item{...}{additional arguments affecting the function \code{trainControl})}

}

\details{Based on the \code{train} function of the \code{caret} package, this function sets up a grid of tuning parameters for a number of random forest routines, fits each model and calculates a resampling based performance measure to choose the best mtry value.}

\value{A list is returned of class train as in the function \code{train} in the \code{caret} package.}

\references{Kuhn, M. (2008). Building Predictive Models in R Using the caret Package. Journal of Statistical Software, 28(5), 1 - 26. doi: <http://dx.doi.org/10.18637/jss.v028.i05>

Georganos, S. and Kalogirou, S. (2022) A Forest of Forests: A Spatially Weighted and Computationally Efficient Formulation of Geographical Random Forests. ISPRS, International Journal of Geo-Information, 2022, 11, 471. <https://www.mdpi.com/2220-9964/11/9/471>}

\author{Stamatis Kalogirou <stamatis@lctools.science>, Stefanos Georganos <sgeorgan@ulb.ac.be>}

\note{This function is under development.}

\examples{

  \donttest{
      data(Income)
      Coords <- Income[ ,1:2]
      results <- rf.mtry.optim(Income01 ~ UnemrT01 + PrSect01, Income)
      }
}
