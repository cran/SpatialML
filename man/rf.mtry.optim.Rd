\name{rf.mtry.optim}
\alias{rf.mtry.optim}
\title{Optimal mtry for a Random Forest via OOB or Cross-Validation}
\description{
Searches for the value of \code{mtry} that minimises the predictive error
of a Random Forest on a user-supplied grid. Three evaluation strategies
are available:
\itemize{
  \item \code{"oob"} (default): one \code{\link[ranger]{ranger}} fit per
    \code{mtry}; OOB error is the criterion (fast).
  \item \code{"cv"}: a single \code{cv.folds}-fold cross-validation is
    performed for each \code{mtry}.
  \item \code{"repeatedcv"}: \code{cv.repeats} repeats of \code{cv.folds}-fold
    CV are performed for each \code{mtry}.
}
The selected value can then be passed to the \code{mtry} argument of
\code{\link{grf}}.
}
\usage{
rf.mtry.optim(formula, dataset, min.mtry = NULL, max.mtry = NULL,
              mtry.step = 1, num.trees = 500,
              cv.method = c("oob", "repeatedcv", "cv"),
              cv.folds = 10, cv.repeats = 5, num.threads = NULL,
              plot.it = TRUE, verbose = TRUE, ...)
}
\arguments{
  \item{formula}{a model formula (or a character string coercible to
    one).}
  \item{dataset}{a data frame containing the variables of \code{formula}.}
  \item{min.mtry}{integer; lower bound of the \code{mtry} grid. Default
    \code{1}.}
  \item{max.mtry}{integer; upper bound of the \code{mtry} grid. Default
    is the number of predictors in \code{formula}.}
  \item{mtry.step}{integer; step of the \code{mtry} grid. Default
    \code{1}.}
  \item{num.trees}{integer; number of trees in each random forest fit.
    Default \code{500}.}
  \item{cv.method}{character; evaluation strategy. One of \code{"oob"}
    (default), \code{"repeatedcv"} or \code{"cv"}.}
  \item{cv.folds}{integer; number of folds when
    \code{cv.method != "oob"}. Default \code{10}.}
  \item{cv.repeats}{integer; number of repeats when
    \code{cv.method == "repeatedcv"}. Default \code{5}.}
  \item{num.threads}{number of threads passed to
    \code{\link[ranger]{ranger}}. The default (\code{NULL}) lets
    \pkg{ranger} pick a sensible value.}
  \item{plot.it}{logical; if \code{TRUE} (default) the \code{mtry} vs.
    RMSE curve is plotted.}
  \item{verbose}{logical; if \code{TRUE} (default) progress messages are
    printed for every evaluated \code{mtry}.}
  \item{\dots}{additional arguments forwarded to
    \code{\link[ranger]{ranger}} (e.g. \code{importance = "impurity"}).}
}
\details{
The criterion to minimise is the average RMSE across folds (or the OOB
RMSE in the \code{"oob"} case). The average squared correlation between
predictions and observations is also recorded as the \code{Rsquared}
column.

For reproducible results call \code{set.seed()} before
\code{rf.mtry.optim()}.

\strong{Compatibility note:} the function no longer returns a
\code{caret} \code{train} object. The dependency on \pkg{caret} (and
indirectly on \pkg{randomForest}) has been removed and the
cross-validation is now implemented directly on top of
\code{\link[ranger]{ranger}}. Code that read \code{$bestTune$mtry} from
the previous return value should now read \code{$best.mtry}.
}
\value{
A list with the following components:
  \item{best.mtry}{The \code{mtry} value with the lowest average RMSE.}
  \item{results}{A data frame with one row per evaluated \code{mtry} and
    the columns \code{mtry}, \code{RMSE}, \code{Rsquared}, \code{SDRMSE}
    and \code{SDRsq}. Standard deviations are \code{NA} for
    \code{cv.method == "oob"}.}
  \item{cv.method}{The evaluation strategy used.}
  \item{num.trees}{Number of trees used in each fit.}
  \item{call}{The matched call.}
}
\references{
Wright, M. N. and Ziegler, A. (2017) ranger: A Fast Implementation of
Random Forests for High Dimensional Data in C++ and R.
\emph{Journal of Statistical Software}, 77(1), 1-17.
\doi{10.18637/jss.v077.i01}.

Georganos, S. and Kalogirou, S. (2022) A Forest of Forests: A Spatially
Weighted and Computationally Efficient Formulation of Geographical
Random Forests. \emph{ISPRS International Journal of Geo-Information},
11(9), 471. \doi{10.3390/ijgi11090471}.
}
\author{Stamatis Kalogirou \email{stamatis.science@gmail.com},
        Stefanos Georganos \email{stefanos.georganos@kau.se}}
\seealso{\code{\link{grf}}, \code{\link[ranger]{ranger}}}
\examples{
  \donttest{
    data(Income)
    set.seed(123)
    res <- rf.mtry.optim(Income01 ~ UnemrT01 + PrSect01, dataset = Income,
                         num.trees = 200, num.threads = 1L,
                         cv.method = "oob", plot.it = FALSE,
                         verbose = FALSE)
    res$best.mtry
    res$results
  }
}
\keyword{models}
\keyword{optimize}
