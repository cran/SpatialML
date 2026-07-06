\name{grf.bw}
\alias{grf.bw}
\title{Optimal Bandwidth Selection for a Geographically Weighted Random Forest}
\description{
Searches an exhaustive grid of bandwidths and returns the one that
maximises the out-of-bag (OOB) R-squared of the local model produced by
\code{\link{grf}}.
}
\usage{
grf.bw(formula, dataset, kernel = "adaptive", coords, bw.min = NULL,
       bw.max = NULL, step = 1, trees = 500, mtry = NULL,
       importance = "impurity", nthreads = 1, forests = FALSE,
       geo.weighted = TRUE, verbose = TRUE, ...)
}
\arguments{
  \item{formula}{a formula (or character string coercible to one) of the
    local model, as accepted by \code{\link[ranger]{ranger}}.}

  \item{dataset}{a data frame containing the dependent and independent
    variables that appear in \code{formula}.}

  \item{kernel}{the kernel used in the regression. Either
    \code{"adaptive"} (default) or \code{"fixed"}.}

  \item{coords}{a numeric matrix or data frame with the X and Y
    coordinates of the observations in \code{dataset}.}

  \item{bw.min}{numeric; lower bound of the bandwidth search.}

  \item{bw.max}{numeric; upper bound of the bandwidth search.}

  \item{step}{numeric; the step size of the search. Default is \code{1}.}

  \item{trees}{integer; the number of trees to grow for each local random
    forest.}

  \item{mtry}{the number of variables randomly sampled as candidates at
    each split. Default is \code{max(floor(p / 3), 1)} where \code{p} is
    the number of predictors in \code{formula}. Values above \code{p}
    are capped at \code{p}.}

  \item{importance}{variable importance measure passed to
    \code{\link[ranger]{ranger}}. Default is \code{"impurity"}.}

  \item{nthreads}{number of threads passed to
    \code{\link[ranger]{ranger}} and \code{\link[stats]{predict}}.}

  \item{forests}{logical; if \code{TRUE} all local forests are kept in
    memory. Default is \code{FALSE} (recommended for bandwidth search).}

  \item{geo.weighted}{logical; if \code{TRUE} fits a Geographically
    Weighted Random Forest using \code{ranger}'s case-weight option.}

  \item{verbose}{logical; if \code{TRUE} (default) prints one progress
    message per evaluated bandwidth and the selected bandwidth.}

  \item{\dots}{additional arguments passed to \code{\link{grf}} and
    \code{\link[ranger]{ranger}}.}
}
\details{
The function evaluates a sequence of candidate bandwidths from
\code{bw.min} to \code{bw.max} in increments of \code{step}. For each
bandwidth a Geographically Weighted Random Forest is fitted via
\code{\link{grf}} and three goodness-of-fit values are recorded:
\describe{
  \item{Local}{R-squared of the local model alone.}
  \item{Mixed}{R-squared of the equally weighted average of local and
    global predictions.}
  \item{Low.Local}{R-squared of \code{0.25 * local + 0.75 * global}.}
}
The optimal bandwidth is selected as the one that maximises the
\emph{Local} R-squared.

If you require fully reproducible bandwidth selection, call
\code{set.seed()} before \code{grf.bw()}.
}
\value{
A list with two components:
  \item{tested.bandwidths}{A data frame with one row per evaluated
    bandwidth and four columns: \code{Bandwidth}, \code{Local},
    \code{Mixed}, \code{Low.Local}.}
  \item{Best.BW}{The bandwidth that yields the highest \code{Local}
    R-squared.}
}
\references{
Georganos, S., Grippa, T., Niang Gadiaga, A., Linard, C., Lennert, M.,
Vanhuysse, S., Mboga, N., Wolff, E., Kalogirou, S. (2019)
Geographical Random Forests: A Spatial Extension of the Random Forest
Algorithm to Address Spatial Heterogeneity in Remote Sensing and
Population Modelling. \emph{Geocarto International}, \doi{10.1080/10106049.2019.1595177}.

Georganos, S. and Kalogirou, S. (2022) A Forest of Forests: A Spatially
Weighted and Computationally Efficient Formulation of Geographical
Random Forests. \emph{ISPRS International Journal of Geo-Information},
11(9), 471. \doi{10.3390/ijgi11090471}.
}
\author{Stamatis Kalogirou \email{stamatis.science@gmail.com},
        Stefanos Georganos \email{stefanos.georganos@kau.se}}
\note{Future versions may include heuristic search strategies
(e.g. \code{\link[stats]{optim}}).}
\section{Warning}{Bandwidth selection on large data sets is
computationally expensive because a full \code{\link{grf}} is fitted at
every bandwidth in the grid.}
\seealso{\code{\link{grf}}, \code{\link[ranger]{ranger}}}
\examples{
  \donttest{
    set.seed(1)
    RDF <- random.test.data(8, 8, 3)
    Coords <- RDF[, 4:5]
    bw.test <- grf.bw(dep ~ X1 + X2, dataset = RDF, kernel = "adaptive",
                      coords = Coords,
                      bw.min = 12, bw.max = 16, step = 2,
                      trees = 100, mtry = 1, nthreads = 1L,
                      forests = FALSE, geo.weighted = TRUE,
                      verbose = FALSE)
    bw.test$Best.BW
  }
  \dontrun{
    ## Real-world demonstration on the Greek Income dataset.
    ## Not run by R CMD check because an exhaustive bandwidth search on
    ## 325 municipalities is intrinsically slow.
    data(Income)
    Coords <- Income[, 1:2]

    bwe <- grf.bw(Income01 ~ UnemrT01 + PrSect01, dataset = Income,
                  kernel = "adaptive", coords = Coords,
                  bw.min = 30, bw.max = 80, step = 5,
                  forests = FALSE, geo.weighted = TRUE)

    m <- grf(Income01 ~ UnemrT01 + PrSect01, dframe = Income,
             bw = bwe$Best.BW, kernel = "adaptive", coords = Coords)
  }
}
\keyword{spatial}
\keyword{optimize}
