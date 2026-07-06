\name{grf}
\alias{grf}
\title{Geographically Weighted Random Forest Model}
\description{
Fits a local version of the Random Forest algorithm. The function fits one
sub-model for every observation in the data set, where each local model is
trained on the geographically nearest neighbours of the focal observation
(adaptive kernel) or on every observation within a fixed distance
(fixed kernel). Local observations may further be weighted with a
bi-square kernel to reproduce a Geographically Weighted Random Forest.
}
\usage{
grf(formula, dframe, bw, kernel = c("adaptive", "fixed"), coords,
    ntree = 500, mtry = NULL, importance = "impurity", nthreads = NULL,
    forests = TRUE, geo.weighted = TRUE, print.results = TRUE,
    progress = print.results, oob.fallback = c("loo", "inbag"), ...)
}
\arguments{
  \item{formula}{a formula (or character string coercible to one) specifying
    the local model to be fitted, using the syntax of
    \code{\link[ranger]{ranger}}.}

  \item{dframe}{a data frame containing the dependent and independent
    variables that appear in \code{formula}.}

  \item{bw}{a positive number representing either the number of nearest
    neighbours (when \code{kernel = "adaptive"}) or the bandwidth in the
    units of \code{coords} (when \code{kernel = "fixed"}).}

  \item{kernel}{the type of kernel used to define the neighbourhood of each
    focal observation. One of \code{"adaptive"} (default) or
    \code{"fixed"}.}

  \item{coords}{a numeric matrix or data frame with two columns giving the
    X and Y coordinates of the observations in \code{dframe}. The number of
    rows must match \code{nrow(dframe)}.}

  \item{ntree}{integer; the number of trees to grow for each (global and
    local) random forest. Default is \code{500}.}

  \item{mtry}{the number of variables randomly sampled as candidates at
    each split. The default is \code{max(floor(p / 3), 1)} where \code{p}
    is the number of predictors in \code{formula}.}

  \item{importance}{variable importance measure. Default is
    \code{"impurity"} (see \code{\link[ranger]{ranger}}).}

  \item{nthreads}{number of threads used by \code{\link[ranger]{ranger}}
    and by \code{\link[stats]{predict}}. The default (\code{NULL}) lets
    \pkg{ranger} pick a sensible value.}

  \item{forests}{logical; if \code{TRUE} (default) all local forests are
    stored in the returned object. Required if you intend to call
    \code{\link{predict.grf}} on new data.}

  \item{geo.weighted}{logical; if \code{TRUE} (default) the local random
    forests are fitted using bi-square geographical case weights. If
    \code{FALSE}, local forests are fitted on the same neighbourhood but
    without weighting (i.e. local but unweighted).}

  \item{print.results}{logical; if \code{TRUE} (default) the function
    prints global and local model summaries to the console.}

  \item{progress}{logical; if \code{TRUE} a text progress bar is shown
    while the local random forests are fitted (only in interactive
    sessions). Defaults to \code{print.results}.}

  \item{oob.fallback}{character; fallback used when \pkg{ranger} cannot
    produce an out-of-bag prediction for the focal observation. The
    default \code{"loo"} fits one additional local forest after removing
    the focal observation and predicts it out-of-sample. \code{"inbag"}
    uses the in-sample local prediction.}

  \item{\dots}{additional arguments passed to \code{\link[ranger]{ranger}}.}
}
\details{
Geographically Weighted Random Forest (GRF) is a spatial analysis method
that fits a local version of the Random Forest algorithm to investigate
spatial non-stationarity in the relationship between a dependent variable
and a set of predictors. A sub-model is fitted for every observation in
space, considering only its neighbouring observations. The technique is
inspired by Geographically Weighted Regression (Kalogirou, 2003) but
replaces the linear local model by a non-linear ensemble that is robust
in high dimensions and less prone to over-fitting thanks to its
bootstrapping nature, which relaxes the distributional assumptions of
classical Gaussian statistics.

For each focal observation \eqn{i} the function:
\enumerate{
  \item identifies the local neighbourhood of \eqn{i} according to
    \code{kernel} and \code{bw};
  \item computes bi-square geographical weights
    \eqn{w_{ij} = (1 - (d_{ij} / h)^2)^2} where \eqn{d_{ij}} is the
    distance between observations \eqn{i} and \eqn{j} and \eqn{h} is the
    local bandwidth;
  \item fits a \pkg{ranger} random forest on the local sub-sample (with
    case weights when \code{geo.weighted = TRUE});
  \item stores the local variable importance, the OOB and in-sample
    predictions for \eqn{i}, and (optionally) the entire local forest.
}

For very small neighbourhoods or highly concentrated case weights,
\pkg{ranger} can occasionally have no out-of-bag trees for the focal
observation. In that case the OOB prediction is undefined. The default
\code{oob.fallback = "loo"} avoids noisy retry warnings by fitting one
extra local forest with the focal row removed and using that
out-of-sample prediction for the focal row.

The procedure is described in Georganos et al. (2019) and Georganos and
Kalogirou (2022).
}
\value{
A list of class \code{"grf"} with the following components:
  \item{Global.Model}{A \code{ranger} object containing the global random
    forest fitted on \code{dframe}.}
  \item{Locations}{The \code{coords} object passed to the function.}
  \item{Local.Variable.Importance}{A data frame with one row per
    observation and one column per predictor giving the local feature
    importance.}
  \item{LGofFit}{A data frame with the observed values and local
    goodness-of-fit statistics: OOB and in-sample predicted values,
    residuals, mean squared error, R-squared and the number of attempts
    needed to obtain the OOB or fallback prediction (\code{LPerm}).}
  \item{Forests}{When \code{forests = TRUE}, a list of length
    \code{nrow(dframe)} containing all local \code{ranger} forests.}
  \item{LocalModelSummary}{A list with the average local variable
    importance, mean squared error (OOB and predicted), R-squared (OOB and
    predicted), and AIC / AICc statistics.}
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

Kalogirou, S. (2003) The Statistical Analysis and Modelling of Internal
Migration Flows within England and Wales. PhD Thesis, School of
Geography, Politics and Sociology, University of Newcastle upon Tyne.
}
\author{Stamatis Kalogirou \email{stamatis.science@gmail.com},
        Stefanos Georganos \email{stefanos.georganos@kau.se}}
\note{For tutorials, publications and contact information, see the
package web site.}
\section{Warning}{Large data sets may take a long time to calibrate. A
high number of observations combined with \code{forests = TRUE} can
produce a very large output object.}
\seealso{\code{\link{predict.grf}}, \code{\link{grf.bw}},
         \code{\link[ranger]{ranger}}}
\examples{
  \donttest{
    set.seed(1)
    RDF <- random.test.data(8, 8, 3)
    Coords <- RDF[, 4:5]
    m <- grf(dep ~ X1 + X2, dframe = RDF, bw = 12,
             kernel = "adaptive", coords = Coords,
             ntree = 100, mtry = 1, nthreads = 1,
             print.results = FALSE, progress = FALSE)
  }
  \dontrun{
    ## Real-world demonstration on the Greek Income dataset.
    ## Not run by R CMD check because a full GRF on 325 municipalities
    ## with 500 trees per local fit takes longer than a CRAN example
    ## should.
    data(Income)
    Coords <- Income[, 1:2]
    m <- grf(Income01 ~ UnemrT01 + PrSect01, dframe = Income, bw = 60,
             kernel = "adaptive", coords = Coords)
  }
}
\keyword{spatial}
\keyword{models}
\keyword{regression}
