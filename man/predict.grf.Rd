\name{predict.grf}
\alias{predict.grf}
\title{Predict Method for a Geographically Weighted Random Forest}
\description{
Generates predictions for new observations using a fitted
\code{\link{grf}} object. For each new observation the local random
forest fitted at the geographically nearest training location is used.
The user can blend local and global predictions through the
\code{local.w} and \code{global.w} weights.
}
\usage{
\method{predict}{grf}(object, new.data, x.var.name, y.var.name,
                      local.w = 1, global.w = 0, ...)
}
\arguments{
  \item{object}{an object of class \code{"grf"} created by
    \code{\link{grf}} with \code{forests = TRUE}.}
  \item{new.data}{a data frame containing the predictors required by the
    model formula and the X and Y coordinates of the new observations.}
  \item{x.var.name}{the name of the column in \code{new.data} that
    contains the X coordinate.}
  \item{y.var.name}{the name of the column in \code{new.data} that
    contains the Y coordinate.}
  \item{local.w}{numeric weight applied to the local model prediction.
    Default is \code{1}.}
  \item{global.w}{numeric weight applied to the global model prediction.
    Default is \code{0}.}
  \item{\dots}{further arguments passed to the \code{ranger}
    \code{\link[stats]{predict}} method (e.g. \code{num.threads}).}
}
\details{
For every row \code{i} in \code{new.data} the function computes the
Euclidean distance between
\code{(new.data[i, x.var.name], new.data[i, y.var.name])} and every
training location in \code{object$Locations}, picks the local random
forest fitted at the nearest training location, and combines its
prediction with the prediction of the global random forest using the
weights \code{local.w} and \code{global.w}.
}
\value{
A numeric vector of predictions of length \code{nrow(new.data)}.
}
\references{
Georganos, S., Grippa, T., Niang Gadiaga, A., Linard, C., Lennert, M.,
Vanhuysse, S., Mboga, N., Wolff, E., Kalogirou, S. (2019)
Geographical Random Forests: A Spatial Extension of the Random Forest
Algorithm to Address Spatial Heterogeneity in Remote Sensing and
Population Modelling. \emph{Geocarto International}, \doi{10.1080/10106049.2019.1595177}.
}
\author{Stamatis Kalogirou \email{stamatis.science@gmail.com},
        Stefanos Georganos \email{stefanos.georganos@kau.se}}
\note{For tutorials, publications and contact information, see the
package web site.}
\seealso{\code{\link{grf}}}
\examples{
  \donttest{
    set.seed(1)
    RDF <- random.test.data(8, 8, 3)
    Coords <- RDF[, 4:5]
    m <- grf(dep ~ X1 + X2, dframe = RDF, bw = 12,
             kernel = "adaptive", coords = Coords,
             ntree = 100, mtry = 1, nthreads = 1,
             print.results = FALSE, progress = FALSE)

    set.seed(2)
    RDF.Test <- random.test.data(2, 2, 3)
    predict(m, RDF.Test, x.var.name = "X", y.var.name = "Y",
            local.w = 1, global.w = 0)
  }
  \dontrun{
    ## Real-world demonstration on the Greek Income dataset.
    ## Not run by R CMD check because fitting a full GRF on 325
    ## municipalities exceeds the few-second example budget.
    data(Income)
    Coords <- Income[, 1:2]

    m <- grf(Income01 ~ UnemrT01 + PrSect01, dframe = Income, bw = 60,
             kernel = "adaptive", coords = Coords)

    set.seed(123)
    x <- runif(20, min = 142498, max = 1001578)
    y <- runif(20, min = 3855768, max = 4606754)
    u <- runif(20, min = 5,      max = 50)
    p <- runif(20, min = 0,      max = 100)
    f <- runif(20, min = 2,      max = 30)
    df2 <- data.frame(X = x, Y = y, UnemrT01 = u, PrSect01 = p, Foreig01 = f)

    predict(m, df2, x.var.name = "X", y.var.name = "Y",
            local.w = 1, global.w = 0)
  }
}
\keyword{spatial}
\keyword{models}
