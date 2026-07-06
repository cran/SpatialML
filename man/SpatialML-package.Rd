\name{SpatialML-package}
\alias{SpatialML-package}
\alias{SpatialML}
\docType{package}
\title{
Spatial Machine Learning: Geographically Weighted Random Forest
}
\description{
\pkg{SpatialML} implements a spatial extension of the Random Forest
algorithm. The main functions are \code{\link{grf}} for fitting a
Geographically Weighted Random Forest (GRF), \code{\link{grf.bw}} for
selecting an optimal bandwidth, \code{\link{predict.grf}} (the S3
method dispatched by \code{predict()}) for generating spatial
predictions, and \code{\link{rf.mtry.optim}} for tuning the global
\code{mtry} parameter through cross-validation. The package also ships a
small synthetic-data generator \code{\link{random.test.data}} and a real
data set \code{\link{Income}}.
}
\details{
GRF fits one Random Forest per observation in space, using a local
neighbourhood defined by an \emph{adaptive} kernel (the \eqn{k}{k}
nearest neighbours) or a \emph{fixed} kernel (a Euclidean radius).
Observations within the neighbourhood are optionally weighted by the
bi-square kernel
\eqn{w_{ij} = (1 - (d_{ij} / h)^2)^2}, where \eqn{d_{ij}} is the
distance between observations \eqn{i} and \eqn{j} and \eqn{h} is the
local bandwidth. The approach combines the flexibility and high
predictive accuracy of Random Forests with a treatment of spatial
non-stationarity inspired by Geographically Weighted Regression.

The package uses \code{\link[ranger]{ranger}} as its random-forest
back-end. Undefined local OOB predictions are handled with a quiet
leave-one-out fallback. See
\code{vignette("SpatialML")} for a complete walk-through of a
\code{mtry} tuning -> bandwidth search -> GRF fit -> prediction
workflow.

For tutorials, related publications and contact information visit the
maintainer's website at \url{https://stamatisgeoai.eu/}.
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
\author{
Stamatis Kalogirou \email{stamatis.science@gmail.com} (maintainer,
\url{https://stamatisgeoai.eu/}),
Stefanos Georganos \email{stefanos.georganos@kau.se}.
}
\seealso{
\code{\link{grf}}, \code{\link{grf.bw}}, \code{\link{predict.grf}},
\code{\link{rf.mtry.optim}}, \code{\link{random.test.data}},
\code{\link{Income}}, \code{\link[ranger]{ranger}}.
}
\keyword{package}
\keyword{spatial}
