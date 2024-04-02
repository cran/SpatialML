\name{grf}
\alias{grf}
\title{Geographically Weighted Random Forest Model}
\description{Fit a local version of the Random Forest algorithm, accounting for spatial non-stationarity.}
\usage{grf(formula, dframe, bw, kernel, coords, ntree=500, mtry=NULL,
           importance="impurity", nthreads = NULL, forests = TRUE,
           geo.weighted = TRUE, print.results=TRUE, ...)}
\arguments{

  \item{formula}{a formula specifying the local model to be fitted, using the syntax of the \code{\link{ranger}} package's \code{\link{ranger}} function.}

  \item{dframe}{a numeric data frame with at least two suitable variables (one dependent and one independent).}

  \item{bw}{a positive number representing either the number of nearest neighbors (for "adaptive kernel") or bandwidth in meters (for "fixed kernel").}

  \item{kernel}{the type of kernel to use in the regression: "adaptive" or "fixed".}

  \item{coords}{a numeric matrix or data frame containing X and Y coordinates of observations.}

  \item{ntree}{an integer referring to the number of trees to grow for each local random forest.}

  \item{mtry}{the number of variables randomly sampled as candidates at each split. Default is \code{p/3}, where \code{p} is the number of variables in the formula.}

  \item{importance}{feature importance measure for the dependent variables used as input in the random forest. Default is "impurity", which refers to the Gini index for classification and the variance of the responses for regression.}

  \item{nthreads}{number of threads for parallel processing. Default is the number of available CPUs. The argument passes to both \code{\link{ranger}} and \code{\link{predict}} functions.}

  \item{forests}{a option to save and export (TRUE) or not (FALSE) all local forests.}

  \item{geo.weighted}{if TRUE, calculate Geographically Weighted Random Forest using case weights. If FALSE, calculate local random forests without weighting each observation.}

  \item{print.results}{a option to print the summary of the analysis (TRUE) or not (FALSE).}

  \item{...}{additional arguments passed to the \code{\link{ranger}} function.}

}
\details{
Geographically Weighted Random Forest (GRF) is a spatial analysis method using a local version of the famous Machine Learning algorithm. It allows for the investigation of the existence of spatial non-stationarity, in the relationship between a dependent and a set of independent variables. The latter is possible by fitting a sub-model for each observation in space, taking into account the neighbouring observations. This technique adopts the idea of the Geographically Weighted Regression, Kalogirou (2003). The main difference between a tradition (linear) GWR and GRF is that we can model non-stationarity coupled with a flexible non-linear model which is very hard to overfit due to its bootstrapping nature, thus relaxing the assumptions of traditional Gaussian statistics. Essentially, it was designed to be a bridge between machine learning and geographical models, combining inferential and explanatory power. Additionally, it is suited for datasets with numerous predictors, due to the robust nature of the random forest algorithm in high dimensionality.


Geographically Weighted Random Forest (GRF) is a spatial analysis method that fits a local version of the Random Forest algorithm for investigating spatial non-stationarity, in the relationship between a dependent and a set of independent variables. The latter is possible by fitting a sub-model for each observation in space, taking into account the neighbouring observations. This technique adopts the idea of the Geographically Weighted Regression, Kalogirou (2003). It models non-stationarity with a flexible non-linear approach, bridging the gap between machine learning and geographical models. The main difference between a tradition (linear) GWR and GRF is that we can model non-stationarity coupled with a flexible non-linear model which is very hard to overfit due to its bootstrapping nature, thus relaxing the assumptions of traditional Gaussian statistics.GRF is suitable for datasets with numerous predictors due to the robustness of the random forest algorithm in high dimensionality.}


\value{
   \item{Global.Model}{A ranger object of the global random forest model.}

   \item{Locations}{a numeric matrix or data frame with X and Y coordinates of observations.}

   \item{Local.Variable.Importance}{anumeric data frame with local feature importance for each predictor in each local random forest model.}

   \item{LGofFit}{a numeric data frame with residuals and local goodness of fit statistics.}

  \item{Forests}{all local forests.}

  \item{lModelSummary}{Local Model Summary and goodness of fit statistics.}
}


\references{Stefanos Georganos, Tais Grippa, Assane Niang Gadiaga, Catherine Linard, Moritz Lennert, Sabine Vanhuysse, Nicholus Odhiambo Mboga, Eléonore Wolff & Stamatis Kalogirou (2019) Geographical Random Forests: A Spatial Extension of the Random Forest Algorithm to Address Spatial Heterogeneity in Remote Sensing and Population Modelling, Geocarto International, DOI: 10.1080/10106049.2019.1595177

Georganos, S. and Kalogirou, S. (2022) A Forest of Forests: A Spatially Weighted and Computationally Efficient Formulation of Geographical Random Forests. ISPRS, International Journal of Geo-Information, 2022, 11, 471. <https://www.mdpi.com/2220-9964/11/9/471>}

\author{Stamatis Kalogirou <stamatis.science@gmail.com>, Stefanos Georganos <sgeorgan@ulb.ac.be>}

\note{
This function is under development, and improvements are expected in future versions of the package \code{SpatialML}. Any suggestions are welcome!
}

\section{Warning}{Large datasets may take long to calibrate. A high number of observations may result in a voluminous forests output.}

\seealso{
\code{\link{predict.grf}}
}

\examples{
  \dontrun{
      RDF <- random.test.data(10,10,3)
      Coords<-RDF[ ,4:5]
      grf <- grf(dep ~ X1 + X2, dframe=RDF, bw=10,
                kernel="adaptive", coords=Coords)
  }
  \donttest{
      data(Income)
      Coords<-Income[ ,1:2]
      grf <- grf(Income01 ~ UnemrT01 + PrSect01, dframe=Income, bw=60,
                kernel="adaptive", coords=Coords)
  }
}

\keyword{spatial random forest}
\keyword{predictive analytics}
