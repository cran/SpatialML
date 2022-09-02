\name{grf}
\alias{grf}
\title{Geographically Weighted Random Forest}
\description{This function refers to a local version of the popular Random Forest algorithm.}
\usage{grf(formula, dframe, bw, kernel, coords, ntree=500, mtry=NULL,
           importance="impurity", nthreads = NULL, forests = TRUE,
           weighted = TRUE, print.results=TRUE, ...)}
\arguments{

  \item{formula}{the local model to be fitted using the same syntax used in the \code{ranger} function of the R package \code{\link{ranger}}. This is a string that is passed to the sub-models' \code{ranger} function. For more details look at the class \code{\link{formula}}.}

  \item{dframe}{a numeric data frame of at least two suitable variables (one dependent and one independent)}

  \item{bw}{a positive number that may be an integer in the case of an "adaptive kernel" or a real in the case of a "fixed kernel". In the first case, the integer denotes the number of nearest neighbours, whereas in the latter case the real number refers to the bandwidth (in meters if the coordinates provided are Cartesian).}

  \item{kernel}{the kernel to be used in the regression. Options are "adaptive" or "fixed".}

  \item{coords}{a numeric matrix or data frame of two columns giving the X,Y coordinates of the observations}

  \item{ntree}{an integer referring to the number of trees to grow for each of the local random forests.}

  \item{mtry}{Number of variables randomly sampled as candidates at each split. Note that the default values is p/3, where p is number of variables in the formula}

  \item{importance}{Feature importance of the dependent variables used as input at the random forest. Default value is "impurity" which refers to the Gini index for classification and the variance of the responses for regression.}

  \item{nthreads}{Number of threads. Default is number of CPUs available. The argument passes to both rnager and predict functions.}

  \item{forests}{a option to save and export (TRUE) or not (FALSE) all the local forests}

  \item{weighted}{if TRUE the algorithm calculates Geographically Weighted Random Forest using the case.weights option of the packare ranger. If FALSE it will calculate local random forests without weighting each observation in the local data set.}

  \item{print.results}{a option to print in the console (TRUE) or not (FALSE) the summary of the analysis}

  \item{...}{further arguments passed to the ranger function}

}
\details{
Geographically Weighted Random Forest (GRF) is a spatial analysis method using a local version of the famous Machine Learning algorithm. It allows for the investigation of the existence of spatial non-stationarity, in the relationship between a dependent and a set of independent variables. The latter is possible by fitting a sub-model for each observation in space, taking into account the neighbouring observations. This technique adopts the idea of the Geographically Weighted Regression, Kalogirou (2003). The main difference between a tradition (linear) GWR and GRF is that we can model non-stationarity coupled with a flexible non-linear model which is very hard to overfit due to its bootstrapping nature, thus relaxing the assumptions of traditional Gaussian statistics. Essentially, it was designed to be a bridge between machine learning and geographical models, combining inferential and explanatory power. Additionally, it is suited for datasets with numerous predictors, due to the robust nature of the random forest algorithm in high dimensionality. }


\value{
   \item{Global.Model}{A ranger object of the global random forest model}

   \item{Locations}{a numeric matrix or data frame of two columns giving the X,Y coordinates of the observations}

   \item{Local.Variable.Importance}{a numeric data frame with the local feature importance for each predictor in each local random forest model}

   \item{LGofFit}{a numeric data frame with residuals and local goodness of fit statistics.}

  \item{Forests}{all local forests.}

  \item{lModelSummary}{Local Model Summary and goodness of fit statistics.}
}


\references{Stefanos Georganos, Tais Grippa, Assane Niang Gadiaga, Catherine Linard, Moritz Lennert, Sabine Vanhuysse, Nicholus Odhiambo Mboga, Eléonore Wolff & Stamatis Kalogirou (2019) Geographical Random Forests: A Spatial Extension of the Random Forest Algorithm to Address Spatial Heterogeneity in Remote Sensing and Population Modelling, Geocarto International, DOI: 10.1080/10106049.2019.1595177

Georganos, S. and Kalogirou, S. (2022) A Forest of Forests: A Spatially Weighted and Computationally Efficient Formulation of Geographical Random Forests. ISPRS, International Journal of Geo-Information, 2022, 11, 471. <https://www.mdpi.com/2220-9964/11/9/471>}

\author{Stamatis Kalogirou <stamatis@lctools.science>, Stefanos Georganos <sgeorgan@ulb.ac.be>}

\note{
This function is under development. There should be improvements in future versions of the package \code{SpatialML}. Any suggestion is welcome!
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
