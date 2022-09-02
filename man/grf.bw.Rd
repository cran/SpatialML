\name{grf.bw}
\alias{grf.bw}
\title{Geographically Weighted Random Forest optimal bandwidth selection}
\description{This function finds the optimal bandwidth for the Geographically Weighted Random Forest algorithm using an exhaustive approach.}
\usage{grf.bw(formula, dataset, kernel="adaptive", coords, bw.min = NULL,
              bw.max = NULL, step = 1, trees=500, mtry=NULL, importance="impurity",
              nthreads = 1, forests = FALSE, weighted = TRUE, ...)}
\arguments{

  \item{formula}{the local model to be fitted using the same syntax used in the \code{ranger} function of the R package \code{\link{ranger}}. This is a string that is passed to the sub-models' \code{ranger} function. For more details look at the class \code{\link{formula}}.}

  \item{dataset}{a numeric data frame of at least two suitable variables (one dependent and one independent)}

  \item{kernel}{the kernel to be used in the regression. Options are "adaptive" (default) or "fixed".}

  \item{coords}{a numeric matrix or data frame of two columns giving the X,Y coordinates of the observations}

  \item{bw.min}{an integer referring to the minimum bandwidth that evaluation starts.}

  \item{bw.max}{an integer referring to the maximum bandwidth that evaluation ends.}

  \item{step}{an integer referring to the step for each iteration of the evaluation between the min and the max bandwidth. Default value is 1.}

  \item{trees}{an integer referring to the number of trees to grow for each of the local random forests.}

  \item{mtry}{the number of variables randomly sampled as candidates at each split. Note that the default values is p/3, where p is number of variables in the formula}

  \item{importance}{feature importance of the dependent variables used as input at the random forest. Default value is "impurity" which refers to the Gini index for classification and the variance of the responses for regression.}

  \item{nthreads}{Number of threads. Default is number of CPUs available. The argument passes to both ranger and predict functions.}

  \item{forests}{a option to save and export (TRUE) or not (FALSE) all the local forests. Default value is FALSE.}

  \item{weighted}{if TRUE the algorithm calculates Geographically Weighted Random Forest using the case.weights option of the package ranger. If FALSE it will calculate local random forests without weighting each observation in the local data set.}

  \item{...}{further arguments passed to the grf and ranger functions}

}
\details{
Geographically Weighted Random Forest (GRF) is a spatial analysis method using a local version of the famous Machine Learning algorithm. It allows for the investigation of the existence of spatial non-stationarity, in the relationship between a dependent and a set of independent variables. The latter is possible by fitting a sub-model for each observation in space, taking into account the neighbouring observations. This technique adopts the idea of the Geographically Weighted Regression, Kalogirou (2003). The main difference between a tradition (linear) GWR and GRF is that we can model non-stationarity coupled with a flexible non-linear model which is very hard to over-fit due to its bootstrapping nature, thus relaxing the assumptions of traditional Gaussian statistics. Essentially, it was designed to be a bridge between machine learning and geographical models, combining inferential and explanatory power. Additionally, it is suited for datasets with numerous predictors, due to the robust nature of the random forest algorithm in high dimensionality.

This function is a first attempt to find the optimal bandwidth for the grf. It uses an exhaustive approach, i.e. it tests sequential nearest neighbour bandwidths within a range and with a user defined step, and returns a list of goodness of fit statistics. It chooses the best bandwidth based on the maximum R2 value of the local model. Future versions of this function will include heuristic methods to find the optimal bandwidth using algorithms such as optim.}


\value{
   \item{tested.bandwidths}{A table with the tested bandwidths and the corresponding R2 of three model configurations: Local that refers to predictions based on the local (grf) model only; Mixed that refers to predictions that equally combine local (grf) and global (rf) model predictors; and Low.Local that refers to a prediction based on the combination of the local model predictors with a weight of 0.25 and the global model predictors with a weight of 0.75).}

  \item{best.bw}{Best bandwidth based on the local model predictions.}
}

\references{Stefanos Georganos, Tais Grippa, Assane Niang Gadiaga, Catherine Linard, Moritz Lennert, Sabine Vanhuysse, Nicholus Odhiambo Mboga, Eléonore Wolff and Stamatis Kalogirou (2019) Geographical Random Forests: A Spatial Extension of the Random Forest Algorithm to Address Spatial Heterogeneity in Remote Sensing and Population Modelling, Geocarto International, DOI: 10.1080/10106049.2019.1595177

Georganos, S. and Kalogirou, S. (2022) A Forest of Forests: A Spatially Weighted and Computationally Efficient Formulation of Geographical Random Forests. ISPRS, International Journal of Geo-Information, 2022, 11, 471. <https://www.mdpi.com/2220-9964/11/9/471>}

\author{Stamatis Kalogirou <stamatis@lctools.science>, Stefanos Georganos <sgeorgan@ulb.ac.be>}

\note{
This function is under development. There should be improvements in future versions of the package \code{SpatialML}. Any suggestion is welcome!
}

\section{Warning}{Large datasets may take long time to evaluate the optimal bandwidth.}

\seealso{
\code{\link{grf}}
}

\examples{
  \dontrun{
      RDF <- random.test.data(8,8,3)
      Coords<-RDF[ ,4:5]
      bw.test <- grf.bw(dep ~ X1 + X2, RDF, kernel="adaptive",
      coords=Coords, bw.min = 20, bw.max = 23, step = 1,
      forests = FALSE, weighted = TRUE)
  }
  \donttest{
      data(Income)
      Coords<-Income[ ,1:2]

      bwe <-grf.bw(Income01 ~ UnemrT01 + PrSect01, Income, kernel="adaptive",
                   coords=Coords, bw.min = 30, bw.max = 80, step = 1,
                   forests = FALSE, weighted = TRUE)

      grf <- grf(Income01 ~ UnemrT01 + PrSect01, dframe=Income, bw=bwe$Best.BW,
                 kernel="adaptive", coords=Coords)
  }
}

\keyword{spatial random forest}
\keyword{predictive analytics}
