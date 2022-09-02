\name{predict.grf}
\alias{predict.grf}
\title{Predict Method for Geographical Random Forest}
\description{Prediction of test data using the geographical random forest.}

\usage{\method{predict}{grf}(object, new.data, x.var.name, y.var.name, local.w=1, global.w=0,...)
}
\arguments{
  \item{object}{an object that created by the function grf that includes all local forests.}
  \item{new.data}{a data frame containing new data.}
  \item{x.var.name}{the name of the variable with X coordinates.}
  \item{y.var.name}{the name of the variable with Y coordinates.}
  \item{local.w}{weight of the local model predictor allowing semi-local predictions. Default value is 1.}
  \item{global.w}{weight of the global model predictor allowing semi-local predictions. Default value is 0.}
  \item{...}{for other arguments passed to the generic predict functions. For example you may pass here the number of threats}}
\details{A Geographical Random Forest prediction on unknown data. The nearest local random forest model in coordinate space is used to predict in each unknown y-variable location.
}

\value{vector of predicted values}

\references{Stefanos Georganos, Tais Grippa, Assane Niang Gadiaga, Catherine Linard, Moritz Lennert, Sabine Vanhuysse, Nicholus Odhiambo Mboga, Eléonore Wolff & Stamatis Kalogirou (2019) Geographical Random Forests: A Spatial Extension of the Random Forest Algorithm to Address Spatial Heterogeneity in Remote Sensing and Population Modelling, Geocarto International, DOI: 10.1080/10106049.2019.1595177}

\author{Stamatis Kalogirou <stamatis@lctools.science>, Stefanos Georganos <sgeorgan@ulb.ac.be>}

\note{This function is under development. There should be improvements in future versions of the package \code{SpatialML}. Any suggestion is welcome!}

\seealso{
\code{\link{grf}}
}

\examples{
  \dontrun{
      RDF <- random.test.data(10,10,3)
      Coords<-RDF[ ,4:5]
      grf <- grf(dep ~ X1 + X2, dframe=RDF, bw=10,
                kernel="adaptive", coords=Coords)

      RDF.Test <- random.test.data(2,2,3)

      predict.grf(grf, RDF.Test, x.var.name="X", y.var.name="Y", local.w=1, global.w=0)
  }
  \donttest{
      #Load the sample data
      data(Income)

      #Create the vector of XY coordinates
      Coords<-Income[,1:2]

      #Fit local model
      grf <- grf(Income01 ~ UnemrT01 + PrSect01, dframe=Income, bw=60,
                kernel="adaptive", coords=Coords)

      #Create New Random Data - XY coordinates inside the sample data map extend
      x<-runif(20, min = 142498, max = 1001578)
      y<-runif(20, min = 3855768, max = 4606754)
      u<-runif(20, min = 5, max = 50)
      p<-runif(20, min = 0, max = 100)
      f<-runif(20, min = 2, max = 30)
      df2<-data.frame(X=x, Y= y,  UnemrT01=u,  PrSect01=p,  Foreig01=f)

      #Make predictions using the local model
      predict.grf(grf, df2, x.var.name="X", y.var.name="Y", local.w=1, global.w=0)
  }
}

\keyword{local random forest}
