\name{Income}
\alias{Income}
\docType{data}
\title{Mean household income at lcoal authorities in Greece in 2011}
\description{Municipality centroids and socioeconomic variables aggregated to the new local authority geography in Greece (Programme Kallikratis).}

\usage{data(Income)}
\format{
  A data frame with 325 observations on the following 5 variables.
  \describe{
    \item{\code{X}}{a numeric vector of x coordinates}
    \item{\code{Y}}{a numeric vector of y coordinates}
    \item{\code{UnemrT01}}{a numeric vector of total unemployment rate in 2001 (Census)}
    \item{\code{PrSect01}}{a numeric vector of the proportion of economically active working in the primary financial sector (mainly agriculture; fishery; and forestry in 2001 (Census))}
    \item{\code{Foreig01}}{a numeric vector of proportion of people who do not have the Greek citizenship in 2001 (Census)}
    \item{\code{Income01}}{a numeric vector of mean recorded household income (in Euros) earned in 2001 and declared in 2002 tax forms}

  }
}

\details{
The X,Y coordinates refer to the geometric centroids of the new 325 Municipalities in Greece (Programme Kallikratis) in 2011.}
\source{
The original shapefile of the corresponding polygons is available from the Hellenic Statistical Authority (EL.STAT.) at \url{http://www.statistics.gr/el/digital-cartographical-data}. The population, employment, citizenship and employment sector data is available from the Hellenic Statistical Authority (EL.STAT.) at \url{http://www.statistics.gr/en/home} but were aggregated to the new municipalities by the author. The income data are available from the General Secretariat of Information Systems in Greece at \url{http://www.gsis.gr/} at the postcode level of geography and were aggregated to the new municipalities by the author.
}
\references{
Kalogirou, S., and Hatzichristos, T. (2007). A spatial modelling framework for income estimation. Spatial Economic Analysis, 2(3), 297-316. \url{http://www.tandfonline.com/doi/abs/10.1080/17421770701576921}

Kalogirou, S. (2010). Spatial inequalities in income and post-graduate educational attainment in Greece. Journal of Maps, 6(1), 393-400.\url{http://www.tandfonline.com/doi/abs/10.4113/jom.2010.1095}

Kalogirou, S. (2013) Testing geographically weighted multicollinearity diagnostics, GISRUK 2013, Department of Geography and Planning, School of Environmental Sciences, University of Liverpool, Liverpool, UK, 3-5 April 2013. \url{http://gisc.gr/?mdocs-file=1140&mdocs-url=false}

}
\examples{
data(Income)
boxplot(Income$Income01)
hist(Income$PrSect01)
}
\keyword{datasets}
\keyword{Greek Municipalities}
\keyword{Income}
