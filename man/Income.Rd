\name{Income}
\alias{Income}
\docType{data}
\title{Mean Household Income at the Local Authorities of Greece in 2011}
\description{
Municipality centroids and socio-economic variables aggregated to the
new local authority geography in Greece (Programme \emph{Kallikratis}).
}
\usage{data("Income")}
\format{
A data frame with 325 observations on the following 6 variables:
  \describe{
    \item{\code{X}}{numeric; X coordinate of the municipality centroid.}
    \item{\code{Y}}{numeric; Y coordinate of the municipality centroid.}
    \item{\code{UnemrT01}}{numeric; total unemployment rate in 2001
      (Census).}
    \item{\code{PrSect01}}{numeric; proportion of the economically active
      population working in the primary sector (mainly agriculture,
      fishery and forestry) in 2001 (Census).}
    \item{\code{Foreig01}}{numeric; proportion of inhabitants without
      Greek citizenship in 2001 (Census).}
    \item{\code{Income01}}{numeric; mean household income, in Euros,
      earned in 2001 and declared in 2002 tax forms.}
  }
}
\details{
The X, Y coordinates refer to the geometric centroids of the 325
municipalities in Greece established by the \emph{Kallikratis} reform in
2011.
}
\source{
The original shapefile of the corresponding polygons was published by
the Hellenic Statistical Authority (EL.STAT.) at
\url{https://www.statistics.gr/}. The population, employment, citizenship
and employment-sector data were also published by EL.STAT.; the income
data were obtained from the General Secretariat of Information Systems
in Greece at the postcode level. All data were aggregated to the new
municipalities by the author.
}
\references{
Kalogirou, S. and Hatzichristos, T. (2007) A spatial modelling framework
for income estimation. \emph{Spatial Economic Analysis}, 2(3), 297-316.
\doi{10.1080/17421770701576921}.

Kalogirou, S. (2010) Spatial inequalities in income and post-graduate
educational attainment in Greece. \emph{Journal of Maps}, 6(1), 393-400.
\doi{10.4113/jom.2010.1095}.
}
\examples{
data(Income)
boxplot(Income$Income01)
hist(Income$PrSect01)
}
\keyword{datasets}
