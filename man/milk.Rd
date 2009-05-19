\name{milk}
\alias{milk}
\docType{data}
\title{Milk production}
\description{
  Records of the milk production of 3397 lactations from first through fifty parity Holsteins.
  These were 1,359 cows, daughters of 38 sires in 57 herds. 
  The data was downloaded from the USDA internet site. 
  All lactation records represent cows with at least 100 days in milk, with an average of 347 days. 
  Milk yield ranged from 4,065 to 19,345 kg estimated for 305 days, averaging 11,636 kg. 
  There were 1,314, 1,006, 640, 334 and 103 records were from first thorough fifth lactation animals. 
  
}
\usage{data(milk)}
\format{
  A data frame with 3397 observations on the following 9 variables.
  \describe{
    \item{\code{id}}{numeric identifier of cow}
    \item{\code{lact}}{number of lactation for which production is measured}
    \item{\code{herd}}{a factor indicating the herd}
    \item{\code{sire}}{a factor indicating the sire}
    \item{\code{dim}}{number of days in milk for that lactation}
    \item{\code{milk}}{milk production estimated at 305 days}
    \item{\code{fat}}{fat production estimated at 305 days}
    \item{\code{prot}}{protein production estimated at 305 days}
    \item{\code{scs}}{the somatic cell score}
  }
}
 \source{
    USDA web site.
    http://www.aipl.arsusda.gov/
 }
% \references{
% }
\examples{
str(milk)
}
\keyword{datasets}
