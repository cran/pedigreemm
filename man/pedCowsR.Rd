\name{pedCowsR}
\alias{pedCowsR}
\docType{data}
\title{Pedigree of the cows in milk with 0.70 of the information in pedCows}
\description{
  A \code{\linkS4class{pedigree}} object giving (part of) the pedigree
  of the cows in the \code{\link{milk}} data frame.
  This pedigree allows the example with 'milk' to run faster. 
}
%\usage{data(pedCowsR)}
\format{
  The format is:
Formal class 'pedigree' [package "pedigreemm"] with 3 slots
  ..@ sire : int [1:6547] NA NA NA NA NA NA NA NA NA NA ...
  ..@ dam  : int [1:6547] NA NA NA NA NA NA NA NA NA NA ...
  ..@ label: chr [1:6547] "1" "2" "3" "4" ...
}
\examples{
str(pedCowsR)
}
\seealso{\code{\link{milk}}}
\keyword{datasets}
