\name{pedCows}
\alias{pedCows}
\docType{data}
\title{Pedigree of the cows in milk}
\description{
  A \code{\linkS4class{pedigree}} object giving (part of) the pedigree
  of the cows in the \code{\link{milk}} data frame.
}
%\usage{data(pedCows)}
\format{
  The format is:
Formal class 'pedigree' [package "pedigreemm"] with 3 slots
  ..@ sire : int [1:6547] NA NA NA NA NA NA NA NA NA NA ...
  ..@ dam  : int [1:6547] NA NA NA NA NA NA NA NA NA NA ...
  ..@ label: chr [1:6547] "1" "2" "3" "4" ...
}
\examples{
str(pedCows)
}
\seealso{\code{\link{milk}}}
\keyword{datasets}
