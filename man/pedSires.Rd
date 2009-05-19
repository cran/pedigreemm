\name{pedSires}
\alias{pedSires}
\docType{data}
\title{Pedigree of the sires from mastitis}
\description{
  A \code{\linkS4class{pedigree}} object giving (part of) the pedigree
  of the sires from the \code{\link{mastitis}} data frame.  The
  pedigree is traced back on sires only.
}
\format{
  The format is:
Formal class 'pedigree' [package "pedigreemm"] with 3 slots
  ..@ sire : int [1:352] NA NA NA NA NA NA NA NA NA NA ...
  ..@ dam  : int [1:352] NA NA NA NA NA NA NA NA NA NA ...
  ..@ label: chr [1:352] "1" "2" "3" "4" ...
}
\seealso{\code{\link{mastitis}}}
\examples{
str(pedSires)
}
\keyword{datasets}
