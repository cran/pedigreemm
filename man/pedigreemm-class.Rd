\name{pedigreemm-class}
\docType{class}
\alias{pedigreemm-class}
\title{Pedigree-based mixed-effects model fits}
\description{A mixed-effects model fit by \code{\link{pedigreemm}}.
  This class extends class \code{"\link[lme4:mer-class]{mer}"} class and includes one
  additional slot, \code{relfac}, which is a list of (left) Cholesky
  factors of the relationship matrices derived from 
  \code{"\linkS4class{pedigree}"} objects.}
\section{Objects from the Class}{Objects are created by calls to the
  \code{\link{pedigreemm}} function.}
\section{Slots}{
  \describe{
    \item{\code{relfac}:}{A list of relationship matrix factors.  All
      other slots are inherited from class \code{"\link[lme4:mer-class]{mer}"}.}
    \item{\code{env}, \code{nlmodel}, \code{frame},
      \code{call}, \code{flist}, \code{X}, \code{Zt},
      \code{pWt}, \code{offset}, \code{y}, \code{Gp},
      \code{dims}, \code{ST}, \code{V}, \code{A}, \code{Cm},
      \code{Cx}, \code{L}, \code{deviance}, \code{fixef}, \code{ranef},
      \code{u}, \code{eta}, \code{mu}, \code{muEta}, \code{var}, \code{resid},
      \code{sqrtXWt}, \code{sqrtrWt}, \code{RZX}, \code{RX}, \code{Xst},
      \code{ghw}, \code{ghx}:}{all other slots are inherited from the
      superclass \code{"\link[lme4:mer-class]{mer}"} in the
      \code{lme4} package.}
  }
}
\section{Extends}{
Class \code{"\link[lme4:mer-class]{mer}"}, directly.
}
\section{Methods}{
No methods defined with class "pedigreemm" in the signature.
}
\seealso{\code{\link{pedigreemm}}
}
\examples{
showClass("pedigreemm")
}
\keyword{classes}
