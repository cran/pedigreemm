\name{pedigreemm}
\alias{pedigreemm}
\title{Fit mixed-effects models incorporating pedigrees}
\description{
  Fit linear or generalized linear mixed models incorporating the
  effects of a pedigree.
}
\usage{
pedigreemm(formula, data, family = NULL, REML = TRUE,
           pedigree = list(), control = list(),
           start = NULL, verbose = FALSE, subset, weights,
           na.action, offset, contrasts = NULL, model = TRUE,
           x = TRUE, ...)
}
\arguments{
  \item{pedigree}{a named list of \code{\linkS4class{pedigree}} objects.
    The names must correspond to the names of grouping factors for
    random-effects terms in the \code{formula} argument.}
  \item{formula}{as in \code{\link[lme4]{lmer}}}
  \item{data}{as in \code{\link[lme4]{lmer}}}
  \item{family}{as in \code{\link[lme4]{lmer}}}
  \item{REML}{as in \code{\link[lme4]{lmer}}}
  \item{control}{as in \code{\link[lme4]{lmer}}}
  \item{start}{as in \code{\link[lme4]{lmer}}}
  \item{verbose}{as in \code{\link[lme4]{lmer}}}
  \item{subset}{as in \code{\link[lme4]{lmer}}}
  \item{weights}{as in \code{\link[lme4]{lmer}}}
  \item{na.action}{as in \code{\link[lme4]{lmer}}}
  \item{offset}{as in \code{\link[lme4]{lmer}}}
  \item{contrasts}{as in \code{\link[lme4]{lmer}}}
  \item{model}{as in \code{\link[lme4]{lmer}}}
  \item{x}{as in \code{\link[lme4]{lmer}}}
  \item{\dots}{as in \code{\link[lme4]{lmer}}}
}
\details{
  All arguments to this function are the same as those to the function
  \code{\link[lme4]{lmer}} except \code{pedigree} which must be a named list of
  \code{\linkS4class{pedigree}} objects.  Each name (frequently there is
  only one) must correspond to the name of a grouping factor in a
  random-effects term in the \code{formula}.  The observed levels
  of that factor must be contained in the pedigree.  For each pedigree
  the (left) Cholesky factor of the
  relationship matrix restricted to the observed levels is calculated
  using \code{\link{relfactor}} and applied to the model matrix for that
  term.
}
\value{
  a \code{\linkS4class{pedigreemm}} object.
}
\seealso{\code{\linkS4class{pedigreemm}}, \code{\linkS4class{pedigree}},
  \code{\link{relfactor}}.}
%\examples{}
\keyword{models}
