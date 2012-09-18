\name{zScores}
\alias{zScores}

\title{
Function to calculate z-scores from p-values
}

\description{
Function to calculate z-score from a normal distribution from a
two-tailed p-value and sign vector (e.g. log2 fold change). To avoid
-Inf/Inf z-scores, p-values < 'limit' are set to 'limit'.
}
\usage{
zScores(pval, direction=NULL, limit=.Machine$double.xmin)
}
\arguments{
  \item{pval}{
Vector with p-values
}
  \item{direction}{
Vector that will be used to determine the sign of
the z-scores. Only the sign of the values is considered, so any
suitable vectors (e.g. log2 fold change) can be supplied.
}

\item{limit}{
  Numerical (default: .Machine$double.xmin). pvalues < 'limit' will be
  set to 'limit' to avoid Inf/-Inf z-scores. Set to NULL to disable.
  }
}

\value{
  A vector of z-scores
}

\author{
  Thomas Sandmann
}

\seealso{
  \code{\link{qnorm}}
}
