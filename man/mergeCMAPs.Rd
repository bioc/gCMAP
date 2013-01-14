\name{mergeCMAPs}
\alias{mergeCMAPs}

\title{
  This function merged two eSets.
}
\description{
  This function merges two eSet objects, if all of the following
  conditions are met:
  \itemize{
  \item{Both objects 'x' and 'y' have to be instances of the same
    class.}
  \item{'x' and 'y' must be annotated with the same character string in
    their 'annotation' slots.}
  \item{'x' and 'y' must have the same AssayDataElementNames /
    channels.}
  \item{'x' and 'y' must have disctinct sampleNames.}
  \item{'x' and 'y' must have the same varLabels / pData columns.}
}}
\usage{
mergeCMAPs(x, y)
}
\arguments{
  \item{x}{
    An eSet.
}
  \item{y}{
    An eSet of the same class as 'x'
}
}

\value{
  An eSet of the same class as 'x' and 'y'.
}

\author{
Thomas Sandmann, thomas.sandmann@gene.com
}

\examples{
library(Biobase)
data( sample.ExpressionSet)

\dontrun{
## this doesn't work, because 'x' and 'y' have identical sampleNames
mergeCMAPs( sample.ExpressionSet, sample.ExpressionSet)
}

y <- sample.ExpressionSet
sampleNames( y ) <- paste( sampleNames( y ), "y", sep=".")
mergeCMAPs( sample.ExpressionSet, y )
}
