\name{signedRankSumTest}
\alias{signedRankSumTest}

\title{
  An implementation of the Wilcox rank sum test / Mann-Whitney test that
  takes into account the direction / sign of gene set members and
  possibly the correlation between cases}

\description{
  This test evaluates whether the mean rank of statistics of gene set
  members is greater or less than the mean rank of the remaining
  statistic values. It extends the rankSumTestWithCorrelation function
  from the 'limma' package by taking into account the 'sign' of gene
  set members by reversing the ranks of down-regulated genes.
}

\usage{
  signedRankSumTest(statistics, index.up, index.down = NULL,
  input.is.ranks=FALSE, correlation=0, df = Inf, adjust.ties=TRUE) }

\arguments{
  \item{statistics}{
    numeric vector giving values of the test statistic.
}
  \item{index.up}{
    an index vector such that statistics[index.up] contains the values of
    the statistic for the up-regulated genes.
}
  \item{index.down}{
    an index vector such that statistics[index.down] contains the values of
    the statistic for the down-regulated genes.
  }
  \item{correlation}{
    numeric scalar, average correlation between cases in the test group. Cases in the second group are assumed independent of each other
    and the first group.
  }
  \item{df}{
    degrees of freedom which the correlation has been estimated.
  }
  \item{adjust.ties}{
    logical: correct for ties ?
    }
  \item{input.is.ranks}{
    logical: is 'statistics' a vector of ranks ? If FALSE (default),
    ranks are computed. If FALSE, 'statistics' is assumed to represent
    ranks and is used directly.
  }
}
\details{
Please see the rankSumTestWithCorrelation function from the limma
    package for details.
}
\value{
  Numeric vector containing U-statistic, z-score and p-value.
  
}
\references{
  Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test
  accounting for inter-gene correlation. Submitted.
  
  Barry, W.T., Nobel, A.B., and Wright, F.A. (2008). A statistical
  framework for testing functional categories in microarray
  data. Annals of Applied Statistics 2, 286-315.
  
  Zar, JH (1999). Biostatistical Analysis 4th Edition. Prentice-Hall
  International, Upper Saddle River, New Jersey.
}
\author{
Thomas Sandmann
}
\seealso{
\code{\link{rankSumTestWithCorrelation}}
}
\examples{
genes.up   <- c(1:10)
genes.down <- c(21:30)

set.seed(123)
scores <- matrix(rnorm(200), ncol=2)

## the first gene set receives increased /
## decreased scores in the first experiment
scores[genes.up,1]   <- scores[genes.up  ,1] + 1 
scores[genes.down,1] <- scores[genes.down,1] - 1

## significantly greater
signedRankSumTest( statistics = scores[,1],
                   index.up   = genes.up,
                   index.down = genes.down)

## not signficant
signedRankSumTest( statistics = scores[,2],
                   index.up   = genes.up,
                   index.down = genes.down)
}
