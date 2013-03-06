\name{center_eSet}
\alias{center_eSet}
\title{
A function to to center columns of eSet channels on either their kernel density peak, their mean or their median.}
\description{
This function works on the eSet assayDataElement specified as 'channel' and sweeps out either the 'peak', ( max of the kernel density), 'mean' or 'median' statistic from each column. A modified eSet containing the centered assayDataElement is returned, with an additional .shift column included in the pData slot recording the shift statistic for each sample.
}
\usage{
center_eSet(eset, channel, center = "peak", report.shifts=FALSE)
}
\arguments{
  \item{eset}{
  An eSet object
}
  \item{channel}{
 A valid channel / AssayDataElementName of 'eset'
}
  \item{center}{
 One of 'peak', 'mean', 'median' or 'none', specifying the statistic to sweep from each column of 'channel' in 'eset'. If 'peak', the max of the kernel density is determined and used a statistic in sweep. If 'none', the original 'eset' is returned.
}
  \item{report.shifts}{Logical, include the shift applied to 'channel' in the pData slot of the returned NChannelSet ?}

}

\value{
An eSet of the same class as 'eset' with the centered 'channel' assayData slot. The swept-out statistic is recorded in the 'channel'.shift column of the phenoData slot. In addition, the median absolute deviation around the center is returned.
}
\author{
Thomas Sandmann
}
\seealso{
\code{\link{sweep}}
}
\examples{
data( gCMAPData )

## column means of uncentered z-scores
round( apply( assayDataElement( gCMAPData, "z"), 2, mean, na.rm=TRUE), 2)

## column means of centered z-scores
centered <- center_eSet( gCMAPData, "z", "mean")
round( apply( assayDataElement( centered, "z"), 2, mean, na.rm=TRUE), 2)
}