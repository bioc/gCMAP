\name{memorize}
\alias{memorize}

\title{
  Create a new NChannelSet instance by selecting specific channels and
  load BigMatrix assayData into memory.
}
\description{
  This function converts BigMatrix objects stored in the assayData slot
  of NChannelSets into standard matrices, loading them fully into
  memory. Standard in-memory objects are returned unchanged. Requires \code{bigmemoryExtras} package to process BigMatrix objects.
}
\usage{
memorize(object, names, ...)
}

\arguments{
  \item{object}{
    An NChannelSet object.
}
  \item{names}{
    Character vector of named channels (default: all channels are returned).
}
  \item{\dots}{
    Additional arguments.
}
}

\value{
Instance of the same class as 'object'.  
}

\author{
Thomas Sandmann
}
\note{
This function can be applied to any class inheriting from the virtual
eSet class. For non NChannelSets, meta data may not be transferred correctly.
}

\seealso{
  \code{\link{eSet}}
  \code{\link{memorize}}
  \code{\link{selectChannels}}
}
\examples{
if( suppressWarnings(require("bigmemoryExtras", 
     quietly=TRUE, character.only=TRUE))){
## load ExpressionSet
data("sample.ExpressionSet") ## from Biobase
sample.ExpressionSet ## two assayDataElements: exprs, se

## 'exprs' data matrix
class( assayDataElement( sample.ExpressionSet, "exprs" ) )

## convert assayData to BigMatrix objects
storage.file <- tempfile() ## create path & basemane for BigMatrices
storage.file

on.disk <- eSetOnDisk( sample.ExpressionSet, storage.file )
on.disk ## ExpressionSet
dir(dirname( storage.file )) ## created 3 files per channel

class( assayDataElement( on.disk, "exprs" ) ) ## BigMatrix object

## BigMatrix objects are loaded only upon subsetting
assayDataElement( on.disk, "exprs") ## retrieves BigMatrix, NOT matrix
head( assayDataElement( on.disk, "exprs")[,] ) ## retrieves matrix

## convert back to standard in-memory ExpressionSet
in.memory <- memorize( on.disk ) ## all channels
class( assayDataElement( in.memory, "exprs" ) ) ## matrix object
assayDataElementNames( in.memory )

in.memory <- memorize( on.disk, names="exprs" ) ## channel "exprs" only
assayDataElementNames( in.memory )

## remove tempfiles generated for this example
unlink(paste(storage.file,"*", sep=""))  
}
}
