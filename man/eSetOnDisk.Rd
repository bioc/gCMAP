\name{eSetOnDisk}
\alias{eSetOnDisk}

\title{
A function to store the assayData of an eSet object as BigMatrix files
on disk.
}
\description{
This function accepts and eSet object and generates separate file-backed
BigMatrix objects for each assayDataElement. Data is only loaded into
memory upon subsetting, allowing the retrieval of selected data without
loading the (potentially large) object into memory in full.
}
\usage{
eSetOnDisk(eset, out.file = NULL)
}

\arguments{
  \item{eset}{
An object inheriting from eSet.
}
  \item{out.file}{
The path and basename of the output file. Three files will be generated
  for each eSet assayDataElement, identified by extending 'out.file' by
  suffices.
}
}

\value{
An object of the same class as 'eset', with BigMatrix elements in the
assayData slot.
}
\author{
Thomas Sandmann
}

\note{
  Please see the BigMatrix package for more details on BigMatrix objects.
}

\seealso{
  \code{\link{BigMatrix}}
  \code{\link{eSet}}
  \code{\link{memorize}}
}
\examples{

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
assayDataElement( on.disk, "exprs")[1:10,1:10] ## loads subset only
dim( assayDataElement( on.disk, "exprs")[,] ) ## retrieves full matrix

## convert back to standard in-memory ExpressionSet
in.memory <- memorize( on.disk )
class( assayDataElement( in.memory, "exprs" ) ) ## standard matrix object

## remove tempfiles generated for this example
unlink(paste(storage.file,"*", sep=""))  
}
