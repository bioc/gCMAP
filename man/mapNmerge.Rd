\name{mapNmerge}
\alias{mapNmerge}

\title{
A function to map eSet featureNames and calculate summaries for
many-to-one mapping features
}
\description{
This function converts the featureNames of an eSet-derived object,
either by applying a user-specified translation function (e.g. to remove
pre- or suffices) or by refering to the annotation slot of the object to
locate the corresponding Bioconductor annotation package.

In cases where multiple features map to the same target identifier,
scores are summarized by applying 'summary.fun' (default: mean).
For eSet-like object with multiple assayDataElements, each element is
summarized separately.
}
\usage{
mapNmerge(
  eset,
  translation.fun = NULL,
  get = "ENTREZID",
  verbose = FALSE,
  summary.fun = function(x) mean(x, na.rm = TRUE))
}

\arguments{
  \item{eset}{
An eSet-like object.
}
  \item{translation.fun}{
A function that will be applied to the results of applying the 'featureNames'
  method to the eSet. If not 'NULL', this parameter takes precendence
  and the 'get' parameter will be ignored.
}
  \item{get}{
A character vector specifying the gene identifier universe to be
  retrieved from the Bioconductor annotation package.
}
  \item{verbose}{
Logical, should basid mapping statistics be returned ?
}
  \item{summary.fun}{
A function that will be applied to the scores after featureName
  mapping (default: mean). 
}
}
\note{
For large eSet objects, applying 'summary.fun' can be
  time-consuming. Other strategies, e.g. based on selecting a single
  probe for each gene based on cross-sample variability are available in
  the genefilter package.
}
\value{
  An eSet object with the same number of samples as the original and one
  row for each unique new featureName (after mapping & summary).
}
\author{
Thomas Sandmann, sandmann.thomas@gene.com
}

\examples{
\dontrun{
## requires hgu95av2.db annotation package

if( require( "hgu95av2.db" )) {
  data(sample.ExpressionSet) ## from Biobase
  dim(sample.ExpressionSet)
  head(featureNames(sample.ExpressionSet))
  entrez <- mapNmerge(sample.ExpressionSet)
  dim(entrez)
  head(featureNames(entrez))
}
}
}
