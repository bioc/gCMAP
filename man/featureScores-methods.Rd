\name{featureScores-methods}
\docType{methods}
\alias{featureScores-methods}
\alias{featureScores}

\alias{featureScores,CMAPCollection,eSet-method}
\alias{featureScores,CMAPCollection,matrix-method}
\alias{featureScores,eSet,CMAPCollection-method}
\alias{featureScores,matrix,CMAPCollection-method}
\alias{featureScores,CMAPCollection,numeric-method}
\alias{featureScores,numeric,CMAPCollection-method}

\title{ Methods to obtain scores for CMAPCollection gene sets from a matrix or eSet}
\description{
These methods extract the scores for CMAPCollection gene set members from eSet or matrix objects and return them as a list (argument 'query') of lists (argument 'dat') with score vectors. Argument order determines the organization of the list, e.g. if 'query' is a CMAPCollection, one list element is returned for each gene set, containing all score vectors for the respective set. If 'simplify' is set to TRUE, score vectors are combined and a list of matrices is returned instead. Score vectors and matrices carries an additional 'sign' attribute corresponding to the sign annotated in the CMAPCollection.
}
\section{Methods}{
\describe{

\item{\code{signature(query = "CMAPCollection", dat = "eSet")}}{\}}

\item{\code{signature(query = "CMAPCollection", dat = "matrix")}}{\}}

\item{\code{signature(query = "eSet", dat = "CMAPCollection")}}{\}}

\item{\code{signature(query = "matrix", dat = "CMAPCollection")}}{\}}

\item{\code{signature(query = "CMAPCollection", dat = "numeric")}}{\}}

\item{\code{signature(query = "numeric", dat = "CMAPCollection")}}{\}}

}}

\usage{

\S4method{featureScores}{CMAPCollection,eSet}(query, dat,element="z",simplify=TRUE)

\S4method{featureScores}{CMAPCollection,matrix}(query, dat,simplify=TRUE)

\S4method{featureScores}{eSet,CMAPCollection}(query, dat,element="z")

\S4method{featureScores}{matrix,CMAPCollection}(query, dat)

\S4method{featureScores}{CMAPCollection,numeric}(query, dat)

\S4method{featureScores}{numeric,CMAPCollection}(query, dat)

}
\arguments{
  
  \item{query}{A \code{\linkS4class{CMAPCollection}}, \code{\linkS4class{eSet}} or matrix.}

  \item{dat}{A \code{\linkS4class{CMAPCollection}}, \code{\linkS4class{eSet}} or matrix.}

  \item{element}{Character string specifying which assayDataElement to extract from eSet objects.}

  \item{simplify}{Logical: when possible, should score columns for each gene set collected in a matrix ? This optoin is only available if the results are organized by }

}
\keyword{methods}
\value{
 A nested list: one list element for each 'query', containing a list with score vectors for each 'dat'.
}

\author{
Thomas Sandmann
}

\examples{
  data(gCMAPData)
  ## generate CMAPCollection with two sets (drug1, drug2)
  sets <- induceCMAPCollection(gCMAPData, "z", higher=-2, lower=2)[,1:2]
  sampleNames(sets) <- c("set1", "set2")
  
  ## extract per-gene scores as matrices
  res <- featureScores(sets, gCMAPData)
  class(res) ## list
  names(res) ## one element per set
  class(res[["set1"]]) ## matrix
  dim(res[["set1"]])

  ## or as lists of score vectors
  res2 <- featureScores(sets, gCMAPData, simplify=FALSE)
  class(res2[["set1"]]) ## list
  length(res2[["set1"]])

  ## stripplot for set2, colored by annotated sign
  require( lattice )
  m <- res[["set2"]][,"drug2"]
  colors <- ifelse( attr(res[["set2"]], "sign") == "up", "red", "blue")
  stripplot(m, col=colors, xlab="z-score", main="set2, drug2")
    
  ## lattice plots
  if( require( reshape )) {
  y <- melt( res[["set1"]] )

  ## gene signs for "set1"
  signs <- attr(res[["set1"]], "sign")
            
  bwplot(value ~ X2, data=y, 
    xlab="Perturbation", ylab="z-score",
    main="set1",    
    panel = function(..., box.ratio) {
    panel.violin(..., col = "transparent",
    varwidth = FALSE, box.ratio = box.ratio)
    panel.xyplot(..., jitter.x=TRUE, fill = NULL,
    col=ifelse( signs == "up", "red", "blue"))
    })
}
}
