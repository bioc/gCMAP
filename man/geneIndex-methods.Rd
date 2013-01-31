\name{geneIndex-methods}
\docType{methods}
\alias{geneIndex}
\alias{geneIndex-methods}
\alias{geneIndex,CMAPCollection,character-method}
\alias{geneIndex,GeneSet,character-method}
\alias{geneIndex,GeneSetCollection,character-method}

\title{  Methods for Function \code{geneIndex} in Package \pkg{gCMAP} ~~}

\description{
These methods match a character vector of gene ids to the members of a
GeneSet, GeneSetCollection or CMAPCollection and return the match indices.
}

\usage{
\S4method{geneIndex}{CMAPCollection,character}(gene.sets,gene.ids,remove.empty=TRUE)

\S4method{geneIndex}{GeneSetCollection,character}(gene.sets,gene.ids,remove.empty=TRUE)

\S4method{geneIndex}{GeneSet,character}(gene.sets,gene.ids,remove.empty=TRUE)
}

\arguments{
  
  \item{gene.sets}{ A \code{\linkS4class{CMAPCollection}},
    \code{\linkS4class{GeneSetCollection}} or
    \code{\linkS4class{GeneSet}} to match the 'gene.ids' against.}

  \item{gene.ids}{ A character string of gene identifiers whose position
    (if any) in the 'gene.sets' is to be determined.}

  \item{remove.empty}{ Logical parameter specifying whether gene sets
  without any matching gene.ids should be removed from the output.}

}

\value{
  An integer vector or (if a collection was searched) a list of
  integer vectors with the matching positions of gene.ids in the
  gene.sets.  }

\keyword{methods}

\examples{

## induce CMAPCollection
data(gCMAPData)
gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)

gene.ids <- geneIds(gene.set.collection[,2]) ## geneIds of the second set
geneIndex(gene.set.collection, gene.ids)
}
