\name{fisher_score-methods}
\docType{methods}

\alias{fisher_score-methods}
\alias{fisher_score}
\alias{fisher_score,CMAPCollection,CMAPCollection,character-method}
\alias{fisher_score,CMAPCollection,NChannelSet,character-method}
\alias{fisher_score,SignedGeneSet,CMAPCollection,character-method}
\alias{fisher_score,SignedGeneSet,NChannelSet,character-method}
\alias{fisher_score,GeneSet,CMAPCollection,character-method}
\alias{fisher_score,GeneSet,NChannelSet,character-method}
\alias{fisher_score,GeneSetCollection,CMAPCollection,character-method}
\alias{fisher_score,GeneSetCollection,NChannelSet,character-method}
\alias{fisher_score,GeneSetCollection,GeneSetCollection,character-method}
\alias{fisher_score,CMAPCollection,GeneSetCollection,character-method}
\alias{fisher_score,GeneSet,GeneSet,character-method}
\alias{fisher_score,GeneSet,GeneSetCollection,character-method}

\title{Hypergeometric probability of gene set enrichment}

\description{ A method for computing enrichment probilities based on the
  hypergeometric distribution. This method performs an over-representation analysis by generating 2x2 incidence matrices for gene sets provided as 'query' and 'sets' as GeneSet, SignedGeneSet, GeneSetCollection or CMAPCollection objects. If 'sets' is an NChannelSet object with quantitative data, gene sets are induced on the fly from the channel specified by the 'element' parameter.}

\section{Methods}{
  \describe{
    
\item{\code{signature(query = "CMAPCollection", sets = "CMAPCollection", universe = "character")}}{
}

\item{\code{signature(query = "CMAPCollection", sets = "NChannelSet", universe = "character")}}{
}

\item{\code{signature(query = "SignedGeneSet", sets = "CMAPCollection", universe = "character")}}{
}

\item{\code{signature(query = "SignedGeneSet", sets = "NChannelSet", universe = "character")}}{
}

\item{\code{signature(query = "GeneSet", sets = "CMAPCollection", universe = "character")}}{
}

\item{\code{signature(query = "GeneSet", sets = "NChannelSet", universe = "character")}}{
}

\item{\code{signature(query = "GeneSetCollection", sets = "CMAPCollection", universe = "character")}}{
}

\item{\code{signature(query = "GeneSetCollection", sets = "NChannelSet", universe = "character")}}{
}

\item{\code{signature(query = "GeneSet", sets = "GeneSetCollection", universe = "character")}}{
}

\item{\code{signature(query = "CMAPCollection", sets = "GeneSetCollection", universe = "character")}}{
}

\item{\code{signature(query = "GeneSetCollection", sets = "GeneSetCollection", universe = "character")}}{
}

\item{\code{signature(query = "GeneSet", sets = "GeneSet", universe = "character")}}{
}    
}
}

\arguments{
  
  \item{query}{A \code{\linkS4class{CMAPCollection}},
    \code{\linkS4class{GeneSet}}, or
    \code{\linkS4class{GeneSetCollection}} object containing the 'query'
    gene sets to compare against the 'sets'}
  
  \item{sets}{ A \code{\linkS4class{CMAPCollection}},
    \code{\linkS4class{GeneSetCollection}} or
    \code{\linkS4class{GeneSet}} object}
  
  \item{universe}{ A character string of gene ids for
    all genes that could potentially be of interest, e.g. all genes
    represented on a microarray, all annotated genes, etc.}
  
  \item{keep.scores}{ Logical: store the identifiers for the genes
    detected in 'query' and 'sets' ? (Default: FALSE)
    The size of the generated CMAPResults object increases with the
    number of contained gene sets. For very large collections,
    setting this parameter to 'TRUE' may require large amounts
    of memory.}
  \item{element}{ A character string corresponding to the assayDataElementName of the NChannelSet object to be thresholded on the fly with the \code{\link{induceCMAPCollection}}.}
  \item{lower}{ The lower threshold for the \code{\link{induceCMAPCollection}}.}
  \item{higher}{ The 'higher' threshold for the \code{\link{induceCMAPCollection}}.}
  \item{min.set.size}{ Number of genes a gene set induced by \code{\link{induceCMAPCollection}} needs to contain to be included in the analysis (Default:5).}
  \item{\dots}{ Additional arguments passed to downstream methods.}
  
}

\keyword{methods}

\value{

  A CMAPResults object 
}

\note{p-values are corrected for multiple testing separately for each
  query set, but not across multiple queries.}

\seealso{ \code{\link{fisher.test}}}

\examples{
data(gCMAPData)

gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)

## compare all gene sets in the gene.set.collection to each other
universe = featureNames(gCMAPData)
fisher_score(gene.set.collection, gene.set.collection, universe = universe)
}
