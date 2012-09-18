\name{fisher_score-methods}
\docType{methods}

\alias{fisher_score-methods}
\alias{fisher_score}
\alias{fisher_score,CMAPCollection,CMAPCollection,character-method}
\alias{fisher_score,SignedGeneSet,CMAPCollection,character-method}
\alias{fisher_score,GeneSet,CMAPCollection,character-method}
\alias{fisher_score,GeneSetCollection,CMAPCollection,character-method}
\alias{fisher_score,GeneSetCollection,GeneSetCollection,character-method}
\alias{fisher_score,CMAPCollection,GeneSetCollection,character-method}
\alias{fisher_score,GeneSet,GeneSet,character-method}
\alias{fisher_score,GeneSet,GeneSetCollection,character-method}

\title{Hypergeometric probability of gene set enrichment}

\description{ A method for computing enrichment probilities based on the
  hypergeometric distribution. }

\section{Methods}{
  \describe{
    
\item{\code{signature(query = "CMAPCollection", sets = "CMAPCollection", universe = "character")}}{
}

\item{\code{signature(query = "SignedGeneSet", sets = "CMAPCollection", universe = "character")}}{
}

\item{\code{signature(query = "GeneSet", sets = "CMAPCollection", universe = "character")}}{
}

\item{\code{signature(query = "GeneSetCollection", sets = "CMAPCollection", universe = "character")}}{
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
