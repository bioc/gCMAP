\name{wilcox_score-methods}
\docType{methods}

\alias{wilcox_score}
\alias{wilcox_score-methods}
\alias{wilcox_score,matrix,CMAPCollection-method}
\alias{wilcox_score,numeric,CMAPCollection-method}
\alias{wilcox_score,eSet,CMAPCollection-method}
\alias{wilcox_score,matrix,GeneSet-method}
\alias{wilcox_score,numeric,GeneSet-method}
\alias{wilcox_score,matrix,GeneSetCollection-method}
\alias{wilcox_score,numeric,GeneSetCollection-method}
\alias{wilcox_score,eSet,GeneSetCollection-method}
\alias{wilcox_score,eSet,GeneSet-method}

\alias{wilcox_score-methods}
\alias{wilcox_score,CMAPCollection,matrix-method}
\alias{wilcox_score,CMAPCollection,numeric-method}
\alias{wilcox_score,CMAPCollection,eSet-method}
\alias{wilcox_score,GeneSet,matrix-method}
\alias{wilcox_score,GeneSet,numeric-method}
\alias{wilcox_score,GeneSetCollection,matrix-method}
\alias{wilcox_score,GeneSetCollection,numeric-method}
\alias{wilcox_score,GeneSetCollection,eSet-method}
\alias{wilcox_score,GeneSet,eSet-method}

\title{ Methods for Function \code{wilcox_score} in Package \pkg{gCMAP} }
\description{
  These methods provide a wrapper for the Mean-rank Gene Set Test function 
 \link[limma]{wilcoxGST}
 
  wilcox_score is a synonym for gst_score with ranks.only=TRUE. This
  test procedure was developed by Michaud et al (2008), who called it
  mean-rank gene-set enrichment.
 
}
\usage{

\S4method{wilcox_score}{matrix,CMAPCollection}(experiment, sets, adjust.ties=FALSE, keep.scores=FALSE, ...)

\S4method{wilcox_score}{numeric,CMAPCollection}(experiment, sets,...)

\S4method{wilcox_score}{eSet,CMAPCollection}(experiment, sets, element="z",...)

\S4method{wilcox_score}{matrix,GeneSet}(experiment, sets,...)

\S4method{wilcox_score}{numeric,GeneSet}(experiment, sets,...)

\S4method{wilcox_score}{eSet,GeneSet}(experiment, sets, element="z",...)

\S4method{wilcox_score}{matrix,GeneSetCollection}(experiment, sets,...)

\S4method{wilcox_score}{numeric,GeneSetCollection}(experiment, sets,...)

\S4method{wilcox_score}{eSet,GeneSetCollection}(experiment, sets, element="z",...)

\S4method{wilcox_score}{CMAPCollection,eSet}(experiment, sets, element="z",adjust.ties=FALSE, keep.scores=FALSE,...)

\S4method{wilcox_score}{CMAPCollection,numeric}(experiment, sets,...)

\S4method{wilcox_score}{CMAPCollection,matrix}(experiment, sets,...)

\S4method{wilcox_score}{GeneSet,numeric}(experiment, sets,...)

\S4method{wilcox_score}{GeneSet,matrix}(experiment, sets,...)

\S4method{wilcox_score}{GeneSet,eSet}(experiment, sets,element="z",...)

\S4method{wilcox_score}{GeneSetCollection,numeric}(experiment, sets,...)

\S4method{wilcox_score}{GeneSetCollection,matrix}(experiment, sets,...)

\S4method{wilcox_score}{GeneSetCollection,eSet}(experiment, sets,element="z",...)

}

\arguments{
  
  \item{sets}{ A \code{\linkS4class{CMAPCollection}},
  \code{\linkS4class{GeneSetCollection}} or \code{\linkS4class{GeneSet}}
  object containing gene sets, with which to query the \code{experiment}
  object. }
  
  \item{experiment}{ An \code{\linkS4class{eSet}} or matrix or vector with
    numeric data to compare the \code{query} object to.}
  
  \item{element}{ Character vector specifying which channel of an eSet
    to extract (defaults to "exprs", alternatives may be e.g. "z", etc.)}
  
  \item{\dots}{ Additional arguments passed on to downstream methods.}
  
  \item{adjust.ties}{ Logical: adjust Wilcox-Mann-Whitney statistic in
    the presence of ties ? (Default: FALSE)}
    
\item{keep.scores}{ Logical: keep gene-level scores for all gene sets (Default: FALSE) ? The size of the generated CMAPResults object increases with the number of contained gene sets. For very large collections, setting this parameter to 'TRUE' may require large amounts of memory.}
}
\keyword{methods}

\examples{
data(gCMAPData)
gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)

profile <- assayDataElement(gCMAPData[,1], "z")
## one profile versus three sets
wilcox_score(profile, gene.set.collection)

## three sets versus three profiles
wilcox_score(gene.set.collection, gCMAPData)
}
