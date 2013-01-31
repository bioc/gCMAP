\name{romer_score-methods}
\docType{methods}

\alias{romer_score}
\alias{romer_score-methods}
\alias{romer_score,matrix,CMAPCollection-method}
\alias{romer_score,eSet,CMAPCollection-method}
\alias{romer_score,eSet,GeneSetCollection-method}
\alias{romer_score,matrix,GeneSetCollection-method}
\alias{romer_score,eSet,GeneSet-method}
\alias{romer_score,matrix,GeneSet-method}

\title{  Methods for Function \code{romer_score} in Package
  \pkg{gCMAP} }

\description{
    These methods provide a wrapper for the Rotation Gene Set
   Enrichment Analysis function \link[limma]{romer}. Romer performes
   a competitive test in that the different gene sets are pitted against
   one another. Instead of permutation, it uses rotation, a parametric
   resampling method suitable for linear models (Langsrud, 2005).
}

\usage{
\S4method{romer_score}{eSet,CMAPCollection}(experiment,sets,predictor=NULL,
design.matrix=NULL, element="exprs", keep.scores=FALSE, ...)

\S4method{romer_score}{matrix,CMAPCollection}(experiment, sets,...)

\S4method{romer_score}{matrix,GeneSet}(experiment,sets,...)

\S4method{romer_score}{eSet,GeneSet}(experiment, sets,...)

\S4method{romer_score}{matrix,GeneSetCollection}(experiment,sets,...)

\S4method{romer_score}{eSet,GeneSetCollection}(experiment, sets,...)
}

\keyword{methods}

\arguments{
  
  \item{sets}{ A \code{\linkS4class{CMAPCollection}},
  \code{\linkS4class{GeneSetCollection}} or \code{\linkS4class{GeneSet}}
  object containing gene sets, with which to query the \code{experiment}
  object. }
  
  \item{experiment}{ An \code{\linkS4class{eSet}} or data matrix with
    numeric data to compare the \code{query} object to.}
  
  \item{predictor}{ A character vector or factor indicating the phenotypic
    class of the \code{experiment} data columns. Either the 'predictor' or
    'design' parameter must be supplied.}
  
  \item{design.matrix}{ A design matrix for the experiment. Either the 'predictor'
    or 'design' parameter must be supplied. If both are supplied, the
    'design' is used.}
  
  \item{element}{ Character vector specifying which channel of an eSet
    to extract (defaults to "exprs", alternatives may be e.g. "z", etc.)}
  
\item{keep.scores}{ Logical: keep gene-level scores for all gene sets (Default: FALSE) ? The size of the generated CMAPResults object increases with the number of contained gene sets. For very large collections, setting this parameter to 'TRUE' may require large amounts of memory.}

\item{\dots}{ Additional arguments passed to downstream methods.}
}

\value{
  A \code{\linkS4class{CMAPResults}} object.
}
\references{
  
  Langsrud, O, 2005. Rotation tests. Statistics and Computing 15, 53-60

  Doerum G, Snipen L, Solheim M, Saeboe S (2009). Rotation testing in
  gene set enrichment analysis for small direct comparison
  experiments. Stat Appl Genet Mol Biol, Article 34.

  Majewski, IJ, Ritchie, ME, Phipson, B, Corbin, J, Pakusch, M, Ebert,
  A, Busslinger, M, Koseki, H, Hu, Y, Smyth, GK, Alexander, WS, Hilton,
  DJ, and Blewitt, ME (2010). Opposing roles of polycomb repressive
  complexes in hematopoietic stem and progenitor cells. Blood, published
  online 5 May 2010. http://www.ncbi.nlm.nih.gov/pubmed/20445021

  Subramanian, A, Tamayo, P, Mootha, VK, Mukherjee, S, Ebert, BL,
  Gillette, MA, Paulovich, A, Pomeroy, SL, Golub, TR, Lander, ES and
  Mesirov JP, 2005. Gene set enrichment analysis: a knowledge-based
  approach for interpreting genome-wide expression profiles. Proc Natl
  Acad Sci U S A 102, 15545-15550
}

\examples{
data(gCMAPData)
gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)
sampleNames( gene.set.collection ) <- c("set1", "set2", "set3")

## random score matrix
y <- matrix(rnorm(1000*6),1000,6, dimnames=list(featureNames(gCMAPData), 1:6))

## set1 is differentially regulated
effect <- as.vector(members(gene.set.collection[,1]) * 2)
y[,4:6] <- y[,4:6] + effect

predictor <- c( rep("Control", 3), rep("Case", 3))

res <- romer_score(y, gene.set.collection, predictor = predictor, keep.scores=TRUE)
res 

## heatmap of expression scores for set1
set1.expr <- geneScores(res)[["set1"]]
heatmap(set1.expr, scale="none", Colv=NA, labCol=predictor,
        RowSideColors=ifelse( attr(set1.expr, "sign") == "up", "red", "blue"),
        margin=c(7,5))
legend(0.35,0,legend=c("up", "down"), fill=c("red", "blue"), title="Annotated sign", horiz=TRUE, xpd=TRUE)
}
