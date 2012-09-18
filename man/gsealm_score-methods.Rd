\name{gsealm_score-methods}
\docType{methods}

\alias{gsealm_score}
\alias{gsealm_score-methods}
\alias{gsealm_score,eSet,GeneSet-method}
\alias{gsealm_score,ExpressionSet,GeneSet-method}
\alias{gsealm_score,eSet,GeneSetCollection-method}
\alias{gsealm_score,ExpressionSet,GeneSetCollection-method}
\alias{gsealm_score,eSet,CMAPCollection-method}
\alias{gsealm_score,ExpressionSet,CMAPCollection-method}
\alias{gsealm_score,matrix,CMAPCollection-method}
\alias{gsealm_score,matrix,GeneSet-method}
\alias{gsealm_score,matrix,GeneSetCollection-method}
\title{ ~~ Methods for Function \code{gsealm_score} in Package \pkg{gCMAPQuery} ~~}

\description{

  This method extends functions from the GSEAlm package to perform label-permutation based 
  differential expression analysis. In addition to gene set membership,
  information about the gene sign (up- or down-regulated) is taken into consideration. 
}

\usage{

\S4method{gsealm_score}{ExpressionSet,CMAPCollection}(query, set,
removeShift=TRUE, predictor=NULL, formula=NULL, nPerm=1000, parametric=FALSE,respect.sign=TRUE, keep.scores=FALSE,...)

\S4method{gsealm_score}{eSet,CMAPCollection}(query, set, element="exprs", ... )

\S4method{gsealm_score}{matrix,CMAPCollection}(query, set, predictor=NULL, ...)

\S4method{gsealm_score}{eSet,GeneSetCollection}(query, set, element="exprs",...)

\S4method{gsealm_score}{matrix,GeneSetCollection}(query, set, ...)

\S4method{gsealm_score}{ExpressionSet,GeneSet}(query, set,...)

\S4method{gsealm_score}{ExpressionSet,GeneSetCollection}(query, set,...)

\S4method{gsealm_score}{eSet,GeneSet}(query, set, element="exprs", ...)

\S4method{gsealm_score}{matrix,GeneSet}(query, set, ...)
}

\arguments{
  
  \item{query}{An \code{\linkS4class{ExpressionSet}} or matrix with
  normalized expression data.}

\item{set}{ A \code{\linkS4class{CMAPCollection}},
  \code{\linkS4class{GeneSetCollection}} or \code{\linkS4class{GeneSet}}
  object containing gene sets. Gene ids must match those of the 'query'}

\item{removeShift}{ logical: should normalization begin with a
  column-wise removal of the mean shift? }

\item{predictor}{ A character string identifying one column in the pData
  slot of a 'query' ExpressionSet from which to construct the formula
  for the linear model. Ignored if 'formula' is provided.}

\item{formula}{ The formula to be used in the linear model. See
  \code{\link{gsealmPerm}} for details.}

\item{nPerm}{ The number of sample-label permutations to perform.}

\item{parametric}{ Logical, if set to 'TRUE', no label permutations are
  performed. Instead, p-values are calculated based on a parametric
  approximation.}

\item{respect.sign}{ Logical, if set to 'FALSE', gene sign information
  is ignored, considering up- and down-regulated genes to be equal.}

\item{element}{ Character string specifying which element to extract
  when coercing an ExpressionSet from an eSet subclass.}

\item{keep.scores}{ Logical: keep gene-level scores for all gene sets (Default: FALSE) ? The size of the generated CMAPResults object increases with the number of contained gene sets. For very large collections, setting this parameter to 'TRUE' may require large amounts of memory.}

\item{\dots}{ Additional arguments passed on to downstream functions.}
}
\value{
  This method returns a CMAPResults object.
}

\keyword{methods}

\seealso{
  \code{\link{gsealmPerm}}
  \code{\link{lmPerGene}}
}

\examples{
data(gCMAPData)

## induce gene sets from a collection of z-scores
gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)
sampleNames(gene.set.collection) <- c("set1", "set2", "set3")

## random score matrix
y <- matrix(rnorm(1000*6),1000,6, dimnames=list(featureNames(gCMAPData), 1:6))

## set1 is differentially regulated
effect <- as.vector(members(gene.set.collection[,1]) * 2)
y[,4:6] <- y[,4:6] + effect

predictor <- c( rep("Control", 3), rep("Case", 3))

## run analysis and keep gene-level expression scores
res <- gsealm_score(y, gene.set.collection, predictor=predictor, nPerm=100, keep.scores=TRUE)
res

## heatmap of expression scores for set1
set1.expr <- geneScores(res)[["set1"]]
heatmap(set1.expr, scale="none", Colv=NA, labCol=predictor,
        RowSideColors=ifelse( attr(set1.expr, "sign") == "up", "red", "blue"),
        margin=c(7,5))
legend(0.35,0,legend=c("up", "down"), fill=c("red", "blue"), title="Annotated sign", horiz=TRUE, xpd=TRUE)
}
