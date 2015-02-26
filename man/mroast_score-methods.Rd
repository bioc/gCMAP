\name{mroast_score-methods}
\docType{methods}

\alias{mroast_score}
\alias{mroast_score-methods}
\alias{mroast_score,matrix,CMAPCollection-method}
\alias{mroast_score,eSet,CMAPCollection-method}
\alias{mroast_score,eSet,GeneSetCollection-method}
\alias{mroast_score,matrix,GeneSetCollection-method}
\alias{mroast_score,eSet,GeneSet-method}
\alias{mroast_score,matrix,GeneSet-method}

\title{  Methods for Function \code{mroast_score} in Package
  \pkg{gCMAP} }

\description{
   These methods provide a wrapper for the Rotation Gene Set Tests function
   \link[limma]{mroast} 
   mroast tests whether any of the genes in the set are differentially
   expressed. 
}

\usage{
\S4method{mroast_score}{eSet,CMAPCollection}(
  experiment,
  sets,
  predictor=NULL,
  design.matrix=NULL,
  element="exprs",
  keep.scores=FALSE,
  ...)

\S4method{mroast_score}{matrix,CMAPCollection}(experiment, sets,...)

\S4method{mroast_score}{matrix,GeneSet}(experiment,sets,...)

\S4method{mroast_score}{eSet,GeneSet}(experiment, sets,...)

\S4method{mroast_score}{matrix,GeneSetCollection}(experiment,sets,...)

\S4method{mroast_score}{eSet,GeneSetCollection}(experiment,sets,...)
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
  
  Goeman, JJ, and Buhlmann, P (2007). Analyzing gene expression data in
  terms of gene sets: methodological issues. Bioinformatics 23, 980-987.

  Langsrud, O (2005). Rotation tests. Statistics and Computing 15,
  53-60.

  Phipson B, and Smyth GK (2010). Permutation P-values should never be
  zero: calculating exact P-values when permutations are randomly
  drawn. Statistical Applications in Genetics and Molecular Biology,
  Volume 9, Article 39.

  Routledge, RD (1994). Practicing safe statistics with the
  mid-p. Canadian Journal of Statistics 22, 103-110.

  Wu, D, Lim, E, Francois Vaillant, F, Asselin-Labat, M-L, Visvader, JE,
  and Smyth, GK (2010). ROAST: rotation gene set tests for complex
  microarray experiments. Bioinformatics 26,
  2176-2182. http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btq401?
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

res<- mroast_score(y, gene.set.collection, predictor = predictor,keep.scores=TRUE)
res 

## heatmap of expression scores for set1
set1.expr <- geneScores(res)[["set1"]]
heatmap(set1.expr, scale="none", Colv=NA, labCol=predictor,
        RowSideColors=ifelse( attr(set1.expr, "sign") == "up", "red", "blue"),
        margin=c(7,5))
legend(0.35,0,legend=c("up", "down"),
  fill=c("red", "blue"),
  title="Annotated sign",
  horiz=TRUE, xpd=TRUE)
}
