\name{gsealm_jg_score-methods}
\docType{methods}
\alias{gsealm_jg_score}

\alias{gsealm_jg_score-methods}
\alias{gsealm_jg_score,CMAPCollection,eSet-method}
\alias{gsealm_jg_score,CMAPCollection,matrix-method}
\alias{gsealm_jg_score,CMAPCollection,numeric-method}
\alias{gsealm_jg_score,eSet,CMAPCollection-method}
\alias{gsealm_jg_score,eSet,GeneSet-method}
\alias{gsealm_jg_score,eSet,GeneSetCollection-method}
\alias{gsealm_jg_score,GeneSet,eSet-method}
\alias{gsealm_jg_score,GeneSet,matrix-method}
\alias{gsealm_jg_score,GeneSet,numeric-method}
\alias{gsealm_jg_score,GeneSetCollection,eSet-method}
\alias{gsealm_jg_score,GeneSetCollection,matrix-method}
\alias{gsealm_jg_score,GeneSetCollection,numeric-method}
\alias{gsealm_jg_score,matrix_or_big.matrix,CMAPCollection-method}
\alias{gsealm_jg_score,matrix,GeneSet-method}
\alias{gsealm_jg_score,matrix,GeneSetCollection-method}
\alias{gsealm_jg_score,numeric,CMAPCollection-method}
\alias{gsealm_jg_score,numeric,GeneSet-method}
\alias{gsealm_jg_score,numeric,GeneSetCollection-method}

\title{Parametric test for testing normally distributed scores for gene
  set enrichment}

\description{
  This method implements the 'JG' summary method introduced by Oron et
  al, 2008, as a stand-alone method to query a set of
  normally-distributed scores (e.g. t-statistics or z-scores) with gene
  sets (and vice versa).

  Scores for gene-set members are summed, respecting their sign
  (up- or down-regulated) and the combined score is divided by the
  square-root of the number of set members.
  
  To fit linear models to an expression profiling experiment instead,
  please use the \code{gsealm_jg_score} method instead.
}

\section{Methods}{
  \describe{
    \item{\code{signature(query = "CMAPCollection", sets = "eSet")}}{
}

\item{\code{signature(query = "CMAPCollection", sets = "matrix")}}{
}

\item{\code{signature(query = "CMAPCollection", sets = "numeric")}}{
}

\item{\code{signature(query = "eSet", sets = "CMAPCollection")}}{
}

\item{\code{signature(query = "eSet", sets = "GeneSet")}}{
}

\item{\code{signature(query = "eSet", sets = "GeneSetCollection")}}{
}

\item{\code{signature(query = "GeneSet", sets = "eSet")}}{
}

\item{\code{signature(query = "GeneSet", sets = "matrix")}}{
}

\item{\code{signature(query = "GeneSet", sets = "numeric")}}{
}

\item{\code{signature(query = "GeneSetCollection", sets = "eSet")}}{
}

\item{\code{signature(query = "GeneSetCollection", sets = "matrix")}}{
}

\item{\code{signature(query = "GeneSetCollection", sets = "numeric")}}{
}

\item{\code{signature(query = "matrix_or_big.matrix", sets = "CMAPCollection")}}{
}

\item{\code{signature(query = "matrix", sets = "GeneSet")}}{
}

\item{\code{signature(query = "matrix", sets = "GeneSetCollection")}}{
}

\item{\code{signature(query = "numeric", sets = "CMAPCollection")}}{
}

\item{\code{signature(query = "numeric", sets = "GeneSet")}}{
}

\item{\code{signature(query = "numeric", sets = "GeneSetCollection")}}{
}
}}

\keyword{methods}

\arguments{

  \item{query}{ An \code{\linkS4class{eSet}},
    \code{\linkS4class{CMAPCollection}},
    \code{\linkS4class{GeneSetCollection}},\code{\linkS4class{GeneSet}},
    matrix or numeric vector with data and gene ids. If a matrix is
    provided, gene ids must be provided as row-names. If a vector is
    provided, gene ids must be provided as names.}
    
  \item{sets}{ See 'query'}
    
  \item{removeShift}{ Optional parameter indicating that the aggregated
  test statistic should be centered by subtracting the column means (default=FALSE).Note: this option is not
  available for analysis of big.matrix backed eSet objects.}

  \item{element}{ For \code{\linkS4class{eSet}} objects, which
    assayDataElement should be extracted ?
    (Default="exprs")}
 
\item{keep.scores}{ Logical: keep gene-level scores for all gene sets (Default: FALSE) ? The size of the generated CMAPResults object increases with the number of contained gene sets. For very large collections, setting this parameter to 'TRUE' may require large amounts of memory.}

\item{\dots}{ Additional arguments to be passed on to downstream
    methods.}
}

\value{
  A CMAPResults object or, in case of multi-dimensional queries, a list
  of CMAPResults objects.
}

\references{ Gene set enrichment analysis using linear models and
  diagnostics. Oron AP, Jiang Z, Gentleman R. Bioinformatics. 2008 Nov
  15;24(22):2586-91. Epub 2008 Sep 11.

  Extensions to gene set
  enrichment. Jiang Z, Gentleman R. Bioinformatics. 2007 Feb
  1;23(3):306-13. Epub 2006 Nov 24.  }

\examples{
data(gCMAPData)
gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)

## comparison with a single user-provided profile of z-scores
profile <- assayDataElement(gCMAPData, "z")[,1]
gsealm_jg_score(profile, gene.set.collection)

## comparison with of multiple profiles of z-scores to the CMAPCollection
res <- gsealm_jg_score(assayDataElement(gCMAPData, "z"), gene.set.collection)

## first CMAPResult object
res[[1]]

## adjusted p-values from all CMAPResult objects
sapply(res, padj)

## inverted query: CMAPCollection is compared to z-score profiles
gsealm_jg_score(gene.set.collection, assayDataElement(gCMAPData, "z"))[[1]]
}
