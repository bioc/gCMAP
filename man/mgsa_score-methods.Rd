\name{mgsa_score-methods}
\docType{methods}

\alias{mgsa_score-methods}
\alias{mgsa_score}
\alias{mgsa_score,GeneSet,CMAPCollection,character-method}
\alias{mgsa_score,GeneSet,NChannelSet,character-method}
\alias{mgsa_score,SignedGeneSet,CMAPCollection,character-method}
\alias{mgsa_score,SignedGeneSet,NChannelSet,character-method}
\alias{mgsa_score,GeneSetCollection,CMAPCollection,character-method}
\alias{mgsa_score,GeneSetCollection,NChannelSet,character-method}
\alias{mgsa_score,GeneSetCollection,GeneSetCollection,character-method}
\alias{mgsa_score,GeneSet,GeneSetCollection,character-method}
\alias{mgsa_score,GeneSet,GeneSet,character-method}
\alias{mgsa_score,CMAPCollection,CMAPCollection,character-method}
\alias{mgsa_score,CMAPCollection,GeneSetCollection,character-method}

\title{Model-based gene set analysis (MGSA) }

\description{ This method is a wrapper for the mgsa methods from the Bioconductor package mgsa, which must be available on the system for the methods to run. The model-based gene set analysis (MGSA) analyzes all categories at once by embedding them in a Bayesian network, naturally taking overlap between categories into account and avoiding the need for multiple testing correction. Please consult the mgsa help page for more details. }

\section{Methods}{
  \describe{
    
\item{\code{signature(query = "GeneSet", sets = "CMAPCollection", universe = "character")}}{
}

\item{\code{signature(query = "GeneSet", sets = "NChannelSet", universe = "character")}}{
}

\item{\code{signature(query = "SignedGeneSet", sets = "CMAPCollection", universe = "character")}}{
}

\item{\code{signature(query = "SignedGeneSet", sets = "NChannelSet", universe = "character")}}{
}

\item{\code{signature(query = "GeneSetCollection", sets = "CMAPCollection", universe = "character")}}{
}

\item{\code{signature(query = "GeneSetCollection", sets = "NChannelSet", universe = "character")}}{
}

\item{\code{signature(query = "GeneSetCollection", sets = "GeneSetCollection", universe = "character")}}{
}

\item{\code{signature(query = "GeneSet", sets = "GeneSetCollection", universe = "character")}}{
}

\item{\code{signature(query = "GeneSet", sets = "GeneSet", universe = "character")}}{
}

\item{\code{signature(query = "CMAPCollection", sets = "CMAPCollection", universe = "character")}}{
}
\item{\code{signature(query = "CMAPCollection", sets = "GeneSetCollection", universe = "character")}}{
}}
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
  \item{\dots}{ Additional arguments passed to mgsa function from the mgsa package, including the following:
   \itemize{
    \item{alpha: }{Grid of values for the parameter alpha. Values represent probabilities of false-positive events and hence must be in [0,1]. numeric}
    \item{beta: }{Grid of values for the parameter beta. Values represent probabilities of false-negative events and hence must be in [0,1]. numeric.}
    \item{p: }{Grid of values for the parameter p. Values represent probabilities of term activity and therefore must be in [0,1]. numeric.}
    \item{steps: }{The number of steps of each run of the MCMC sampler. integer of length 1. A recommended value is 1e6 or greater.}
    \item{restarts: }{The number of different runs of the MCMC sampler. integer of length 1. Must be greater or equal to 1. A recommended value is 5 or greater.}
    \item{threads: }{The number of threads that should be used for concurrent restarts. A value of 0 means to use all available cores. Defaults to 'getOption(mc.cores, default=0)', which will instruct mgsa to use all available cores.}
    }
  }
}

\keyword{methods}

\value{
  A CMAPResults object. The reported p-values represent '1-marginal posterior probability'. For the 'effect' column, the p-values have been transformed to z-scores using a standard normal distribution.
}

\note{ This Bayesian approach does not require any additional correction of p-values for multiple testing. For consistency, the returned CMAPResults object contains a padj column duplicating the content of the pval column.}

\seealso{ mgsa }

\examples{
if( is.element("mgsa", installed.packages()[,1])){
   suppressMessages( require( "mgsa", character.only = TRUE ))
   
   data(gCMAPData)
   gene.set.collection <- induceCMAPCollection(gCMAPData, "z", 
   higher=2, lower=-2)

    ## compare all gene sets in the gene.set.collection 
    ## to each other
    universe = featureNames(gCMAPData)
     mgsa_score(gene.set.collection, gene.set.collection, 
     universe = universe)
  }
}
