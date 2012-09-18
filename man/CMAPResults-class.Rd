\name{CMAPResults-class}
\Rdversion{1.1}
\docType{class}

\alias{CMAPResults}

\alias{cmapTable}
\alias{docs}
\alias{docs<-}
\alias{effect}
\alias{effect<-}
\alias{errors}
\alias{errors<-}
\alias{labels}
\alias{labels<-}
\alias{nFound}
\alias{nFound<-}
\alias{nSet}
\alias{nSet<-}
\alias{padj}
\alias{padj<-}
\alias{pval}
\alias{pval<-}
\alias{set}
\alias{set<-}
\alias{trend}
\alias{trend<-}
\alias{zscores}
\alias{geneScores}

\alias{CMAPResults-class}
\alias{cmapTable,CMAPResults-method}
\alias{initialize,CMAPResults-method}
\alias{docs,CMAPResults-method}
\alias{docs<-,CMAPResults-method}
\alias{effect,CMAPResults-method}
\alias{effect<-,CMAPResults-method}
\alias{errors,CMAPResults-method}
\alias{errors<-,CMAPResults-method}
\alias{labels,CMAPResults-method}
\alias{labels<-,CMAPResults-method}
\alias{nFound,CMAPResults-method}
\alias{nFound<-,CMAPResults-method}
\alias{nSet,CMAPResults-method}
\alias{nSet<-,CMAPResults-method}
\alias{padj,CMAPResults-method}
\alias{padj<-,CMAPResults-method}
\alias{plot,CMAPResults,ANY-method}
\alias{pval,CMAPResults-method}
\alias{pval<-,CMAPResults-method}
\alias{set,CMAPResults-method}
\alias{set<-,CMAPResults-method}
\alias{show,CMAPResults-method}
\alias{trend,CMAPResults-method}
\alias{trend<-,CMAPResults-method}
\alias{zscores,CMAPResults-method}
\alias{geneScores,CMAPResults-method}
\title{Class \code{"CMAPResults"}}

\description{ This class serves as a
container for the output of different gene-set enrichment analysis
methods. It directly extends the \code{\linkS4class{AnnotatedDataFrame}}
class by adding two additional slots ('docs' and 'errors' to store
information about the analysis raun. Data for each queried gene set are
stored in the 'data' slot of the
\code{\linkS4class{AnnotatedDataFrame}}. Additional information about
the data columns, e.g. the definition of 'effect' for the chosen
analysis method, is available in the varMetadata slot of the
AnnotatedDataFrame and can also be accessed through the 'labels'
accessor function.

}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("CMAPResults", ...)}.
  CMAPResults objects are usually created as output by gene set enrichment analysis methods.

}
\section{Slots}{
  \describe{
   \item{\code{data}:}{A data.frame containing results for different
     gene sets (rows), with method-specific output stored in the
     columns.}
   \item{\code{dimLabels}:}{A character vector of length 2 that provides
     labels for the rows and columns.}
   \item{\code{varMetadata}:}{A data.frame with the number of rows equal
     to the number of columns in 'data' and at least one column, named
     'labelDescription', containing additional information about each
     result column.}
   \item{\code{.__classVersion__}:}{A 'Versions' object describing the R
     and Biobase version numbers used to created the instance. Intended
     for developer use.}
    \item{\code{docs}:}{Object of class \code{"character"} Additional information about the analysis run, usually populated by the gene set enrichment method.
    }
    \item{\code{errors}:}{Object of class \code{"list"} Intented for
      warnings or error messages associated with the results. 
    }
    }
}
\section{Methods}{
  \describe{
    \item{cmapTable}{\code{signature(object = "CMAPResults")}:
      Returns data and labels stored in the 'table'
    AnnotatedDataFrame. If no additional parameters are supplied, this
    method is synonymous with \code{pData(object@table)}.

    Optional parameters:
    \code{n} (integer): the number of rows to return.
    \code{columns} (character): indicating which columns of the 'table'
      slot to include in the output.
    }
    
    \item{docs}{\code{signature(object = "CMAPResults")}:
      Accessor method for the 'docs' slot.
    }
    \item{docs<-}{\code{signature(x = "CMAPResults")}: 
      Replacement method for the 'docs' slot.
    }
    \item{effect}{\code{signature(object = "CMAPResults")}: 
      Accessor method for the 'effect' column of the 'table' slot.
    }
    \item{effect<-}{\code{signature(x = "CMAPResults")}:
      Replacement method for the 'effect' column of the 'table' slot.}
    \item{errors}{\code{signature(object = "CMAPResults")}:
      Accessor method for the 'docs' slot.
    }
    \item{errors<-}{\code{signature(x = "CMAPResults")}:
      Replacement method for the 'docs' slot.}
    \item{labels}{\code{signature(object = "CMAPResults")}:
      Returns information about the data columns of the 'table'
    slot. Synonymous with \code{varMetadata(object@table)}.
    }
    \item{labels<-}{\code{signature(x = "CMAPResults")}:
      Replacement method for the varMetdata slot of the 'table' AnnotatedDataFrame. Replacement value must be a data.frame with as many rows as there are columns in 'table' and contain the column 'labelDescription'. See \code{\linkS4class{AnnotatedDataFrame}} for details.
    }
    \item{nFound}{\code{signature(object = "CMAPResults")}: 
      Accessor method for the 'nFound' column of the 'table' slot.
    }
    \item{nFound<-}{\code{signature(x = "CMAPResults")}:
      Replacement method for the 'nFound' column of the 'table' slot.
    }
    \item{nSet}{\code{signature(object = "CMAPResults")}: 
      Accessor method for the 'nSet' column of the 'table' slot.}
    
    \item{nSet<-}{\code{signature(x = "CMAPResults")}: 
      Replacement method for the 'nSet' column of the 'table' slot.}
    
    \item{padj}{\code{signature(object = "CMAPResults")}:
      Accessor method for the 'padj' column of the 'table' slot.}
    
    \item{padj<-}{\code{signature(x = "CMAPResults")}:
      Replacement method for the 'padj' column of the 'table' slot.}
    
    \item{plot}{\code{signature(x = "CMAPResults", y = "ANY")}: Returns
      a plot with three-columns: 1. the distribution of scores
      across all results in the CMAPResults object. 2. a heatmap of
      rank-ordered effect sizes, 3. a stripplot of ranks for selected
      gene sets (by default, the top 5 These plots are most meaningful
      when relatively large sets of gene sets were evaluated.
        \itemize{
    \item{strip.effect: }{String specifying the CMAPResults column to retrieve scores from. Default:"effect"}
    \item{strip.pval: }{String specifying the CMAPResults column to transform into unsigned z-scores. Only evaluated if 'density.effect' column is not present or is set to 'NULL'.Default:"padj"}
    \item{strip.cutoffs: }{Numeric vector of length 2. Scores between strip.cutoffs[1] and strip.cutoffs[2] will be displayed in strip.col[2]. Default:c(-2,2)}
    \item{strip.bounds: }{Numeric vector of length 2 specifying the end points of the color gradient. Scores < strip.cutoffs[1] or > strip.cutoffs[2] will not be distinguishable.Default:c(-4,4)}
    \item{strip.col: }{Vector of length 3, specifying the colors used in the heatmap strip: strip.col[1] = low scores, strip.col[2] = medium score (excluded from gradient), strip.col[3] = high scores. Default:c("blue","white","red")}
    \item{strip.anno: }{Character string identifying the CMAPResults columns used to group samples in the stripplot. For example, strip.anno="trend" plots all sets showing the same direction of regulation together. Default:"set"}
    \item{strip.subset: }{Index vector specifying a subset of rows to consider as input. Default:c(1:5)}
    \item{strip.labels: }{Character vector with alternative labels for the groups defined by the 'strip.anno' parameter. Default:NULL}
    \item{strip.layout: }{Vector specifying the fraction of the plot occoputied by the densityplot, the levelplot and the stripplot components. Defaultc(0.45, 0.1, 0.45)}
    \item{set.inf: }{Numerical replacing Inf/-Inf scores in the density plot (default:+/-20)}
    \item{\dots: }{Further arguments passed on to lattice plot functions.}
    } 
 }
      \item{pval}{\code{signature(object = "CMAPResults")}:
	Accessor method for the 'pval' column of the 'table' slot.}
      \item{pval<-}{\code{signature(x = "CMAPResults")}:
	Replacemant method for the 'pval' column of the 'table' slot}
      \item{set}{\code{signature(object = "CMAPResults")}:
	Accessor method for the 'set' column of the 'table' slot.}  
      \item{set<-}{\code{signature(x = "CMAPResults")}:
	Replacement method for the 'set' column of the 'table' slot.}  
      \item{show}{\code{signature(object = "CMAPResults")}:
	Returns a summary of the CMAPResult object, including the number rows
	in the 'table' slot and shows the top five results with the smallest p-values.}
      \item{trend}{\code{signature(object = "CMAPResults")}:
	Accessor method for the trend' column of the 'table' slot.}  
      \item{trend<-}{\code{signature(x = "CMAPResults")}:
	Replacement method for the geneScores' column of the 'table' slot.} 
      \item{geneScores}{\code{signature(object = "CMAPResults")}:
  Accessor method for the geneScores' column of the 'table' slot. When available, this column stores a list of matrices, one for each row of the CMAPResults object, with raw per-gene scores for all members of the gene set. While the 'show' method displays only a brief summary of the available data, the geneScores method retrieves the full list of score matrices.}  
      \item{zscores}{\code{signature(x = "CMAPResults")}:
	Transforms adjusted p-values stored in a CMAPResults into
	z-scores based on the standard normal distribution.}
    }
  }
  
  \details{
  The AnnotatedDataFrame 'table' is populated by gene set enrichment
  analysis methods as needed. Explicit accessor and replacement methods
  exist for the following columns:
  \itemize{
    \item{set: }{Identifiers of the tested gene sets (e.g. obtained from the sampleNames of an analyzed CMAPCollection object).}
    \item{trend: }{The direction of the detected effect, e.g. 'upregulated', 'overrepresented', etc.}
    \item{pval: }{The raw p-value of the observed effect. Default}
    \item{effect: }{The detected effect size, e.g. log odds ratio (returned by 'fisher_score) or summay t-statistic (returned by gsealm_jg_score), etc.}
    \item{nSet: }{The number genes in the query gene sets}
    \item{nFound: }{The number of query set genes detected in the target set, e.g. genes common to both sets. Default:NULL}
  }
  In addition, gene set annotations can be included in further columns of the 'table' AnnotatedDataFrame, e.g. retrieved from the phenoData slot of a \code{\linkS4class{CMAPCollection}}.
}
\author{
  Thomas Sandmann, sandmann.thomas@gene.com
}

\seealso{
  \code{\linkS4class{AnnotatedDataFrame}}
}

\examples{
## create random score profile
set.seed(123)
z <- rnorm(1000)
names(z) <- paste("g", 1:1000, sep="")

## generate random incidence matrix of gene sets
m <-replicate(1000, {
  s <- rep(0,1000)
  s[ sample(1:1000, 20)] <-  1
  s[ sample(1:1000, 20)] <- -1
  s
  })
dimnames(m) <- list(names(z), paste("set", 1:1000, sep=""))

## Set1 is up-regulated
z <- z + m[,1]*2

## create CMAPCollection
cmap <- CMAPCollection(m, signed=rep(TRUE,1000))

## gene-set enrichment test
res <- gsealm_jg_score(z, cmap)
class(res)
res

## overview plot
plot(res, strip.subset=1:5)

## rerun, this time store gene-level scores
res <- gsealm_jg_score(z, cmap, keep.scores=TRUE)
res
m <- geneScores(res)
m[["set1"]] ## scores for set1

## stripplot for set1, colored by annotated sign
if( require( lattice ) ){
colors <- ifelse( attr(m[["set1"]], "sign") == "up", "red", "blue")
stripplot(m[["set1"]], xlab="z-score", main="set1", pch=23, col=colors)
}

## See examples for method 'featureScores'
## for additional gene-level plots
}

\keyword{classes}
