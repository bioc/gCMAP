\name{generate_gCMAP_NChannelSet}

\alias{generate_gCMAP_NChannelSet}

\title{Generate a perturbation profile library from expression sets of control/treatment pairs}

\description{ 
When provided with a list of ExpressionSet or countDataSet objects, comparisons are made between
control and perturbation samples on a set basis. For countDataSets, a
moderated log2 fold change for each set is calculated after variance-stabilizing
transformation of the count data is performed globally across all
countDataSets in the list.
}

\usage{
generate_gCMAP_NChannelSet(
                           data.list,
                           uids=1:length(data.list),
                           sample.annotation=NULL,
                           platform.annotation="",
                           control_perturb_col="cmap",
                           control="control",
                           perturb="perturbation",
                           limma=TRUE,
                           limma.index=2,
                           big.matrix=NULL,
                           center.z="peak",
                           center.log_fc="peak"
                           )
}

\arguments{
  
	\item{data.list}{List of \code{\linkS4class{ExpressionSet}} or
	\code{\linkS4class{CountDataSet}} objects. 
	Each element includes all array / RNAseq data for a single instance, 
	plus metadata on which samples are perturbation and control.}

	\item{uids}{Vector of unique identifiers for the instances in \code{data.list}}

	\item{sample.annotation}{An optional \code{data.frame} of additional annotation for instances,
	each row corresponds to one instance, ordered to correspond with the \code{data.list}.
	This is not used for the control/perturbation comparisons, instead it is simply attached
	to the \code{\linkS4class{NChannelSet}} for future reference.}

	\item{platform.annotation}{The name of the platform as used by
	the annotation package. }

	\item{control_perturb_col}{See \code{\link{pairwise_compare}}.}

	\item{control}{See \code{\link{pairwise_compare}}.}

	\item{perturb}{See \code{\link{pairwise_compare}}.}
	
	\item{limma}{Use limma package to perform moderated t-tests
	  (Default: TRUE) instead of a standard t-test ?}
	  
	\item{limma.index}{ Integer specifying the index of the
	  parameter estimate for which we to extract t and other
	  statistics. The default corresponds to a two-class comparison
	  with the standard parameterization. The function assumes that
	  there was no missing data, so that test for all genes were
	  performed on the same sample size.}
	
	\item{big.matrix}{Character string providing the path and
	  filename to store the NChannelSets on disk instead of in
	  memory.
	  If 'NULL' (default), an NChannelSet is returned.
	  If not 'NULL', the BigMatrix package will create (or
	  overwrite !) three binary files for each channel of the
	  NChannelSet at the location provided as 'big.matrix',
	  distinguishing files for the different channes by their
	  suffices. To load the NChannelSet into a different R session,
	  the binary files must be accessible.}
  \item{center.z}{One of 'none', 'peak', 'mean', 'median', selecting whether / how to center the z-scores for each experiment. Option 'peak' (default) will center on the peak of the z-score kernel density. Options 'mean' and 'median' will center on their respective values instead.}
  \item{center.log_fc}{One of 'none', 'peak', 'mean', 'median', selecting whether / how to center the log2 fold-change distribution for each experiment. Option 'peak' (default) will center on the peak of the z-score kernel density. Options 'mean' and 'median' will center on their respective values instead.}
}

\value{

  The function returns an \code{\linkS4class{NChannelSet}} with one
  channel for each of the columns returned by
  \code{\link{pairwise_compare}}. This can be worked with directly (e.g,
  \code{assayData(obj)$z}), or specific channels can be converted to
  regular \code{\linkS4class{ExpressionSet}} objects (e.g.,\code{es <-
  channel(obj, "z")}). In the latter case, one would access \code{z} by
  \code{exprs(es)}.

}

\examples{
## list of ExpressionSets
data("sample.ExpressionSet") ## from Biobase

es.list <- list( sample.ExpressionSet[,1:4], sample.ExpressionSet[,5:8], sample.ExpressionSet[,9:12])
names(es.list) <- paste( "Instance", 1:3, sep=".")

de <- generate_gCMAP_NChannelSet(
    es.list,
    1:3,
    platform.annotation = annotation(es.list[[1]]),
    control_perturb_col="type",
    control="Control",
    perturb="Case") 

assayDataElementNames(de)
head( assayDataElement(de, "z") ) 

\dontrun{
## list of CountDataSets
set.seed( 123 )
cds.list <- lapply( 1:3, function(n) {
   cds <- makeExampleCountDataSet()
   featureNames(cds) <- paste("gene",1:10000, sep="_")
   cds
})

cde <- generate_gCMAP_NChannelSet(cds.list,
                           uids=1:3,
                           sample.annotation=NULL,
                           platform.annotation="Entrez",
                           control_perturb_col="condition",
                           control="A",
                           perturb="B")

assayDataElementNames(cde)
}
}
