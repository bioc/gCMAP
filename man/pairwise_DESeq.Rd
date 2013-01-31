\name{pairwise_DESeq}
\alias{pairwise_DESeq}

\title{Generate statistics associated with pairwise differential
  expression from RNAseq count data}

\description{ When provided with an CountDataSet, comparisons are made between control and perturbation samples. }

\usage{ pairwise_DESeq(cds, vst, control_perturb_col = "condition",
control="control", perturb="perturbation", try.hard=FALSE) }

\arguments{
  
  \item{cds}{ CountDataSet with all count data
  for a single instance, plus metadata on which samples are perturbation
  and control. }

  \item{vst}{ Matrix of variance-stabilized count data that must include
  columns with colnames matching the sampleNames of the cds object. The
  vst matrix may contain additional columns / samples, which will be ignored. }

  \item{control_perturb_col}{ Column name in \code{phenoData} of
    \code{cds} where control/perturbation designations are stored. }

  \item{control}{ String designating control samples in the
    \code{control_perturb_col} column. }

  \item{perturb}{ String designating perturbation samples in the
    \code{control_perturb_col} column. }

  \item{try.hard}{ Logical parameter indicating how to proceed when
    DESeq's parametric estimation of the dispersion parameter
    fails. If set to FALSE (default), the function exits with an
    error. If set to TRUE, the function will try a non-parametric
    approach instead. } 
}

\value{
  
  The function returns a data frame with the following columns:

  \item{log_fc}{ Moderated log2 fold-change between perturbed and control data. (A
    positive value denotes higher expression in the perturbed samples.)
    The change was calculated from the (mean) counts after variance
    stabilizing transformation. Please consult the \pkg{DESeq} vignette for details on
    the transformation.}

  \item{z}{ For ease of comparison across instances with different
    numbers of samples, the (uncorrected) DESEq p-value is converted to
    the standard normal scale. The result is reported here. As for
    \code{log_fc}, positive values denote higher expression in perturbed
    samples.  }

  \item{p}{ p-value for differential expression calculated by the
    \code{nbinomTest} function from the \pkg{DESeq} package. In the absence of replicates, the dispersion parameter
    is estimated across all samples, ignoring the class labels, by using
    the blind method of the \code{estimateDispersions} function. When
    replicates are available, the pooled method is used instead. Note
    that this p-value can also be computed from \code{z}, via
    \code{\link{pnorm}} (doubling for two tails). }
  
}

\note{
To use this function, please install the suggested Bioconductor package 'DESeq'.
}