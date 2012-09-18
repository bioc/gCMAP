\name{pairwise_compare}

\alias{pairwise_compare}
\alias{pairwise_compare_limma}

\title{Generate statistics associated with pairwise differential expression}

\description{ When provided with an \code{\linkS4class{ExpressionSet}},
comparisons are made between control and perturbation samples. }

\usage{ 
pairwise_compare(eset, control_perturb_col = "cmap", control="control", perturb="perturbation")
pairwise_compare_limma(eset, control_perturb_col = "cmap",
control="control", perturb="perturbation", limma.index=2)
}

\arguments{
  
  \item{eset}{ \code{\linkS4class{ExpressionSet}} with all array data
  for a single instance, plus metadata on which arrays are perturbation
  and control. }

  \item{control_perturb_col}{ Column name in \code{phenoData} of
    \code{eset} where control/perturbation designations are stored. }

  \item{control}{ String designating control samples in the
    \code{control_perturb_col} column. }

  \item{perturb}{ String designating perturbation samples in the
    \code{control_perturb_col} column. }
  
  \item{limma.index}{ Integer specifying the index of the
    parameter estimate for which we to extract t and other
    statistics. The default corresponds to a two-class comparison
    with the standard parameterization. The function assumes that
    there was no missing data, so that test for all genes were
    performed on the same sample size.}
  
}

\note{The pairwise_compare functions returns p-values from a standard
  t-test. The pairwise_compare_limma functions uses the \pkg{limma}
  package instead to perform a moderated t-test.}

\value{
  
  The function returns a data frame with the following columns:

  \item{log_fc}{ Log fold-change between perturbed and control data. (A
    positive value denotes higher expression in the perturbed samples.) }

  \item{z}{ When at least one condition has two or more samples, the
    pairwise_compare_limma functions uses \code{\link{lmFit}},
    \code{eBayes} and \code{topTable} to compare the two classes and
    compute an (uncorrected) limma p-value. The pairwise_compare
    functions performs a standard t-test instead. For ease of comparison
    across instances with different numbers of samples, either p-value is
    converted to the standard normal scale. The result is reported
    here. As for \code{fc}, positive values denote higher expression in
    perturbed samples.  }

  \item{p}{ When at least one condition has two or more samples, the
    two-tailed standard (pairwise_compare) or \pkg{limma} p-value
    (pairwise_compare_limma), as computed by \code{\link{eBayes}}. Note
    that this p-value can also be computed from \code{z}, via
    \code{\link{pnorm}} (doubling for two tails). }
  
}
