\name{DESeq_nbinom}
\alias{.DESeq_nbinom}

\title{
Function to perform a DESeq analysis to detect differential expression
between perturbation and control groups.}

\description{ This function is a wrapper for a standard DESeq analyis
with two classes, perturbation and control, annotated in the
'conditions' column of the cds phenoData slot.  First , the size factors
are determined using default parameters. Next, a dispersion parameter is
estimated using the default (pooled) method. Finally, p-values are
estimated for differential expression between treatment and control
groups.  }

\usage{
.DESeq_nbinom(cds, control = "control", perturb = "perturbation",
try.hard = FALSE, control_perturb_col = "cmap",...)
}

\arguments{
  \item{cds}{
    A \code{\linkS4class{CountDataSet}} with perturbation and control
    samples identified in the pData condition slot.
}
\item{control}{
  Character string corresponding to the control factor level of the
  condition phenoData slot.
}
  \item{perturb}{
  Character string corresponding to the perturbation factor level of the
  condition phenoData slot.
}
\item{try.hard}{
    Logical parameter indicating the function's behavior in case the
  parametric (default) dispersion estimation fails. If FALSE (default),
  the function exits with an error. If TRUE, a non-parametric (loess)
  esimation is attempted instead. 
}
\item{control_perturb_col}{ Column name in \code{phenoData} of
  \code{cds} where control/perturbation designations are stored.
}
  
  \item{\dots}{
    Any additional parameters passed on to \link[DESeq]{estimateDispersions}
}
}
\value{
  See \link[DESeq]{nbinomTest} for details.
}
