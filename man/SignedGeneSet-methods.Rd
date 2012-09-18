\name{SignedGeneSet}

\docType{methods}

\alias{SignedGeneSet}
\alias{SignedGeneSet,ANY-method}
\alias{SignedGeneSet,missing-method}

\title{Construtor for SignedGeneSet}

\description{
  The construtor is largely identical to \code{\link{GeneColorSet}}, but
  also handles an optional \code{geneSign} argument, which is an alias
  for \code{geneColor}.
}

\section{Methods}{

  \describe{

    \item{\code{signature(type = "ANY")}}{ Construtor which uses a
      template object. See all methods for the
      \code{\link{GeneColorSet}} constructor. If a \code{geneSign}
      argument is included by name, it will be used to populate the
      \code{geneColor} slot of the returned object. }

    \item{\code{signature(type = "missing")}}{ Basic method with no
      template object. }

  }

}

\keyword{methods}
