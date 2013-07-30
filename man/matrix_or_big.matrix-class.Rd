\name{matrix_or_big.matrix-class}
\Rdversion{1.1}
\docType{class}
\alias{matrix_or_big.matrix-class}

\title{Class \code{"matrix_or_big.matrix"}}
\description{
Union of  \code{base}'s 'matrix' and \code{bigmemory}'s 'big.matrix' objects. 
If the \code{bigmemory} package is unavailable, an empty, dummy \code{big.matrix} 
class is created first.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}

\author{
Thomas Sandmann, sandmann.thomas@gene.com
}
\note{
Helper class used to dispatch methods on either matrix or big.matrix objects.
}

\keyword{classes}
