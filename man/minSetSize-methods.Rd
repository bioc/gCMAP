\name{minSetSize-methods}
\docType{methods}
\alias{minSetSize}
\alias{minSetSize,CMAPCollection-method}

\title{
GeneSetCollection length filtering
}

\description{
  This function filteres a GeneSetCollection by removing all contained
  GeneSets that do not include at least the user-specified number of
  genes also found in the user-specified universe.
}

\usage{
\S4method{minSetSize}{CMAPCollection}(sets, universe=NULL, min.members = 5)
}

\arguments{
  \item{sets}{
    A \code{\linkS4class{CMAPCollection}} object.
  }
  \item{universe}{
    Optional character vector of gene identfiers to be considered as the
    universe. Only geneIds included in the universe will count toward the
    gene set membership counts. If 'NULL' (default), all featureNames of the
    CMAPCollection will be considered.
  }
  \item{min.members}{
    Number of genes (in the univerese) a gene set needs to contain to be retained.
  }
}

\value{
  A CMAPCollection with all gene sets containing more than the specified
  number of members.
}

\author{
Thomas Sandmann
}

\examples{
data(gCMAPData)
gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)

minSetSize(gene.set.collection, min.members=100)

}
