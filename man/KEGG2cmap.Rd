\name{KEGG2cmap}
\alias{KEGG2cmap}
\alias{reactome2cmap}
\alias{wiki2cmap}
\alias{go2cmap}

\title{
Functions to generate species-specific CMAPCollections from Bioconductor KEGG.db, reactome.db or GO.db annotation packages or the wikipathways \url{http://www.wikipathways.org/index.php/Download_Pathways} project.
}
\description{
These functions extract the gene sets defined for a species of interest from the KEGG,  Reactome or GO annotation packages or download the latest gene sets from the Wikipathways website and create a CMAPCollections. Wikipathways provides gene identifiers using different annotation sources; the wiki2cmap function only considers Entrez and Ensembl gene identifiers and return a CMAPCollection with Entrez ids. Please note that the GO graph structure will be used to associate every gene with all upstream annotation terms. The relationships between GO categories is not represented in the output CMAPCollection.
}
\note{
The reactome.db package is currently \href{https://support.bioconductor.org/p/63024/#64379}{broken}. Unfortunately, the \code{reactome2cmap} function is defunct until the annotation package has been fixed.
}

\usage{
KEGG2cmap(species, annotation.package, min.size=5, max.size=200)
reactome2cmap(species, annotation.package, min.size=5, max.size=200)
wiki2cmap(species, annotation.package, min.size=5, max.size=200 )
go2cmap(annotation.package="org.Hs.eg.db", ontology="BP", evidence=NULL, min.size=5, max.size=200)
}

\arguments{
  \item{species}{
  Character, the identifier used by the KEGG or Reactome databases to identify the species of interest. For example, human categories are identified by 'Homo sapiens' in the Reactome and 'hsa' in the KEGG database.
  }
  \item{annotation.package}{
  Character, the name of the Bioconductor annotation package for the species of interest, e.g. 'org.Hs.eg.db' for human gene annotations.
}
  \item{ontology}{
  One of 'BP', 'MF' and 'CC', identifying the 'biological process', 'molecular function' and 'cellular component' gene ontology domain of interest.
  }
  \item{evidence}{
  Character string identifying the evidence required for a GO annotation to be included in the result. If 'NULL' all annotated terms will be included.
  }
  \item{min.size}{
  Scalar integer, only gene sets with at least 'min.size' members will be returned.
  }
  \item{max.size}{
  Scalar integer, only gene sets with no more than 'max.size' members will be returned.
  }

}

\value{
A \code{\link{CMAPCollection}} object.
}

\author{
Thomas Sandmann
}

\seealso{
\code{\link{induceCMAPCollection}}
}
\examples{
\dontrun{
KEGG2.hs <- KEGG2cmap( species="hsa",
                       annotation.package="org.Hs.eg.db")
#reactome.hs <- reactome2cmap( species="Homo sapiens",
#                              annotation.package="org.Hs.eg.db")
wikipathways.hs <- wiki2cmap( species="Homo sapiens",
                              annotation.package="org.Hs.eg.db")
}
}
