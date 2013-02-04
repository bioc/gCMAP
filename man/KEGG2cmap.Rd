\name{KEGG2cmap}
\alias{KEGG2cmap}
\alias{reactome2cmap}
\alias{wiki2cmap}

\title{
Functions to generate species-specific CMAPCollections from Bioconductor KEGG.db, reactome.db annotation packages or the wikipathways \url{http://www.wikipathways.org/index.php/Download_Pathways} project.
}
\description{
These functions extract the gene sets defined for a species of interest from the KEGG or Reactome annotation packages or download the latest gene sets from the Wikipathways website and create a CMAPCollections. Wikipathways provides gene identifiers using different annotation sources; the wiki2cmap function only considers Entrez and Ensembl gene identifiers and return a CMAPCollection with Entrez ids.
}
\usage{
KEGG2cmap(species, annotation.package)
reactome2cmap(species, annotation.package)
wiki2cmap(species, annotation.package)
}

\arguments{
  \item{species}{
  Character, the identifier used by the KEGG or Reactome databases to identify the species of interest. For example, human categories are identified by 'Homo sapiens' in the Reactome and 'hsa' in the KEGG database.
  }
  \item{annotation.package}{
  Character, the name of the Bioconductor annotation package for the species of interest, e.g. 'org.Hs.eg.db' for human gene annotations.
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

reactome.hs <- reactome2cmap( species="Homo sapiens",
                              annotation.package="org.Hs.eg.db")
wikipathways.hs <- wiki2cmap( species="Homo sapiens",
                              annotation.package="org.Hs.eg.db")
}
}
