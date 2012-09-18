\name{gCMAPData}
\alias{gCMAPData}

\docType{data}
\title{
Example \code{\linkS4class{NChannelSet}}
}

\description{
The gCMAPData object is an NChannelSet object with 3 samples x 1000
features x 3 channels (p-value, z-score and log_fc).
}

\usage{
data(gCMAPData)
}

\examples{
data(gCMAPData)

gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)

## comparison with a single user-provided profile of z-scores
profile <- assayDataElement(gCMAPData, "z")[,1]
gsealm_jg_score(profile, gene.set.collection)
}

\keyword{datasets}
