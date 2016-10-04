\name{eset_instances}
\alias{eset_instances}

\title{
A function to subset an eSet with expression data into smaller datasets, each corresponding to a single perturbation experiment.
}
\description{
This function takes two parameters, an eSet object (e.g. an ExpressionSet or CountDataSet) containing muliple samples, and a numerical matrix defining how these samples should be compared to investigate perturbations of interest. For each perturbation, a separate eSet object is generated, ready for analysis with the \code{\link{generate_gCMAP_NChannelSet}} function. Samples can be used in multiple instances, e.g. common controls can be specified in each column of the 'instance.matrix'.
}
\usage{
eset_instances(instance.matrix, eset, control_perturb_col = "cmap", 
control = "control", perturb = "perturbation")
}

\arguments{
  \item{instance.matrix}{
  A numeric matrix of -1 and 1's. Each columns defines a contrast of interest and indicates whether a sample (row) corresponds to a control sample (-1) or a perturbation sample (1). The row.names of the instance.matrix correspond to sampleNames of 'eset'. Entries other than -1 or 1 will be ignored.
  }
  \item{eset}{
  An eSet object to be subset into smaller datasets. The row.names of 'eset' must correspond to the row.names of the 'instance.matrix'.
  }
  \item{control_perturb_col}{
  Character, indicating which phenoData column to use to store 'control' and 'perturb' labels.
}
  \item{control}{
  Character, defining the label stored in each new eSet to indicate control samples.
}
  \item{perturb}{
  Character, defining the label stored in each new eSet to indicate perturbation samples.
}
}
\value{
A list of eSet objects, each corresponding to one instance defined by the columns of 'incidence.matrix'.
}

\author{
Thomas Sandmann, sandmann.thomas@gene.com
}
\note{
This function can be used to generate the 'eset.list' required for differential expression analyses with \code{\link{generate_gCMAP_NChannelSet}}.
}

\seealso{
\code{\link{generate_gCMAP_NChannelSet}}
}
\examples{
  library(Biobase)
  data(sample.ExpressionSet)
  
  ## contains Male/Female and Control/Case annotations
  pData( sample.ExpressionSet)
  
  ## separate analysis of Male/Female patients
  male <- ifelse( pData( sample.ExpressionSet )$type == "Control", -1, 1)
  male[which( pData( sample.ExpressionSet )$sex == "Female")] <- 0
  
  female <- ifelse( pData( sample.ExpressionSet )$type == "Control", -1, 1)
  female[which( pData( sample.ExpressionSet )$sex == "Male")]  <- 0
  
  instance.matrix <- cbind( male, female)
  row.names( instance.matrix ) <- sampleNames( sample.ExpressionSet )
  
  eset_instances( instance.matrix, sample.ExpressionSet)
}
