\name{splitPerturbations}
\alias{splitPerturbations}

\title{
Function to split an ExpressionSet downloaded from ArrayExpress based on the experimental factors present in the phenoData slot
}
\description{
The ArrayExpress Bioconductor package provides access to microarray and RNAseq datasets from the EBI's ArrayExpress repository. Sample and experiment annotations are contained in the phenoData slot and can be used to automatically construct single-factor comparisons by subsetting the original ExpressionSet object. This function be used to automatically identify perturbation vs. control comparisons and splits the original dataset in to instance-level objects.
}

\details{
To identify 'perturbation versus control' comparisons, the user needs to inspect the phenoData slot (see examples) and choose the appropriate factor of interest as well as a term in this column that identifies experimental control samples. The 'controlled.factors' character string identifies additional factors (= columns in the phenoData slot), in which control samples must match their corresponding perturbation samples. 

For example (see example code section), an ExpressionSet may be annotated with two different annotated factors, Compound and Solven, corresponding to two columns in the pData slot.

The first column is of interest and therefore 'factor.of.interest should be set to 'Compound'. For each level of 'factor.of.interest' unique experimental conditions are identified based on the remaining pData columns. (To exclude columns, use the 'ignore.factors' parameter.) Separate ExpressionSet objects will be constructed for each unique experimental condition. 

To distinguish control samples from perturbations, the 'control' parameter needs to be provided. For example, if control samples in the 'factor.of.interest' column are annotated as 'vehicle', the 'control' parameter should be set to 'vehicle'. 

The second column in this example,'Solvent', contains additional information about the type of vehicle used for each experiment, e.g. DMSO, ethanol, etc. To ensure that each sample is matched to the correct control condition the 'controlled.factors' parameter is set to 'Vehicle' to include this annotation column when assigning control to perturbation samples. 

To consider all available annotation columns to match controls, the 'controlled.factors' parameters can be set to 'all' instead. (In this example, either setting the parameter to 'Vehicle' or 'all' yields identical results, as there is only one column in addition to the 'factor.of.interest'.)
}

\usage{
splitPerturbations(
  eset,
  control = "none",
  controlled.factors = NULL,
  factor.of.interest = "Compound",
  ignore.factors = NULL,
  cmap.column ="cmap",
  prefix = NULL)
}

\arguments{
  \item{eset}{
  An eSet object with experimental factors annotated in the phenoData slot. Experimental factors are identified by the prefix of the column name, specified in the 'prefix' parameter. In ExpressionSets obtained from the ArrayExpress respository experimental factors can be identified by their "^Factor" prefix.
}
  \item{control}{
 A character string identifying control samples in the 'factor.of.interest' column.
}
  \item{controlled.factors}{
 A character vector specifying which annotation columns should be matched to assign controls to perturbation samples. If the set to NULL shared controls are used for all perturbations. If set to 'all' all experimental factors must be identical between control and perturbation samples. Alternatively, individual factors and their combinations can be specified as a vector, e.g. as c("Vehice", "Time"). Column names can be abbreviated as long as they uniquely identify the pData column.
}
  \item{factor.of.interest}{
 Character string, the name of the pData column containing the factor of interest. Column name can be abbreviated aslong as it uniquely identifies the pData column.
}
  \item{ignore.factors}{
  A character vector with valid pData names specifying annotation columns to exclude. Column names can be abbreviated as long as they uniquely identify the pData column.
}
  \item{cmap.column}{
  Column name for an additional annotation column that will be added to all generated eSets. Used by the \code{\link{generate_gCMAP_NChannelSet}} function.
}
  \item{prefix}{
  String identifying pData columns with experimental factors. Setting the prefix to NULL will include all pData columns. 
}
}
\value{
A list of eSet objects, one for each unique experimental perturbation with perturbation and control samples.
}

\author{
Thomas Sandmann, sandmann.thomas@gene.com
}
\note{
Only experimental instances with valid controls will be returned.
}

\section{Warning }{
Annotations for experiments in public repositories are provided by the uploader and vary widely in quality and notation. This function is expected to handle experiments with clear perturbation / control annotations. Mileage on other datasets may vary.
}

\seealso{
\code{\link{generate_gCMAP_NChannelSet}}
\code{\link{annotate_eset_list}}
}
\examples{
require(Biobase)
data( sample.ExpressionSet )
head(pData( sample.ExpressionSet))
eset.list <- splitPerturbations( eset=sample.ExpressionSet,
                    factor.of.interest="type",
                    control="Control",
                    controlled.factors="sex",
                    ignore.factors="score",
                    prefix=""
)
length( eset.list )
## the first eset contains male Cases & controls
pData( eset.list[[1]])
## the second eset contains female Cases & controls
pData( eset.list[[2]])

## generate data.frame with sample annotations
annotate_eset_list( eset.list)
}
