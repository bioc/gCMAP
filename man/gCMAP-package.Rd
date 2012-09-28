\name{gCMAP-package}
\alias{gCMAP-package}
\alias{gCMAP}
\docType{package}
\title{
Tools for Connectivity Map-like analyses
~~ package title ~~
}
\description{

  This package provides tools for carrying out Connectivity Map-like analyses to
  compare gene sets to a collection of expression profiles, or more generally,
  to compare expression profiles to one another. The original Connectivity Map
  statistical approach is implemented, along with numerous improvements and
  extensions. 
~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab gCMAP\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2012-08-13\cr
Depends: \tab methods, GSEABase, DESeq\cr
Imports: \tab Biobase, limma (>= 3.12.1), GSEAlm, Category, bigmemoryExtras, Matrix, parallel, annotate, genefilter, RColorBrewer\cr
Suggests: \tab KEGG.db, reactome.db\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
Built: \tab R 2.15.0; ; 2012-08-16 21:08:11 UTC; unix\cr
}

Index:
\preformatted{
.DESeq_nbinom           Function to perform a DESeq analysis to detect
                        differential expression between perturbation
                        and control groups.
CMAPCollection-class    Class '"CMAPCollection"'
CMAPResults-class       Class '"CMAPResults"'
GeneSet                 Methods for 'GeneSet' and 'GeneColorSet'
SignedGeneSet           Construtor for SignedGeneSet
SignedGeneSet-class     Class '"SignedGeneSet"'
camera_score-methods    ~~ Methods for Function 'camera_score' in
                        Package 'gCMAPQuery' ~~
connectivity_score      Broad CMAP gene set enrichment metrics
fisher_score-methods    Hypergeometric probability of gene set
                        enrichment
gCMAPData               Example 'NChannelSet' and 'ExpressionSet' for
                        use with the gCMAPQuery package.
gCMAPQuery-package      Tools for Connectivity Map-like analyses
geneIndex-methods       ~~ Methods for Function 'geneIndex' in Package
                        'gCMAPQuery' ~~
generate_gCMAP_NChannelSet
                        Generate a perturbation profile library from
                        expression sets of control/treatment pairs
gsealm_jg_score-methods
                        Parametric test for testing normally
                        distributed scores for gene set enrichment
gsealm_score-methods    ~~ Methods for Function 'gsealm_score' in
                        Package 'gCMAPQuery' ~~
induceCMAPCollection-methods
                        ~~ Methods for Function 'induceCMAPCollection'
                        in Package 'gCMAPQuery' ~~
mapNmerge               A function to map eSet featureNames and
                        calculate summaries for many-to-one mapping
                        features
memorize                Create a new NChannelSet instance by selecting
                        specific channels and load BigMatrix assayData
                        into memory, if present
minSetSize              GeneSetCollection length filtering
mroast_score-methods    ~~ Methods for Function 'mroast_score' in
                        Package 'gCMAPQuery' ~~
pairwise_DESeq          Generate statistics associated with pairwise
                        differential expression from RNAseq count data
pairwise_compare        Generate statistics associated with pairwise
                        differential expression
romer_score-methods     ~~ Methods for Function 'romer_score' in
                        Package 'gCMAPQuery' ~~
signedRankSumTest       An implementation of the Wilcox rank sum test /
                        Mann-Whitney test that takes into account the
                        direction / sign of gene set members and
                        possibly the correlation between cases
wilcox_score-methods    ~~ Methods for Function 'wilcox_score' in
                        Package 'gCMAPQuery' ~~
zScores                 Function to calculate z-scores from p-values
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Thomas Sandmann <sandmann.thomas@gene.com>, Richard Bourgon <bourgon.richard@gene.com> and Sarah Kummerfeld <kummerfeld.sarah@gene.com>

Maintainer: Thomas Sandmann <sandmann.thomas@gene.com>
}
\keyword{ package }
