\name{NEWS}
\title{NEWS file for the gCMAP}

\section{Changes in version 1.11.3}{
\itemize{
\item{BUGFIX: }{The wiki2cmap function now excludes empty files downloaded from wikipathways from parsing.}
\item{BUGFIX: }{The go2cmap function is now exported.}
\item{BUGFIX: }{The reactom2cmap function is currently defunct, because the reactome.db package is broken.}
}
}

\section{Changes in version 1.11.2}{
\itemize{
\item{CHANGE: }{Fixed typo in mroast_score method.}
}
}

\section{Changes in version 1.11.1}{
\itemize{
\item{BUGFIX: }{Fixed the \code{romer_score method} to accomodate the change in argument order in limma's \code{romer} function.}
}
}

\section{Changes in version 1.11.0}{
\itemize{
\item{CHANGE: }{Version bump for BioC update.}
}
}

\section{Changes in version 1.10.0}{
\itemize{
\item{CHANGE: }{Version bump for BioC update.}
}
}

\section{Changes in version 1.9.2}{
\itemize{
\item{BUGFIX: }{Updated geneIndex method to accept both symbols2indices
  and id2indices functions from limma, depending on the limma version.}
}
}

\section{Changes in version 1.8.1}{
\itemize{
\item{BUGFIX: }{Updated extraction of big.max objects while maintaining row.names.}
}
}

\section{Changes in version 1.7.2}{
  \itemize{
 \item{CHANGE: }{Helper functions KEGG2cmap, reactome2cmap, wiki2cmap and go2cmap have two additional parameters, min.size and max.size, restricting the output to gene sets containing members within the specified range.}
 \item{CHANGE: }{Method camera_score now outputs correlations in column 'effect'.}
 \item{BUGFIX: }{Correct numbers of annotated genes are now returned for each set.}
}
}

\section{Changes in version 1.7.1}{
  \itemize{
 \item{CHANGE: }{Vignettes are now stored in vignettes directory.}
}
}

\section{Changes in version 1.5.4}{
  \itemize{
 \item{NEW: }{The bigmemory and bigmemoryExtras packages are not optional, enabling use of gCMAP on Windows}
}
}
\section{Changes in version 1.5.3}{
  \itemize{
 \item{NEW: }{Updated titles of vignettes}
}
}
\section{Changes in version 1.5.2}{
  \itemize{
 \item{BUGFIX: }{Updated eSet construction in mapNmerge function}
}
}
\section{Changes in version 1.5.1}{
   \itemize{
\item{BUGFIX: }{fixed incorrect varMetaData element definition in .process_counts}
}
}
\section{Changes in version 1.3.5}{
   \itemize{
 \item{NEW: }{DESeq package is now listed as 'Suggested' instead of 'Depends'}
}
}

\section{Changes in version 1.3.4}{
  \itemize{
 \item{NEW: }{Removed lattice and reshape dependencies.}
 \item{NEW: }{Removed cmapQAPlot function.}
 \item{NEW: }{Added generate_gCMAP_NChannelSet's 'report.center' parameter.}
}
}
\section{Changes in version 1.3.3}{
  \itemize{
 \item{NEW: }{Experimented with cmapQAPlot function (removed again in version 1.3.4).}
}
}
\section{Changes in version 1.3.2}{
   \itemize{
\item{NEW: }{mgsa_score method wrapping the model-based gene set analysis approach available in the mgsa package.}
}
}

\section{Changes in version 1.3.1}{
  \itemize{
 \item{NEW: }{KEGG2cmap, reactome2cmap, wiki2cmap and go2cmap to functions to generate CMAPCollections from public annotation databases}
}
}

\section{Changes in version 1.3.1}{
  \itemize{
 \item{NEW: }{upIds and downIds methods for CMAPCollection objects.}
}
}

\section{Changes in version 1.3.0}{
  \itemize{
    \item{NEW: }{fisher_score method now accepts NChannelSets as 'sets' and thresholds a channel of choice on the fly.}
    \item{BUGFIX: }{induceCMAPCollection method was updated to fix a bug that prevented processing of certain single-sample eSets.}
    \item{NEW: }{parameters 'center.z' and 'center.log_fc' for the generate_gCMAP_NChannelSet function, centering each column of the score matrix on the column mean, median or density peak. (Internally, the NEW:  'center_eSet' function is called.)}
    \item{NEW: }{function 'center_eSet' to center columns of eSet channels on either their kernel density peak, their mean or their median.}
    \item{NEW: }{function mergeCMAPs to combined two eSets that share class, varLabels, annotation and assayDataElementNames.}
    \item{BUGFIX: }{Fixed missing gene names in featureScores method with CMAPCollections queries containing single-gene sets.}
    \item{NEW: }{splitPerturbations function to split ExpressionSet (e.g. downloaded from ArrayExpress) by annotations provided in the phenoData slot.}
    \item{NEW: }{eset_instances function to generate lists eSets from larger datasets by specifying an instance.matrix. }
    \item{NEW: }{setSizes methods for CMAPCollections returns a data.frame with the number of members for each gene set stored in the Collection. }
    \item{NEW: }{fisher_score method now reports z-scores in the effect column instead of log-odds-ratio. The latter are returned in the NEW: , separate LOR column. }
    \item{NEW: }{featureScores method now return list of matched identfiers when dispatched on two CMAPCollection. }
    \item{NEW: }{Added 'density.only', 'level.only' and 'strip.only' parameters to the CMAPResults plot method to selectively return only one plot component.}
    \item{IMPROVED}{ Reduced memory footprint of gsealm_jg_score and induceCMAPCollection methods.}
    \item{CHANGED: }{By default, 'removeShift' option for gsealm_jg_score and gsealm_score methods is now set to FALSE.}
    \item{BUGFIX: }{Corrected sign of effect size for fisher_score for perfect overlap between set and query.}
    \item{IMPROVED: }{Improved plot method for CMAPResults objects, lattice package is no longer required.}
  }
}	

\section{Changes in version 1.1.4}{
  \itemize{
    \item{BUGFIX: }{The gealm_score method is not defined for queries with missing values and will now return the corresponding error message.}
    \item{BUGFIX: }{GeneSetCollections containing only one set are processed correctly by gsealm_score.}
  }
}	

\section{Changes in version 1.1.3}{
  \itemize{
    \item{BUGFIX: }{Corrected sign of effect size for fisher_score for perfect overlap between set and query.}
  }
}	

\section{Changes in version 1.1.2}{
  \itemize{
    \item{BUGFIX: }{Fixed log-odds calculation in fisher_score: non-overlapping sets now get a log odds ratio of 0}
  }
}	

\section{Changes in version 1.1.1}{
  \itemize{
    \item{IMPROVED: }{Reduced memory footprint of induceCMAPCollection method for bigmemory-backed eSet objects.}
  }
}	

\section{Changes in version 1.1.0}{
  \itemize{
    \item{NEW: }{Version number was increased for the Bioconductor release.}
  }
}	

\section{Changes in version 0.99.6}{
  \itemize{
    \item{IMPROVED: }{Updated mapIdentifiers method for SignedGeneSet class now handles more different annotation types.}
  }
}	

\section{Changes in version 0.99.5}{
  \itemize{
    \item{NEW: }{The BigMatrix package changed its name to bigmemoryExtras. Currently, bigmemoryExtras only supports Mac and Unix platforms, thus gCMAP does not support Windows, either.}
  }
}	
\section{Changes in version 0.99.1}{	
 \itemize{		
    \item{NEW: }{First release of the gCMAP package to Bioconductor.}
  }
}
