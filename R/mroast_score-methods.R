
setMethod(
          "mroast_score",
          signature( experiment = "eSet", sets = "CMAPCollection"),
          function( experiment, sets, predictor=NULL, design.matrix=NULL,element="exprs",keep.scores=FALSE, ... ) {

            signed(sets) <- rep(FALSE, ncol(sets)) ## gene signs are not used

            if(! element %in% assayDataElementNames(experiment)) {
              stop( "AssayDataElement ", element," not found in the eSet.")
            }
            
            if( is.null( predictor ) & is.null( design.matrix ) ) {
              stop("Either the 'predictor' or 'design' parameter has to be provided to run mcamera.")

            } else if( ! is.null( predictor) & is.null( design.matrix )) {
              if( length( predictor) == 1 ) {
                if( predictor %in% colnames(pData(experiment))) {

                  message(paste("Design is based on pData column", predictor))
                  predictor <- as.integer(factor(pData(experiment)[,predictor]))
                }
              } else if( length( predictor ) == ncol(experiment)) {
                predictor <- as.integer(factor(predictor))-1
                message(paste("Design is based on 'predictor' parameter."))
              } else {
                stop("Parameter 'predictor' must either correspond to a pData column of an eSet or provide annotation for each sample column.")
              }
              
              ## create design matrix from predictor factor
              design.matrix <- cbind(Intercept=rep(1, dim(experiment)[2]),
                                     Group=as.integer(factor(predictor))-1)
            }
            
            ## match set and experiment identifiers
            indices <- geneIndex(sets, featureNames(experiment), remove.empty=FALSE)

            ## run mroast
            scores <- mroast(y=assayDataElement(experiment, element), index=indices, design=esign.matrix, ...)

            ## store raw per-gene expression scores
            if( keep.scores == TRUE ) {
            gene.scores <- featureScores(sets, experiment, element=element )
            } else{ 
              gene.scores <- rep( NA, nrow( scores))
            }
            
            ## store results
            res <- CMAPResults(
                               data=data.frame(
                                 set = row.names(scores),
                                 trend = scores$Direction,
                                 pval = scores$PValue,
                                 padj = scores$FDR,
                                 nSet = Matrix::colSums( abs( members (sets)[,row.names(scores)] ) ),
                                 geneScores = I(gene.scores),
                                 pData(sets)[row.names(scores),, drop=FALSE]),
                               docs = "\n All results, including adjusted p-values, were obtained \n with the 'mroast' function from the 'limma' package.."
            )
            
            varMetadata(res)$labelDescription <- 
              c("SetName",
                "Direction",
                "P-value",
                "False-discovery rate",
                "Number of genes annotated in the query set",
                "Per-gene raw expression scores",                                               
                colnames(pData(sets)))

            res
          }
          )

setMethod(
          "mroast_score",
          signature( experiment = "matrix", sets = "CMAPCollection"),
          function( experiment, sets,  ...) {
            mroast_score( ExpressionSet(experiment), sets, ...)
          }
          )


setMethod(
          "mroast_score",
          signature(experiment = "matrix",  sets = "GeneSet"),
          function( experiment, sets, ...) {
            mroast_score( ExpressionSet(experiment), as(sets, "CMAPCollection"), ...)
          }
          )

setMethod(
          "mroast_score",
          signature( experiment = "eSet", sets = "GeneSet"),
          function( experiment, sets, ...) {
            mroast_score( experiment, as(sets, "CMAPCollection"), ...)
          }
          )


setMethod(
          "mroast_score",
          signature( experiment = "matrix", sets = "GeneSetCollection"),
          function( experiment, sets, ...) {
            mroast_score(ExpressionSet( experiment ), as(sets, "CMAPCollection"), ...) 
          }
          )

setMethod(
          "mroast_score",
          signature( experiment = "eSet", sets = "GeneSetCollection"),
          function( experiment, sets, ...) {
            mroast_score(experiment, as(sets, "CMAPCollection"), ...) 
          }
          )
