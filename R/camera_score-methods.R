
setMethod(
          "camera_score",
          signature( experiment = "eSet", sets = "CMAPCollection"),
          function( experiment, sets, predictor=NULL, design.matrix=NULL, element="exprs", keep.scores=FALSE, ... ) {

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
            indices <- geneIndex(sets, featureNames(experiment), remove=FALSE)
            
            ## run camera
            scores <- camera(indices, assayDataElement(experiment, element), design.matrix, ...)
            
            ## store raw per-gene expression scores as matrix
            if( keep.scores == TRUE ) {
              gene.scores <- featureScores(sets, experiment, element=element )
              gene.scores <- I(rep(gene.scores, ncol(scores)-2))
            } else { 
              gene.scores <- NA
           }

            ## store results
            res <- CMAPResults(
                               data=data.frame(
                                 set = rep(sampleNames(sets), ncol(scores)-2),
                                 trend = rep(colnames(scores)[-c(1,2)], each=nrow(scores)),
                                 pval = as.vector(as.vector(scores[,-c(1,2)])),
                                 padj = p.adjust( as.vector(as.vector(scores[,-c(1,2)])), method="BH"),
                                 nSet = rep( Matrix::colSums( abs( members (sets) ) ), ncol(scores)-2),
                                 nFound = rep( as.vector(scores[,1]), ncol(scores)-2),
                                 geneScores = gene.scores,
                                 pData(sets)[match(rep(sampleNames(sets),
                                                       ncol(scores)-2), sampleNames(sets)),,drop=FALSE ]),
                               docs ="\n Results were obtained with the 'camera' function from the 'limma' package.
 P-values were adjusted with the 'p-adjust' function using method 'BH'.")
            
            varMetadata(res)$labelDescription <-
              c("SetName",
                "Direction",
                "p-value",
                "adjusted p-value (BH)",
                "Number of genes annotated in the query set",
                "Number of query genes found in the dataset",
                "Per-gene raw expression scores",                                               
                colnames(pData(sets)))
            res
          }
          )

setMethod(
          "camera_score",
          signature( experiment = "matrix", sets = "CMAPCollection"),
          function( experiment, sets, element = "exprs", ...) {
            camera_score( ExpressionSet(experiment), sets, ...) 
          }
          )

setMethod(
          "camera_score",
          signature(experiment = "matrix",  sets = "GeneSet"),
          function( experiment, sets, ...) {
            camera_score(ExpressionSet(experiment), as(sets, "CMAPCollection"), ...)
          }
          )

setMethod(
          "camera_score",
          signature( experiment = "eSet", sets = "GeneSet"),
          function( experiment, sets, element = "exprs", ...) {
            camera_score( experiment, as(sets, "CMAPCollection"),  ...)
          }
          )


setMethod(
          "camera_score",
          signature( experiment = "matrix", sets = "GeneSetCollection"),
          function( experiment, sets, ...) {
            camera_score(ExpressionSet(experiment), as(sets, "CMAPCollection"),  ...) 
          }
          )

setMethod(
          "camera_score",
          signature( experiment = "eSet", sets = "GeneSetCollection"),
          function( experiment, sets, element = "exprs", ...) {
            camera_score(experiment, as(sets, "CMAPCollection"),  ...) 
          }
          )
