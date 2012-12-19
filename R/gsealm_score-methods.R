
setMethod(
          "gsealm_score",
          signature( query = "ExpressionSet", set="CMAPCollection"),
          function( query, set, removeShift=TRUE, predictor=NULL,
                   formula=NULL, nPerm=1000, parametric=FALSE, respect.sign=TRUE,
                   keep.scores=FALSE, ...) {

            ## create formula
            if( is.null( predictor ) & is.null( formula ) ) {
              stop("Either the 'predictor' or 'formula' parameter has to be provided to run gsealm_score")

            } else if( ! is.null( predictor) & is.null( formula )) {
              if (! predictor %in% colnames( pData( query ) ) ) {
                stop("Predictor parameter does not correspond to a phenoData column.")
              }
              
              formula= as.formula(paste("~", eval(predictor)))
            }
            ## create incidence matrix

            mat <- Matrix::t( members (set) )

            shared.genes <- intersect(featureNames(query), colnames(mat))
            mat <- try(
                       mat[,shared.genes, drop=FALSE]
                       )

            if( respect.sign == FALSE) {
              mat <- abs( mat )
              signed(set) <- rep(FALSE, ncol(set))
            }
            
            if (class(mat) == "try-error" | sum( abs ( mat) ) == 0) {
              stop( "None of the user-provided gene ids were found in the gene set collection.")
            }

            query <- query[shared.genes,]
            if( any( is.na( exprs( query )))){
              stop( "Query ExpressionSet contains NA values. Please replace (e.g. impute) these missing values or remove genes with missing values from the ExpressionSet/")
            }
            
            ## store raw per-gene expression scores
            if( keep.scores == TRUE) {
              gene.scores <- featureScores(set, query, element="exprs" )
              gene.scores <- I(gene.scores)

            } else {
              gene.scores <- NA
            }
            ## fit user-specified linear model
            lm <- lmPerGene(query, formula, na.rm=TRUE)
            notNAN <- sapply(lm$tstat[2,],function(x) !is.nan(x))
            
            GSNorm <- .gsNormalize( lm$tstat[2,,drop=FALSE][notNAN], mat[,notNAN, drop=FALSE], removeShift=removeShift)
            GSNorm <- as.numeric( GSNorm )
            
            if (parametric == TRUE) {
              bestPval <- 1-sapply( abs( GSNorm ), pnorm )
              
            } else {
              pVals  <- .gsealmPerm( query[notNAN,], formula, mat[,notNAN], nperm=nPerm)
              bestPval <- ifelse(GSNorm < 0, pVals[,"Lower"], pVals[,"Upper"]) ## only test in the direction of sign(t)
            }
            
            ## return different 'trend' descriptors for signed and unsigned sets
            trend <- sapply(seq( ncol(set) ), function( n ){
              if( signed(set)[n] == FALSE ) {
                ifelse(GSNorm[n] >=0, "up", "down")
              } else {
                ifelse(GSNorm[n] >=0, "correlated", "anticorrelated")
              }})
            
            if (parametric == FALSE ) {
              docs <- paste("\n GSEAlm analysis with formula ", paste(as.character(formula), collapse=""),
                            " using ", nPerm, " sample label permutations.\n P-values were adjusted with the 'p-adjust' function using method 'BH'."
                            , sep="")
            } else{
              docs <-"\n GSEALm analysis with parametric p-value estimate (JG-score).\n P-values were adjusted with the 'p-adjust' function using method 'BH'."
            }
            
            ## store results for differential expression analysis
            res <- CMAPResults(
                               data=data.frame(
                                 set = sampleNames(set),
                                 trend  = trend,
                                 pval   = bestPval,
                                 padj   = p.adjust( bestPval, method="BH"),
                                 effect = GSNorm,
                                 nSet = Matrix::colSums( abs( members (set) ) ),
                                 nFound = Matrix::rowSums( abs( mat ) ),
                                 geneScores = gene.scores,                                                          
                                 pData(set)),
                               docs = docs)
            varMetadata(res)$labelDescription <- 
              c("SetName",
                "Direction",
                "P-value",
                "Adjusted p-value (BH)",
                "Summarized t-statistic",
                "Number of genes annotated in the query set",
                "Number of query genes found in the dataset",
                "Per-gene raw expression scores",                
                colnames(pData(set)))
            
            return( res )
          }        
          )

setMethod(
          "gsealm_score",
          signature( query = "matrix", set="CMAPCollection"),
          function( query, set, predictor=NULL, ...) {
            if( is.null( predictor ) ){
              stop("To analyse query matrices, please provide the 'predictor'")
            }
            
            if( ncol( query ) < 2) {
              stop("Data matrix with multiple columns is required for sample label permutations.")
            }
            query <- ExpressionSet( query,
                                   phenoData=as( data.frame(predictor=factor( predictor ), 
                                     row.names=colnames( query )), 
                                     "AnnotatedDataFrame"))
            
            gsealm_score( query, set, predictor="predictor", ...)
          }
          )

setMethod(
          "gsealm_score",
          signature( query = "eSet", set="CMAPCollection"),
          function( query, set, element="exprs", ...) {
            query <- ExpressionSet(assayDataElement( query, element),
                                   phenoData=phenoData( query ),
                                   annotation=annotation( query )
                                   )
            
            gsealm_score( query, set, ...)
          }
          )

setMethod(
          "gsealm_score",
          signature( query = "eSet", set="GeneSetCollection"),
          function( query, set, element="exprs", ...) {
            query <- ExpressionSet( assayDataElement( query, element),
                                   phenoData=phenoData( query ),
                                   annotation=annotation( query )
                                   )
            gsealm_score( query, as(set, "CMAPCollection"), ...)
          }
          )

setMethod(
          "gsealm_score",
          signature( query = "eSet", set="GeneSet"),
          function( query, set, element="exprs", ...) {
            query <- ExpressionSet( assayDataElement( query, element),
                                   phenoData=phenoData( query ),
                                   annotation=annotation( query )
                                   )
            
            gsealm_score( query, as(set, "CMAPCollection"), ...)
          }
          )

setMethod(
          "gsealm_score",
          signature( query = "matrix", set="GeneSetCollection"),
          function( query, set, ...) {
            gsealm_score( query, as(set, "CMAPCollection"), ...)
          }
          )

setMethod(
          "gsealm_score",
          signature( query = "matrix", set = "GeneSet"),
          function( query, set,  ...) {
            gsealm_score( query, as(set, "CMAPCollection"), ...)
          }
          )

setMethod(
          "gsealm_score",
          signature( query = "ExpressionSet", set="GeneSet"),
          function( query, set, ...) {
            gsealm_score( query, as(set, "CMAPCollection"), ...)
          }
          )

setMethod(
          "gsealm_score",
          signature( query = "ExpressionSet", set="GeneSetCollection"),
          function( query, set,...){
            gsealm_score(query, as(set,"CMAPCollection"), ...)
          }
          )

## functions duplicated from the GSEAlm package, modified to work on sparce matrices
.gsNormalize <- function (dataset, incidence, fun1 = "/", 
                          fun2 = sqrt, removeShift = FALSE, removeStat = function( x) { mean(x, na.rm=TRUE)} , ...) 
  {
    dataset = as.matrix(dataset)
    if (removeShift) {
      colStats = apply(dataset, 2, removeStat)
      dataset = sweep(dataset, 2, STATS = colStats)
    }
    if (ncol(incidence) != nrow(dataset)) 
      stop(".gsNormalize: non-conforming matrices")
    outm = Matrix::crossprod(Matrix::t(incidence), dataset, ...)  
    rownames(outm) = rownames(incidence)
    colnames(outm) = colnames(dataset)
    normby = fun2(Matrix::rowSums(abs( incidence )))
    outm = sweep(outm, 1, normby, FUN = fun1)
  }


.gsealmPerm <- function (eSet, formula = "", mat, nperm, na.rm = TRUE, pooled = FALSE, 
                         detailed = FALSE, ...) {

  nSamp = ncol(eSet)
  if (formula == "") {
    nvar = 0
  } else {
    xvarnames = all.vars(formula)
    nvar <- length(xvarnames)
  }
  obsRaw = lmPerGene(eSet = eSet, formula = formula, na.rm = na.rm, 
    pooled = pooled)
  if (nvar > 0) {
    observedStats = .gsNormalize(obsRaw$tstat[2, ], incidence = mat, 
      ...)
  } else {
    observedStats = .gsNormalize(t(obsRaw$tstat), incidence = mat, 
      fun2 = identity, ...)
  }
  perm.eset = eSet
  i <- 1L
  if (nvar > 0) {
    permMat <- matrix(0, nrow = nrow(eSet), ncol = nperm)
    while (i < (nperm + 1)) {
      if (nvar >= 2) {
        splitter = pData(eSet)[, xvarnames[2]]
        if (nvar > 2) 
          splitter = as.list(pData(eSet)[, xvarnames[2:nvar]])
        label.perm = unsplit(lapply(split(1:nSamp, splitter), 
          sample), splitter)
        pData(perm.eset)[, xvarnames[1]] <- pData(eSet)[label.perm, 
                                                        xvarnames[1]]
      } else if (nvar == 1) {
        pData(perm.eset)[, xvarnames[1]] <- pData(eSet)[sample(1:nSamp), 
                                                        xvarnames[1]]
      }
      temp.results <- lmPerGene(eSet = perm.eset, formula = formula, 
                                na.rm = na.rm, pooled = pooled)
      permMat[, i] <- temp.results$tstat[2, ]
      i <- i + 1L
    }
    permMat <- .gsNormalize(permMat, incidence = mat, ...)
    rownames(permMat) = rownames(mat)
  }
  else if (nvar == 0) {
    permMat <- matrix(0, nrow = nrow(mat), ncol = nperm)
    rownames(permMat) = rownames(mat)
    for (i in 1:nperm) permMat[, i] = .gsNormalize(t(obsRaw$tstat), 
                                incidence = mat[, sample(1:ncol(mat))], fun2 = identity, 
                                ...)
  }
  ## convert to standard matrix objects
  observedStats <- as.matrix( observedStats )
  permMat <- as.matrix( permMat)
  
  if (!detailed) {
    return(.pvalFromPermMat(observedStats, permMat))
  }
  else return(
              list(
                   pvalues = .pvalFromPermMat(observedStats, 
                     permMat ), lmfit = obsRaw, stats = observedStats, perms = permMat)
              )
}


.pvalFromPermMat <- function (obs, perms) {
  N <- ncol(perms)
  pvals <- matrix(as.double(NA), nrow = nrow(perms), ncol = 2)
  dimnames(pvals) <- list(rownames(perms), c("Lower", "Upper"))
  tempObs <- rep(obs, ncol(perms))
  dim(tempObs) <- dim(perms)
  pvals[, 1] <- (1 + rowSums(perms <= tempObs))/(N + 1)
  pvals[, 2] <- (1 + rowSums(perms >= tempObs))/(N + 1)
  pvals
}
