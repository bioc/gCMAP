##---------- queries are expression scores

setMethod(
  "gsealm_jg_score",
  signature( query = "matrix_or_big.matrix", sets="CMAPCollection"),
  function( query, sets, removeShift=FALSE, respect.sign=TRUE,keep.scores=FALSE, ...) {
   
    if( respect.sign == FALSE) {
      signed(sets) <- rep(FALSE, ncol(sets))
    }

    score <- .jg( query, sets, removeShift=removeShift, respect.sign=respect.sign)
    
    ## replace NaN values produced for empty sets
    score[["t"]][ score[["setSize"]] == 0 ] <- 0
    score[["pval"]][ score[["setSize"]] == 0 ] <- 1
    
    ## store per-gene scores as data-column:gene-set list-of-list
    if( keep.scores == TRUE) {
      gene.scores <- featureScores( query, sets)
    } else {
      gene.scores <- NA
    }
    
    ## create separate CMAPResults object for each query column
    results <- mclapply( seq( ncol( query )), function(n ) {
      ## return different 'trend' descriptors for signed and unsigned sets
      trend <- sapply( sampleNames( sets ), function( set.name) {
        if( signed(sets)[set.name] == FALSE) {
          trend <- ifelse(score[["t"]][set.name, n] >=0, "up", "down")
        } else {
          trend <- ifelse(score[["t"]][set.name,n] >=0, "correlated", "anticorrelated")
        }
      })

      if(! all(is.na( gene.scores ))) { 
        geneScores <- I(gene.scores[[n]])
      } else {
        geneScores <- NA
      }
      
      res <- CMAPResults(
        data=data.frame(
          set = sampleNames(sets),
          trend = trend,
          pval = score[["pval"]][,n],
          padj =  p.adjust( score[["pval"]][,n], method="BH"),
          effect = score[["t"]][,n],
          nSet = Matrix::colSums( abs( members (sets) ) ),
          nFound = score[["setSize"]][,n],
          geneScores = geneScores,
          pData(sets)),
        docs="\n Parametric 'JG' score summary.\n P-values were adjusted with the 'p-adjust' function using method 'BH'."
      )
      
      varMetadata(res)$labelDescription <- 
        c("SetName",
          "Direction",
          "p-value from standard normal distribution",
          "Adjusted p-value (BH)",
          "Summarized t-statistic",
          "Number of genes annotated in the query set",
          "Number of annotated genes found in the data",
          "Per-gene raw scores",
          colnames(pData(sets)))
      
      res
    })
    
    ## return a single CMAPResult object for one-column queries
    names(results) <- colnames( query )
    if( length( results) == 1 ) return( results[[1]] )
    return( results)
  })

setMethod(
  "gsealm_jg_score",
  signature( query = "numeric", sets="CMAPCollection"),
  function( query, sets, ...) {
    if( is.null( names( query ) ) ){
      stop( "The 'query' vector does not provide gene names.")
    }
    query.matrix <- matrix(query, ncol=1, dimnames=list(names(query), "Query"))
    gsealm_jg_score(query.matrix, sets, ...)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "eSet", sets="CMAPCollection"),
  function( query, sets, element="z", ...){
    gsealm_jg_score( as( assayDataElement(query, element), "matrix"), sets, ...)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "matrix", sets="GeneSetCollection"),
  function( query, sets, ...) {
    if( is.null( row.names( query ) ) ) {
      stop( "The query matrix does not contain row.names.")
    }
    sets <- as(sets, "CMAPCollection")
    gsealm_jg_score( query, sets, ...)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "matrix", sets = "GeneSet"),
  function( query, sets, ...){
    sets <- as(sets, "CMAPCollection")
    gsealm_jg_score( query, sets, ...)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "numeric", sets="GeneSet"),
  function( query, sets, ...) {
    if( is.null( names( query ) ) ){
      stop( "The 'query' vector does not provide gene names.")
    }
    query.matrix <- matrix(query, ncol=1)
    row.names(query.matrix) <- names(query)
    colnames(query.matrix) <- "Query"
    gsealm_jg_score(query.matrix, sets,...)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "numeric", sets="GeneSetCollection"),
  function( query, sets,...) {
    if( is.null( names( query ) ) ){
      stop( "The 'query' vector does not provide gene names.")
    }
    query.matrix <- matrix(query, ncol=1)
    row.names(query.matrix) <- names(query)
    colnames(query.matrix) <- "Query"
    gsealm_jg_score(query.matrix, sets, ...)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "eSet", sets="GeneSet"),
  function( query, sets, element="z", ...) {
    gsealm_jg_score( as( assayDataElement(query, element), "matrix"), sets, ...)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "eSet", sets="GeneSetCollection"),
  function( query, sets,element="z", ...){
    gsealm_jg_score(as( assayDataElement(query, element), "matrix"), sets,...)
  }
)

##-------- queries are gene set(s)

setMethod(
  "gsealm_jg_score",
  signature( query = "CMAPCollection", sets="eSet"),
  function( query, sets, element="z", removeShift=FALSE, respect.sign=TRUE, keep.scores=FALSE,...) {
    
    if( respect.sign == FALSE) {
      signed(query) <- rep(FALSE, ncol(query))
    }
    
    set.matrix <- assayDataElement(sets, element)

    score <- .jg( set.matrix, query, removeShift=removeShift, ## rows=query, columns=set
                  respect.sign=respect.sign)

    ## store per-gene scores as gene-set:data-column list-of-list
    if( keep.scores == TRUE) {
      gene.scores <- featureScores( query, sets, element=element)
    } else{
      gene.scores <- NA
    }
    
    ## create separate CMAPResults object for each sample in the eSet
    results <- mclapply( seq( ncol( query )), function( n) {
      
      ## return different 'trend' descriptors for signed and unsigned sets
      if( signed(query)[n] == FALSE ) {
        trend <- ifelse(score[["t"]][n,] >=0, "up", "down")
      } else {
        trend <- ifelse(score[["t"]][n,] >=0, "correlated", "anticorrelated")
      }

      if( all(is.na( gene.scores ))) { 
        geneScores <- lapply( seq(ncol( gene.scores[[n]]) ), function( m ) {
                          g <- gene.scores[[n]][,m]
                          attr( g, "sign") <- attr( gene.scores[[n]], "sign")
                          g
                        })
        names( geneScores ) <- sampleNames(sets)
        geneScores <- I( geneScores )
      } else {
        geneScores <-  NA
      }
                        
      res <- CMAPResults( 
        data = data.frame(
          set = sampleNames(sets),
          trend = trend,
          pval = score[["pval"]][n,],
          padj =  p.adjust( score[["pval"]][n,], method="BH"),
          effect = score[["t"]][n,],
          nSet = rep( Matrix::colSums( abs( members (query) ) )[n], ncol(sets)),
          nFound = score[["setSize"]][n,],
          geneScores = geneScores, 
          pData(sets)),
        docs="\n Parametric 'JG' score summary.\n P-values were adjusted with the 'p-adjust' function using method 'BH'."
      )
      
      varMetadata(res)$labelDescription <- 
        c("SetName",
          "Direction",
          "p-value from standard normal distribution",
          "Adjusted p-value (BH)",
          "Summarized t-statistic",
          "Number of genes annotated in the gene set",
          "Number of annotated genes found in the data",
          "Per-gene raw scores",
          colnames(pData(sets)))
      
      res
    })
    
    ## return a single CMAPResult object for one-column queries
    names(results) <- sampleNames( query )
    if( length( results) == 1 ) return( results[[1]] )
    return( results)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "CMAPCollection", sets="numeric"),
  function( query, sets, ...) {
    set.matrix <- matrix(sets, ncol=1)
    row.names(set.matrix) <- names(sets)
    colnames(set.matrix) <- "Dataset"
    sets <- ExpressionSet(set.matrix)
    gsealm_jg_score(query, sets, element="exprs", ...)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "CMAPCollection", sets="matrix"),
  function( query, sets, ...){
    if( is.null( row.names( sets ) ) ){
      stop( "The 'sets' matrix does not contain row.names.")
    }
    sets <- ExpressionSet(sets)
    gsealm_jg_score( query, sets, element="exprs",...)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "GeneSetCollection", sets="matrix"),
  function( query, sets, ...) {
    if( is.null( row.names( sets ) ) ){
      stop( "The 'sets' matrix does not contain row.names.")
    }
    query <- as(query, "CMAPCollection")
    sets <- ExpressionSet(sets)
    gsealm_jg_score( query, sets, element="exprs",...)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "GeneSet", sets = "matrix"),
  function( query, sets, ...){
    if( is.null( row.names( sets ) ) ){
      stop( "The 'sets' matrix does not contain row.names.")
    }
    query <- as(query, "CMAPCollection")
    sets <- ExpressionSet(sets)
    gsealm_jg_score( query, sets, element="exprs", ...)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "GeneSet", sets="numeric"),
  function( query, sets, ...) {
    set.matrix <- matrix(sets, ncol=1)
    row.names(set.matrix) <- names(sets)
    colnames(set.matrix) <- "Dataset"
    gsealm_jg_score(as(query, "CMAPCollection"), set.matrix, ...)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "GeneSetCollection", sets="numeric"),
  function( query, sets,...) {
    if( is.null( names( sets ) ) ){
      stop( "The 'sets' vector does not provide gene names.")
    }
    set.matrix <- matrix(sets, ncol=1)
    row.names(set.matrix) <- names(sets)
    colnames(set.matrix) <- "Dataset"
    gsealm_jg_score(as(query, "CMAPCollection"), set.matrix, ...)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "GeneSet", sets="eSet"),
  function( query, sets, element="z", ...) {
    gsealm_jg_score( as(query, "CMAPCollection"), sets, ...)
  }
)

setMethod(
  "gsealm_jg_score",
  signature( query = "GeneSetCollection", sets="eSet"),
  function( query, sets,element="z", ...){
    gsealm_jg_score(as(query, "CMAPCollection"), sets,...)
  }
)

###

.jg <- function(query, sets, removeShift=FALSE, respect.sign=TRUE) {
  
  ## intersect query gene ids with gene ids in the collection
  a.matrix <- Matrix::t( members (sets) )
  common.genes <- colnames(a.matrix)[na.omit(match(row.names(query), colnames(a.matrix)))]
  
  a.matrix <- try(
    a.matrix[,common.genes]
  )
  
  if (class(a.matrix) == "numeric") a.matrix <- matrix(a.matrix, nrow=1) ## single gene.sets
  if (class(a.matrix) == "try-error" | sum( abs ( a.matrix) ) == 0) {
    stop( "None of the user-provided gene ids were found in the gene set collection.", call. = FALSE)
  }
  
  if( respect.sign == FALSE) {
    a.matrix <- abs( a.matrix )
    signed(sets) <- rep(FALSE, ncol(sets))
  }

  ## preprocess query matrix / BigMatrix
  ## extract big.matrix component
  if( inherits( query, "BigMatrix")){
    query <- query$bigmat
  }
  
  if( inherits( query, "matrix")){
    ## set missing values to 0
    if( removeShift == TRUE){
      query[ is.na( query ) ] <- 0
      query <- scale(query, center=TRUE, scale=FALSE)
    }

  } else if( inherits( query, "big.matrix")){
    if( removeShift == TRUE){
      stop( "The 'removeShift' option is currently not supported for eSets containing 'BigMatrix' or 'big.matrix' assayData objects.")
    }
  }
  
  rttStats <- matrix( query[ common.genes,], nrow=length( common.genes),
                      dimnames=list(common.genes, colnames(query)))
  rttStats[ is.na( rttStats ) ] <- 0

  ## compute aggregated statistic for each GeneSet
  tA <- matrix( a.matrix %*% rttStats,
                ncol=ncol(query),
                dimnames=list(sampleNames(sets), colnames(query))
  )
  
  setSize <- Matrix::rowSums(abs(a.matrix))
  tAadj <- apply(tA, 2, function(x) x/sqrt(setSize)) ## adjust for GeneSet size
  tAadj <- matrix( tAadj, ncol=ncol(query), dimnames=list(sampleNames(sets), colnames(query))) 
  pval <- 2*pnorm(abs(tAadj), lower.tail = FALSE) ## get p-value from Normal Distribution
  
  setSize <- matrix( rep(setSize, ncol(query)),
                     ncol=ncol(query),
                     dimnames=list(sampleNames(sets), colnames(query))
  )
  
  ## return all results
  list(t= tAadj, pval=pval, setSize=setSize)
}
