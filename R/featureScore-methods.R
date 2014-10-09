
setMethod(
          "featureScores",
          signature( query = "CMAPCollection", dat="matrix"),
          function( query, dat, simplify=TRUE){
            common.genes <- intersect(featureNames(query), row.names(dat))
            if( length( common.genes ) == 0) {
              stop("None of CMAPCollection's featureNames found in data matrix row.names.")
            } else {
              query <- query[common.genes,]
              dat <- dat[common.genes,,drop=FALSE]
            }
            
            ## generate gene-set:data-column list-of-list
            res <- lapply( seq( ncol( query ) ), function( n ) {
              set.members <- which( members(query)[,n] != 0)
              sign <-  members(query)[set.members,n]
              
              query.score <- lapply( seq( ncol( dat )), function( m ){
                score <- dat[set.members,m,drop=TRUE]
                attr(score, "sign") <- ifelse( sign == 1, "up", "down")
                return( score )
              })
              names(query.score) <- colnames( dat )
              query.score
            })
            names( res ) <- sampleNames( query )
            
            ## collect scores for each set in a single matrix
            if( simplify == TRUE ){
              res <- lapply( res, function( x ) {
                if( length( x ) > 1){ ## multiple score vectors ?
                  s <- do.call(cbind, x)
                  attr(s, "sign") <- attr(x[[1]], "sign")
                  s
                } else {
                  x ## single score vector
                }
              })  
            }
            ## for list of length 1, return first element
            return( res )
          })

setMethod(
          "featureScores",
          signature( query = "matrix", dat="CMAPCollection"),
          function( query, dat ){    

            common.genes <- intersect(row.names(query), featureNames(dat))
            if( length( common.genes ) == 0) {
              stop("None of CMAPCollection's featureNames found in data matrix row.names.")
            } else {
              query <- query[common.genes,,drop=FALSE]
              dat <- dat[common.genes,]
            }    
            
            ## generate data-column:gene-set list-of-list
            res <- lapply( seq( ncol( query ) ), function( n ) {
              query.scores <- lapply( seq( ncol( dat) ), function( m ) {
                set.members <- which( members(dat)[,m] != 0)
                sign <-  members(dat)[set.members,m]
                
                ## prevent gene sets with a single member from dropping gene id
                score <- query[set.members,n,drop=FALSE]
                gene.id <- row.names( score )
                score <- as.numeric( score )
                names(score) <- gene.id
                
                attr(score, "sign") <- ifelse( sign == 1, "up", "down")
                return( score )
              })
              names(query.scores) <- sampleNames( dat )
              query.scores
            })
            
            names( res ) <- colnames( query )
            return( res )
          })

setMethod(
          "featureScores",
          signature( query = "CMAPCollection", dat="BigMatrix"),
          function( query, dat, simplify=TRUE){
            common.genes <- intersect(row.names(query), row.names(dat))
            if( length( common.genes ) == 0) {
              stop("None of CMAPCollection's featureNames found in row.names of BigMatrix dat.")
            } else {
              query <- query[common.genes,,drop=FALSE]
              dat <- dat[common.genes,,drop=FALSE]
            }      
            featureScores( query, dat, simplify=simplify )
          })

setMethod(
          "featureScores",
          signature( query = "BigMatrix", dat="CMAPCollection"),
          function( query, dat ){
            common.genes <- intersect(row.names(query), row.names(dat))
            if( length( common.genes ) == 0) {
              stop("None of the query's row.names (BigMatrix) overlap with row.names of dat (CMAPCollection).")
            } else {
              query <- query[common.genes,,drop=FALSE]
              dat <- dat[common.genes,,drop=FALSE]
            }      
            featureScores( query,  dat )
          })

setMethod(
          "featureScores",
          signature( query = "CMAPCollection", dat="eSet"),
          function( query, dat, element="z",simplify=TRUE ){
            if( ! element %in% assayDataElementNames(dat)) {
              stop( sprintf( "Element %s not found in assayDataElementNames.", element))
            }
            dat <- assayDataElement(dat, element)
            featureScores( query, dat, simplify=simplify )
          })

setMethod(
          "featureScores",
          signature( query = "eSet", dat="CMAPCollection"),
          function( query, dat, element="z" ){
            if( ! element %in% assayDataElementNames(query)) {
              stop( sprintf( "Element %s not found in assayDataElementNames.", element))
            }
            query <- assayDataElement(query, element)
            featureScores( query, dat )
          })

setMethod(
          "featureScores",
          signature( query = "CMAPCollection", dat="numeric"),
          function( query, dat ){
            featureScores( query, as.matrix( dat ), simplify=FALSE )
          })

setMethod(
          "featureScores",
          signature( query = "numeric", dat="CMAPCollection"),
          function( query, dat ){
            featureScores( as.matrix( query ), dat )
          })

##-------- methods returning only character vectors

setMethod(
          "featureScores",
          signature( query="CMAPCollection", dat="CMAPCollection"),
          function( query, dat){
            res <- lapply( sampleNames( query ), function( set.name ){
              gene.ids <- geneIds( query[, set.name])
              common.genes <- intersect(gene.ids, featureNames(dat))
              ids <- geneIds( dat[common.genes,] )
              sapply( ids, function( x){
                scores <- rep(NA, length(x))
                names(scores) <- x
                scores
              })
            })
            return( res )
          })

