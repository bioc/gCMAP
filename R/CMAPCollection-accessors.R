
.harmonizeDimnames <- function(object) {
  err <- function(conflicts)
    stop("assayData element dimnames conflict: ",
         paste(names(conflicts), collapse=", "))
  okNames <- list(featureNames(featureData(object)),
                  sampleNames(phenoData(object)))
  dimNames <- .assayDataDimnames(assayData(object))
  dimConflict <- function(dimNames, okNames, dim) {
    nm <- lapply(dimNames, "[[", dim)
    isConflict <- !sapply(nm, function(x, y) {
      is.null(x) || all.equal(x, y, check.attr=FALSE)
    }, okNames[[dim]])
    isNamed <- sapply(lapply(nm, names), length) > 0
    isNull <- sapply(nm, is.null)
    if (all(!isConflict & !isNamed & !isNull))
      return (FALSE)
    if (any(isConflict & !isNull))
      err(isConflict[!isNull])
    TRUE
  }
  if (dimConflict(dimNames, okNames, 1))
    featureNames(assayData(object)) <- okNames[[1]]
  if (dimConflict(dimNames, okNames, 2))
    sampleNames(assayData(object)) <- okNames[[2]]
  object
}

.assayDataDimnames <- function(assayData) {
  switch(storageMode(assayData),
         lockedEnvironment=,
         environment=eapply(assayData, dimnames),
         list=lapply(assayData, dimnames))
}

.annotatedDataFrameFromMatrix <- function(object, byrow=FALSE, ...) {
  ## contract: 'object' is matrix-like, with dim, rownames, colnames
  ## methods. Returns AnnotatedDataFrame with appropriate dimensions.
  dims <- dim(object)
  if (is.null(dims) || all(dims==0))
    annotatedDataFrameFrom(NULL, byrow=byrow, ...)
  else {
    n <- if (byrow) dims[1] else dims[2]
    nms <-
      if(byrow) rownames(object)
      else colnames(object)
    data <- data.frame(numeric(n), row.names=nms)[,FALSE]
    dimLabels <-
      if (byrow) c("featureNames", "featureColumns")
      else c("sampleNames", "sampleColumns")
    new("AnnotatedDataFrame", data=data, dimLabels=dimLabels)
  }
}

setMethod("initialize", "CMAPCollection",
          function(.Object,
                   assayData,
                   phenoData,
                   featureData,
                   members = new("dgCMatrix"),
                   signed,
                   ... ) {
            if (missing(assayData)) {
              if (missing(phenoData))
                phenoData <- annotatedDataFrameFrom(members, byrow=FALSE)
              if (missing(featureData))
                featureData <- annotatedDataFrameFrom(members, byrow=TRUE)
              .Object <- callNextMethod(.Object,
                                        phenoData = phenoData,
                                        featureData = featureData,
                                        members = members,
                                        ...)
            } else if (missing(members)) {
              if (missing(phenoData))
                phenoData <- annotatedDataFrameFrom(assayData, byrow=FALSE)
              if (missing(featureData))
                featureData <- annotatedDataFrameFrom(assayData, byrow=TRUE)
              .Object <- callNextMethod(.Object,
                                        assayData = assayData,
                                        phenoData = phenoData,
                                        featureData = featureData,
                                        ...)
            } else stop("provide at most one of 'assayData' or 'members' to initialize CMAPCollection",
                        call.=FALSE)

            if( missing( signed ) ) {
              pData(.Object)$signed <- rep(NA, ncol(.Object))
            } else {
              stopifnot( ncol( .Object ) == length( signed ) )
              stopifnot( is.logical( signed) )
              pData(.Object)$signed <- signed
            }
            .harmonizeDimnames(.Object)
          })


setMethod("CMAPCollection",
          signature=signature(assayData="missing"),
          function(assayData,
                   phenoData=AnnotatedDataFrame(),
                   featureData=AnnotatedDataFrame(),
                   protocolData=AnnotatedDataFrame(),
                   ...)
          {
            new("CMAPCollection",
                assayData=assayDataNew(members=new("dgCMatrix"),
                  phenoData=phenoData,
                  featureData=featureData,
                  annotation=annotation,
                  protocolData=protocolData, ...)
                )
          })

setMethod("CMAPCollection",
          signature=signature(assayData="environment"),
          function(assayData,
                   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                   annotation=character(),
                   protocolData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   ...)
          {
            new("CMAPCollection", assayData=assayData, phenoData=phenoData,
                featureData=featureData,
                annotation=annotation, protocolData=protocolData, ...)
          })


setMethod("CMAPCollection",
          signature=signature(assayData="Matrix"),
          function(assayData,
                   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                   annotation=character(),
                   protocolData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   ...)
          {
            assayData <- assayDataNew(members=assayData)
            new("CMAPCollection", assayData=assayData,
                phenoData=phenoData,
                featureData=featureData,
                annotation=annotation, protocolData=protocolData, ...)
          })

setMethod("CMAPCollection",
          signature=signature(assayData="matrix"),
          function(assayData,
                   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                   annotation=character(),
                   protocolData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   ...)
          {
            assayData <- assayDataNew(members=as(assayData, "dgCMatrix"))
            new("CMAPCollection", assayData=assayData,
                phenoData=phenoData,
                featureData=featureData,
                annotation=annotation, protocolData=protocolData, ...)
          })


setMethod("annotatedDataFrameFrom",
          signature(object="Matrix"),
          .annotatedDataFrameFromMatrix)

setAs("CMAPCollection", "matrix",
      function (from) as.matrix( members( from ) ) )

setAs("list", "CMAPCollection",
      function (from) {
        cmapData <- incidence( from )
        cmapData <- Matrix::t( cmapData )
        CMAPCollection( cmapData )
      })

setAs("GeneSetCollection", "CMAPCollection",
      function (from) {
        cmapData <- incidence(from)
        cmapData <- Matrix::t( cmapData )
        cmap <- CMAPCollection(
                               cmapData,
                               signed=ifelse( lapply(from, class) == "SignedGeneSet", TRUE, FALSE)
                               )
        
        from.anno <- unique( lapply( from, geneIdType))
        if( length (from.anno ) > 1) {
          annotation(cmap) <- "mixed"
        } else {
          annotation(cmap) <- annotation(from.anno[[1]])
        }
        desc <- sapply( from, description)
        if( ! all( desc == "")) {
          pData(cmap)$description <- desc
        }
        cmap
      })

setAs("GeneSet", "CMAPCollection",
      function (from) {
        if( is.na( setName( from ) ) ) {
          setName(from) <- "1"
        }
        from <- GeneSetCollection(from)
        cmapData <- incidence( from )
        cmapData <- Matrix::t( cmapData )
        cmap <- CMAPCollection(
                                  cmapData,
                                  signed=ifelse( lapply(from, class) == "SignedGeneSet", TRUE, FALSE)
                                  )
        from.anno <- unique( lapply( from, geneIdType))
        if( length (from.anno ) > 1) {
          annotation(cmap) <- "mixed"
        } else {
          annotation(cmap) <- annotation(from.anno[[1]])
        }
        cmap
      })

setAs("CMAPCollection", "GeneSetCollection",
      function (from) {
        
        ## try to identify organism identifier
        organism <- tryCatch({
          pkg <- annotation(from)
          if (length(pkg) == 1 && nchar(pkg) > 0) 
            getAnnMap("ORGANISM", pkg)
          else ""
        }, error = function(err) "")
        if ( length (pkg) == 0 ) annotation(from) <- "None"
        
        ## create individual SignedGeneSets
        set.list <- lapply( sampleNames (from), function( n ) {
          dat <- members(from)[,n]
          ids <- featureNames(from)[dat != 0]
          geneSign <- ifelse( dat[dat != 0 ] == 1, "up", "down")
          SignedGeneSet(ids,
                        geneSign = geneSign, 
                        setName=n,
                        geneIdType = AnnoOrEntrezIdentifier(annotation(from)),
                        shortDescription = experimentData(from)@title, 
                        longDescription = abstract(from), organism = organism, 
                        pubMedIds = pubMedIds(experimentData(from)), urls = experimentData(from)@url, 
                        contributor = experimentData(from)@name 
                        )
        })
        ## convert SignedGeneSets to GeneSets based on 'signed' CMAPCollection entries
        set.list <- lapply( seq( ncol(from) ), function(n) {
          if( signed(from)[n] == FALSE ) {
            as(set.list[[n]], "GeneSet")
          } else {
            set.list[[n]]
          }})
        
        ## return GeneSetCollection
        GeneSetCollection(set.list)
      })

setAs("CMAPCollection", "GeneSet",
      function (from) {
        if( ncol(from) > 1) {
          stop( "Cannot coerce a CMAPCollection with multiple sets into a single GeneSet.\nConsider a GeneSetCollection instead.")
        }
        ## try to identify organism identifier
        organism <- tryCatch({
          pkg <- annotation(from)
          if (length(pkg) == 1 && nchar(pkg) > 0) 
            getAnnMap("ORGANISM", pkg)
          else ""
        }, error = function(err) "")
        if ( length (pkg) == 0 ) annotation(from) <- "None"
        
        ## create GeneSets
        dat <- members(from)
        ids <- featureNames(from)[dat[,1] != 0]
        if( signed( from ) == FALSE) {
          GeneSet(ids,
                  setName=sampleNames(from),
                  geneIdType = AnnoOrEntrezIdentifier(annotation(from)),
                  shortDescription = experimentData(from)@title, 
                  longDescription = abstract(from), organism = organism, 
                  pubMedIds = pubMedIds(experimentData(from)), urls = experimentData(from)@url, 
                  contributor = experimentData(from)@name 
                  )
        } else {
          geneSign <- ifelse( dat[dat[,1] != 0, 1] == 1, "up", "down")
          SignedGeneSet(ids,
                        geneSign = geneSign, 
                        setName=sampleNames(from),
                        geneIdType = AnnoOrEntrezIdentifier(annotation(from)),
                        shortDescription = experimentData(from)@title, 
                        longDescription = abstract(from), organism = organism, 
                        pubMedIds = pubMedIds(experimentData(from)), urls = experimentData(from)@url, 
                        contributor = experimentData(from)@name 
                        )
        }
      })

setMethod(
          "induceCMAPCollection",
          signature( "eSet" ),
          function( eset, element, lower=NULL, higher=NULL, sign.sets=TRUE) {
            if( ! is.null(lower) && ! is.null(higher) && higher == lower) {
              stop("Please specify two different cutoffs as 'higher' and 'lower' parameters.")
            }
            
            if(! element %in% assayDataElementNames(eset) ) stop(paste( "AssayDataElement", element, "not found."))
            eset <- memorize( eset, element)
            ade <- assayDataElement( eset, element )

            gss <- lapply(
                          seq(1:ncol(eset)),
                          function( n ) {
                            if (! is.null( lower )) {
                              down <- which( ade[,n] < lower )
                            } else {
                              down <- NULL
                            }                            
                            if (! is.null( higher )) {
                              up   <- which( ade[,n] > higher )
                            } else {
                              up <- NULL
                            }                            
                            list( j = c(down, up),
                                 x = c(rep(-1, length(down)), rep(1, length(up)))
                                 )
                          })
            
            i <- unlist(
                        sapply( seq( length( gss ) ), function( m ) {
                          rep( m, length( gss[[ m ]]$j ) )
                        }))
            j <- unlist(sapply(gss ,function( m ) {m$j }))
            x <- unlist(sapply(gss ,function( m ) {m$x }))            
            if( sign.sets == TRUE ){
              set.signs <- rep(TRUE, ncol(eset))
            } else {
              set.signs <- rep(FALSE, ncol(eset))
            }
            cmap <- CMAPCollection(
                                      Matrix::t( sparseMatrix(i=i,
                                                              j=j,
                                                              x=x,
                                                              dims=list(ncol(eset), nrow(eset)),
                                                              dimnames = list(sampleNames(eset), featureNames(eset)))
                                                ),
                                      phenoData = as(pData(eset), "AnnotatedDataFrame"),
                                      featureData = as(fData(eset),"AnnotatedDataFrame"),
                                      signed=set.signs)
            notes( cmap ) <- list(CMAPCollection=paste("induced from channel",element,"selecting scores <",lower,"or >",higher))
            cmap
          }
          )

setMethod(
          "induceCMAPCollection",
          signature( "matrix" ),
          function( eset, element, ...) {
            induceCMAPCollection( ExpressionSet( eset ), element="exprs", ... )
          })

setMethod(
          "geneIds",
          signature( "CMAPCollection" ),
          function( object, ... ) {
            dat <- members( object )
            gene.ids <- lapply( seq( ncol( dat ) ), function( n ) {
              featureNames( object )[ which(dat[,n] != 0 ) ]
            })
            names(gene.ids) <- sampleNames( object )
            if( length (gene.ids) == 1) {
              return( gene.ids[[1]] )
            } else {
              return( gene.ids )
            }
          }
          )

setMethod(
          "members",
          signature( "CMAPCollection" ),
          function( object) {
            assayData(object)[["members"]]
          }
          )

setMethod(
          "members",
          signature( "CMAPCollection" ),
          function( object) {
            assayData(object)[["members"]]
          }
          )

setMethod(
          "signed",
          signature( "CMAPCollection" ),
          function( object) {
            signs <- pData(object)$`signed`
            names( signs ) <- sampleNames(object)
            signs
          }
          )


setReplaceMethod("signed", "CMAPCollection",
                 function(x, value) {
                   pData(x)$`signed` <- value
                   validObject( x )
                   x
                 }
                 )

setMethod(
          "mergeCollections",
          signature=signature( x="CMAPCollection", y="CMAPCollection" ),
          function( x, y ) {
            common.genes <- intersect( featureNames( x ), featureNames( y ) )
            c.members <- Matrix::cBind( members( x[common.genes,] ), members( y[common.genes,] ) )

            dupl.sets <- duplicated( colnames( c.members ) )
            if( any( dupl.sets) ){
              message( "Duplicated set / column names detected: attached suffix.")
            }
            
            n <- 2
            while( any( dupl.sets ) ) {
              colnames( c.members )[dupl.sets] <- paste( colnames( c.members )[dupl.sets], n, sep=".")
              dupl.sets <- duplicated( colnames( c.members ) )
              n <- n + 1
            }

            ## try to merge phenoData
            if( ncol( pData( x ) ) == ncol( pData( y ) ) && colnames( pData( x ) ) == colnames( pData( y ) ) ) {
              p.data <- rbind( pData(x), pData(y) )
              row.names(p.data) <- colnames( c.members )
            } else {
              message("Could not merge phenoData frame: colnames don't match.")
              p.data <- data.frame(row.names = colnames( c.members ) )
            }

            p.data <- as(p.data, "AnnotatedDataFrame")

            CMAPCollection(c.members, phenoData=p.data)
          }
          )

