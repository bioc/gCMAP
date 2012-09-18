setMethod("initialize", signature(.Object="CMAPResults"),
          function(.Object, data = data.frame(), varMetadata = data.frame(),
                   docs="", errors=list(), ...) {
            tryCatch({
              if (missing(varMetadata)) {
                if (!missing(data)) stopifnot( class(data) == "data.frame" )
                varMetadata <- data.frame(labelDescription = rep(NA, ncol(data)))
                row.names(varMetadata) <- as.character(colnames(data))
              } else {
                stopifnot( class(varMetadata) == "data.frame")
                if (!"labelDescription" %in% colnames(varMetadata))
                  varMetadata[["labelDescription"]] <- rep(NA, nrow(varMetadata))
                row.names(varMetadata) <- names(data)
              }
              varMetadata[["labelDescription"]] <- as.character(varMetadata[["labelDescription"]])
            }, error=function(err) {
              stop(conditionMessage(err),
                   "\n  AnnotatedDataFrame 'initialize' could not update varMetadata:",
                   "\n  perhaps pData and varMetadata are inconsistent?")
            })

            ## sort results
            if( "pval" %in% colnames( data ) ) {
              data <- data[order( data[,"pval"], decreasing =  FALSE ),]
            } else if( "effect" %in% colnames( data ) ) {
              data <- data[order( abs( data[,"effect"] ), decreasing =  TRUE ),]
            }
            
            callNextMethod(.Object, data=data, varMetadata=varMetadata, docs=docs, errors=errors, ...)
          })

setMethod("CMAPResults",
          signature(data="missing", varMetadata="missing"),
          
          function(data, varMetadata, dimLabels=c("rowNames", "columnNames"), docs="", errors=list(), ...)
          {
            CMAPResults(data=data.frame(), varMetadata=data.frame(),
                        dimLabels=dimLabels, docs=docs, errors=errors, ...)
          })

setMethod("CMAPResults",
          signature(data="data.frame", varMetadata="missing"),
          function(data, varMetadata,
                   dimLabels=c("rowNames", "columnNames"), docs="", errors=list(), ...)
          {
            varMetadata <- data.frame(labelDescription = rep(NA, ncol(data)))
            row.names(varMetadata) <- names(data)
            CMAPResults(data=data, varMetadata=varMetadata,
                        dimLabels=dimLabels, docs=docs, errors=errors,...)
          })

setMethod("CMAPResults",
          signature(data="data.frame", varMetadata="data.frame"),
          function(data, varMetadata,
                   dimLabels=c("rowNames", "columnNames"), docs="", errors=list(),...)
          {
            if (!"labelDescription" %in% colnames(varMetadata))
              varMetadata[["labelDescription"]] <- rep(NA, nrow(varMetadata))
            row.names(varMetadata) <- names(data)
            varMetadata[["labelDescription"]] <-
              as.character(varMetadata[["labelDescription"]])
            new("CMAPResults", data=data, varMetadata=varMetadata,
                dimLabels=dimLabels, docs=docs, errors=errors, ...)
          })


setMethod("show", "CMAPResults",
          function(object) {
            if( all( dim( object) == 0 ) ) {
              cat("Empty",
                  class( object ),
                  "object\n")
            } else {
              cat("\n", class( object ), "object",
                  "with the following data slots:\n",
                  paste( varLabels(object), collapse=", "))
              cat( "\n for",
                   nrow( object ),
                   "gene sets.")
              if( "padj" %in% varLabels( object ) ) {
                cat("\n",length( which( object$padj < 0.05)),
                    "test(s) obtained an adjusted p-value < 0.05\n")
              }
              cat( docs( object ),"\n\n" )              
              dat <- cmapTable(object, n=5)
              
              print(dat)
              cat( "... (only top 5 results shown, use 'cmapTable' function to see all) ... \n")
            }
          })


setMethod("cmapTable", "CMAPResults",
          function(object, columns=NULL, n=NULL) {
            
            if( ! all( dim( object) == 0 ) ) {
              res <- pData(object)[, which( colnames( pData(object )) != "geneScores" ), drop=FALSE]
              
              ## per-gene score
              gene.scores <- geneScores( object )
              if( ! is.null(gene.scores)) {
                if(! all(is.na( gene.scores ))) {
                  geneScores <- unlist( lapply( gene.scores, function(x) {
                    if( inherits( x, "matrix" )){
                      sprintf("%d x %d matrix", nrow(x), ncol(x))
                    } else {                   
                      sprintf("%d scores", length(x))
                    }
                  }
                  )
                  )
                  res <- cbind(res, geneScores)
                }
              }
                            
              ## select user-defined subset of columns and rows
              if (! is.null( columns ) ) res <- subset(res, select=columns)
              if (! is.null( n) ) res <- head( res, n)
              row.names(res) <- NULL
              return( res )
              
            } else {
              cat("CMAPRresult object is empty.\n")
            }            
          })

## accessor and replacement functions


setMethod("set", "CMAPResults",
          function(object) {
            object$set
          }
)

setReplaceMethod("set", "CMAPResults", function(x, value) {
  x$set <- value
  x <- update(x)
  x
})


setMethod("trend", "CMAPResults",
          function(object) {
            object$trend
          }
)

setReplaceMethod("trend", "CMAPResults", function(x, value) {
  x$trend <- value
  x <- updateObject(x)
  x
})

setMethod("pval", "CMAPResults",
          function(object) {
            x <- object$pval
            names(x) <- set(object)
            x
          }
)

setReplaceMethod("pval", "CMAPResults", function(x, value) {
  x$pval <- as.numeric(value)
  x
})

setMethod("padj", "CMAPResults",
          function(object) {
            x <- object$padj
            names(x) <- set(object)
            x
          }
)


setReplaceMethod("padj", "CMAPResults", function(x, value) {
  x$padj <- as.numeric(value)
  x <- update( x )
  x
})

setMethod("zscores", "CMAPResults",
          function(object) {
            
            if( "padj" %in% varLabels( object ) ) {
              x <- zScores(object$padj, object$effect)
            } else { stop("CMAPResults object does not contain adjusted p-value (padj) data.")
            }
            names(x) <- set(object)
            x
          }
)


setMethod("effect", "CMAPResults",
          function(object) {
            x <- object$effect
            names(x) <- set(object)
            x
          }
)

setReplaceMethod("effect", "CMAPResults", function(x, value) {
  x$effect <- as.numeric(value)
  x <- update(x)
  x
})

setMethod("nSet", "CMAPResults",
          function(object) {
            object$nSet
          }
)

setReplaceMethod("nSet", "CMAPResults", function(x, value) {
  x$nSet <- as.integer(value)
  x <- update(x)
  x
})

setMethod("nFound", "CMAPResults",
          function(object) {
            object$nFound
          }
)

setReplaceMethod("nFound", "CMAPResults", function(x, value) {
  x$nFound <- as.integer(value)
  x <- update(x)
  x
})

setMethod("docs", "CMAPResults",
          function(object) {
            object@docs
          }
)


setReplaceMethod("docs", "CMAPResults", function(x, value) {
  x@docs <- value
  x
})

setMethod("errors", "CMAPResults",
          function(object) {
            object@errors
          }
)

setReplaceMethod("errors", "CMAPResults", function(x, value) {
  x@errors <- value
  x
})

setMethod("labels", "CMAPResults",
          function(object) {
            varMetadata(object)
          }
)

setReplaceMethod("labels", "CMAPResults", function(x, value) {
  varMetadata(x) <- value ## must be a data.frame with column 'labelDescription'
  x
})

setMethod("geneScores", "CMAPResults",
          function(object) {
            object$geneScores
          }
)

setAs("CMAPResults", "data.frame",
      function (from) pData(from)
)


zScores <- function(pval, direction=NULL, limit=.Machine$double.xmin) {
  if( !is.null( limit ) ){
    pval[which(pval < limit )] <- limit ## set lower limit to avoid Inf/-Inf zscores
  }
  z <- qnorm( pval/2, lower.tail=FALSE )
  if ( !is.null( direction) ) {
    z <-  z * sign( direction )
  }
  z
}

## plots

setMethod("plot", "CMAPResults",
          function(x, y=NULL,
                   strip.effect="effect",
                   strip.pval="padj",
                   strip.cutoffs=c(-2,2),
                   strip.bounds=c(-4,4),
                   strip.col=c("blue","white","red"),
                   strip.anno="set",
                   strip.subset=1:5,
                   strip.labels=NULL,
                   strip.layout=c(0.45, 0.1, 0.45),
                   set.inf=20,
                   ...
          ) {
            
            ## extract (sorted) data.frame
            dat <- cmapTable( x )
            
            ## retreive or calculate effect
            if( is.null( strip.effect) ) {
              effect <- zScores( dat[,strip.pval] ) ## convert to z-scores from normal distribution
              names( effect ) <- row.names( dat )
              effect.label <- "z-score"
              
            } else if( strip.effect %in% varLabels( x ) ) {
              effect <- dat[,strip.effect] ## retrieve effect sizes
              names(effect) <- row.names( dat )
              effect.label <- labels( x )[strip.effect,]
              
            } else if ( strip.pval %in% varLabels( x ) ){ ## retrieve pvalues
              message(sprintf("Column %s not found in CMAPResults object: Transforming %s to unsigned z-scores.", strip.effect, strip.pval))
              ## calculate zscores from pvalues
              effect <- zScores(dat[,strip.pval]) ## convert to z-scores from normal distribution
              names(effect) <- row.names( dat )
              effect.label <- "z-score"
              
            } else { ## neither effects nor pvalues are found: fail
              stop(sprintf("Neither column %s nor %s found in CMAPResults object.", strip.effect, strip.pval))
            }
            
            ## check if the requested annotation can be found
            if( ! is.null(strip.anno ) ) { ## annotate instances
              if( length( strip.anno ) > 1 | ! inherits(strip.anno, "character") ) {
                stop("Parameter 'strip.anno' must be a character vector of length 1.")
              }
              if( ! strip.anno %in% varLabels( x ) ) {
                stop("The 'strip.anno' column %s was not found in the CMAPResults object.", names(strip.anno))
              } else {
                groups <- dat[,strip.anno]
              }
            }
            
            ## check if the 'strip.subset' parameter is valid
            try( strip.labels <- levels( factor( groups[ strip.subset ] ) ),  silent =TRUE )
            if( class( strip.labels ) == "try-fail") {
              stop("Could not subset the data with the provided 'strip.subset' parameter.")
            }
            
            ## check if strip.layout makes sense
            if(! length( strip.layout != 3)) {
              stop("'strip.layout' must contain three elements.")
            }
            
            density.plot <- .score_density( effect=effect, 
                                            strip.subset=strip.subset, 
                                            strip.labels=strip.labels, 
                                            rug.groups=groups,
                                            xlab=effect.label,
                                            set.inf=set.inf,
                                            ...)
            
            level.plot <-   .score_levelplot(   effect=effect, 
                                                ylab=paste( effect.label, "ranks"))
            
            strip.plot <-   
              .score_stripplot(   effect=effect, 
                                  subset=strip.subset, 
                                  groups=groups, 
                                  strip.labels=strip.labels,
                                  ...)
            
            ## generate composite plot
            ##trellis.device( theme = custom.theme(symbol = brewer.pal(9,"Set1"), new=FALSE)  )
            print( density.plot, position=c(0, 0, strip.layout[1], 1), more=TRUE)
            print( level.plot,   position=c(strip.layout[1], 0, strip.layout[1]+strip.layout[2], 1), more=TRUE)
            print( strip.plot,   position=c(strip.layout[1]+strip.layout[2], 0, 1, 1), more=FALSE)
          })

.score_density <- function(effect,
                           strip.subset=1:5,
                           strip.labels=1:5,
                           rug.groups=NULL,
                           xlab="effect",
                           set.inf=20,
                           ...) {

  rug.groups <- factor(rug.groups)
  
  ## replace Inf/-Inf values with set.inf                                                                                                                     
  if( any( identical( abs( effect), Inf ))) {                                                                                                                 
    message(sprintf("Inf/-Inf values were replaced with %s/-%s, as specified by the 'set.inf' parameter.", set.inf, set.inf))                                 
  }                        
  effect[ effect ==  Inf] <-  set.inf                                                                                                                         
  effect[ effect == -Inf] <- -set.inf

  density.plot <- 
      densityplot(effect, xlab=xlab, strip.subset=strip.subset, rug.groups=rug.groups, panel=function(x,...){
        panel.densityplot(x, plot.points=FALSE, col="black", ...)
        panel.rug(x[strip.subset],
                  col=rep( trellis.par.get()$superpose.symbol$col,
                    length.out=length(levels( rug.groups[strip.subset, drop=TRUE])))[ as.integer(rug.groups[ strip.subset, drop=TRUE])]
                  , ...)
        panel.key(levels(rug.groups[strip.subset,drop = TRUE]), corner = c(1,.98),lines = TRUE, points = FALSE)
      })
  return( density.plot )
}  

.score_levelplot <- function(effect,
                             strip.cutoffs=c(-2,2),
                             strip.bounds=c(-4,4),
                             strip.col=c("red","white","blue"),
                             ylab="effect"){
  
  colorFun <- colorRampPalette(c(strip.col[3],strip.col[2],strip.col[1]))
  breaks <- seq(-max(abs(effect)), max(abs(effect)), length.out=256)
  breaks <- breaks[breaks < strip.cutoffs[1] | breaks > strip.cutoffs[2] ] ## blurr small scores
  breaks <- c( -max(abs(effect)), breaks[ breaks > strip.bounds[1] & breaks < strip.bounds[2] ], max(abs(effect)))  ## set scores with min/max colors
  
  level.plot <- levelplot( t(effect[ order( effect ) ]), 
                          scales=list(draw=FALSE), 
                          colorkey=FALSE, 
                          aspect='fill',
                          xlab=NULL, 
                          col.regions=colorFun(256), 
                          at=breaks, 
                          ylab=list(label=ylab),
                          ylim=range(1,length(effect)),
                          par.settings =list(layout.heights=list(top.padding=4, bottom.padding=6))) 
  return( level.plot )
}

.score_stripplot <- function(effect,
                             subset=1:5,
                             groups=NULL,
                             strip.labels=NULL,
                             ylab=NULL,
                             ...) {

  y.range <- range(rank(effect)) * c(-5, 1.01)
  if( !is.null( subset )){
    effect <- rank(effect)[subset]
    groups <- groups[subset,drop=TRUE]
  }

  strip.plot <-
       stripplot(effect ~ groups,
              groups=groups,
              aspect="fill", 
              scales=list(x=list(labels=strip.labels), y=list(draw=FALSE)), 
              ylim=y.range,
              ylab=ylab,
              par.settings=list(layout.heights=list(top.padding=2, bottom.padding=3)),
              ...)

  return( strip.plot )
}
