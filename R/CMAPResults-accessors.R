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

setAs("MgsaMcmcResults", "CMAPResults",
      function (from) {
        res <- gCMAP:::CMAPResults(
          data=data.frame( 
            set=row.names( from@setsResults),
            trend=sapply( from@setsResults$estimate,
                          function(x){
                ifelse( !is.na(x) && x < 1, "over", "under")
              }),
            pval=sapply( from@setsResults$estimate,
                         function(x){
                           ifelse( is.na(x), 1, 1-x)
                         }),
            padj=sapply( from@setsResults$estimate,
                         function(x){
                           ifelse( is.na(x), 1, 1-x)
                         }),
            effect=sapply( from@setsResults$estimate, 
                           function(x) {
                             ifelse( is.na(x), 0, zScores(1-x)) 
                             }),
            nSet=from@setsResults$inPopulation,
            nFound=from@setsResults$inStudySet
          ),
          docs="\n Results from a Bayesian MGSA analysis.\nThe reported p-values represent '1-marginal posterior probability'."
          )
        
        varMetadata(res)$labelDescription <- 
          c("SetName",
            "Deviation from random expectation",
            "1-estimate of the marginal posterior probability",
            "1-estimate of the marginal posterior probability",
            "z-score based on the standard normal distribution",
            "Number of genes annotated in the reference set",
            "Number of genes found in query and target sets"
          )
        res
      }
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
setMethod("plot", "CMAPResults", function(x,
                                          strip.effect="effect",
                                          strip.pval="padj",
                                          set.inf=20,
                                          main="Distribution of similarity scores",
                                          col.up="red",
                                          col.down="blue",
                                          strip.cutoffs=c(-3,3),
                                          strip.bounds=c(-6,6),
                                          strip.col=c("red","white","blue"),
                                          ...
                                          ) {

  ## extract (sorted) data.frame
  dat <- cmapTable( x )
  
  ## retreive or calculate effect
  if( is.null( strip.effect) ) {
    effect.population <- zScores( dat[,strip.pval] ) ## convert to z-scores from normal distribution
    effect.label <- "z-score"

  } else if( strip.effect %in% varLabels( x ) ) {
    effect.population <- dat[,strip.effect] ## retrieve effect sizes
    effect.label <- labels( x )[strip.effect,]

  } else if ( strip.pval %in% varLabels( x ) ){ ## retrieve pvalues
    message(sprintf("Column %s not found in CMAPResults object: Transforming %s to unsigned z-scores.", strip.effect, strip.pval))
    ## calculate zscores from pvalues
    effect.population <- zScores(dat[,strip.pval]) ## convert to z-scores from normal distribution
    effect.label <- "z-score"

  } else { ## neither effects nor pvalues are found: fail
    stop(sprintf("Neither column %s nor %s found in CMAPResults object.", strip.effect, strip.pval))
  }

  effect.sample <- effect.population[ which(effect.population > max( strip.cutoffs) | effect.population < min(strip.cutoffs))]

  plot.new()
  par(fig=c(0,0.9,0,1), new=TRUE)
  .overview.density( effect.sample, effect.population, main=main, xlab=effect.label, col.set="black", col.up=col.up, col.down=col.down)
  par(fig=c(0.85,1,0,1),mar=c(5, 0, 4, 2) + 0.1, new=TRUE )
  .overview.heatmap( effect.population, strip.cutoffs=strip.cutoffs, strip.bounds=strip.bounds, strip.col=strip.col, set.inf=set.inf )
}
          )


.overview.density <- function(effect.sample, effect.population, main="Distribution of similarity scores", xlab="Effect size", col.set="black", col.up="red", col.down="blue"){

  population.density <- density( effect.population, na.rm =TRUE )

  ## plot densities
  plot(population.density, col="lightgrey", lty=1, type="l", xlab=xlab, main=main,
       ylim=c( min( population.density$y), max( population.density$y, 0.4) )
       )
  polygon( population.density, col="lightgrey", border=NA)
  ## add Gaussian
  x=seq(min( population.density$x), max( population.density$x),length=200)
  lines( x,
        y=dnorm(x,mean=0,sd=1),
        type="l",lty=2,col="darkgrey")

  ## add rug
  rug( effect.sample[effect.sample <= 0], ticksize = 0.05, col=col.down, lwd=1)
  rug( effect.sample[effect.sample > 0], ticksize = 0.05, col=col.up, lwd=1)

  ## add legend
  legend("topright", c("All scores", "Correlated instance", "Anti-correlated instance", "Normal distr."), col=c("lightgrey", col.up, col.down, "darkgrey"), lty=c(1,1,1,2), bg=NULL, cex=0.5)
}

.overview.heatmap <- function(effect.population,
                              strip.cutoffs=c(-3,3),
                              strip.bounds=c(-6,6),
                              strip.col=c("red","white","blue"),
                              set.inf=20
                              ){
  effect.population[is.infinite( effect.population )] <- set.inf * sign( effect.population[is.infinite( effect.population )] )
  min.score <- min( strip.bounds[1], -max(abs(effect.population)))
  max.score <- max( strip.bounds[2], max(abs(effect.population)))
  breaks <- seq(min.score, max.score, length.out=256)
  breaks <- breaks[breaks < strip.cutoffs[1] | breaks > strip.cutoffs[2] ] ## blurr small scores
  breaks <- c( min.score, breaks[ breaks > strip.bounds[1] & breaks < strip.bounds[2] ], max.score )  ## set scores with min/max colors
  map.colors <- colorRampPalette(c(strip.col[3],strip.col[2],strip.col[1]))(length(breaks)-1)
  image(matrix(sort( effect.population), nrow=1),main="", ylab="", xlab="", xaxt="n", yaxt="n", breaks=breaks, col=map.colors)
}
