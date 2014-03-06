signedRankSumTest <- function (statistics, index.up, index.down=NULL,
                               input.is.ranks = FALSE,
                               correlation = 0, df = Inf,
                               adjust.ties=TRUE) {

  if( length( index.up   ) == 0 ) { index.up <- NULL }
  if( length( index.down ) == 0 ) { index.down <- NULL }
  if( is.null(index.up) & is.null(index.down) ) {
    return( c(U=NA, p=1, z=0) )
  } else {
    n <- length(statistics)
    if ( input.is.ranks == FALSE ) {
      r <- rank(statistics)
    } else {
      r <- statistics
    }
    rev.r <- n-r+1
    r1 <- c( r[index.up], rev.r[index.down] )
    n1 <- length(r1)
    n2 <- n - n1
    U <- n1 * n2 + n1 * (n1 + 1)/2 - sum(r1)
    mu <- n1 * n2/2
    if (correlation == 0 || n1 == 1) {
      sigma2 <- n1 * n2 * (n + 1)/12
    }
    else {
      sigma2 <- asin(1) * n1 * n2 + asin(0.5) * n1 * n2 * (n2 - 1)
      + asin(correlation/2) * n1 * (n1 - 1) * n2 * (n2 - 1)
      + asin((correlation + 1)/2) * n1 * (n1 - 1) * n2
      sigma2 <- sigma2/2/pi
    }

    if( adjust.ties == TRUE) {
      TIES <- any( duplicated( r ))
      if (TIES) {
        NTIES <- table(r)
        adjustment <- sum(NTIES * (NTIES + 1) * (NTIES - 1))/(n *
                                                              (n + 1) * (n - 1))
        sigma2 <- sigma2 * (1 - adjustment)
      }
    }

    zlowertail <- (U + 0.5 - mu)/sqrt(sigma2)
    zuppertail <- (U - 0.5 - mu)/sqrt(sigma2)
    p.less <- pt(zuppertail, df = df, lower.tail = FALSE)
    p.greater <- pt(zlowertail, df = df)
    direction <- ifelse( p.less <= p.greater, "less", "greater")

    c(U=U,
      p=ifelse( direction == "less", p.less, p.greater),
      z=ifelse( direction == "less", -zuppertail, -zlowertail)
      )
  }
}

##---------- queries are expression scores

setMethod(
          "wilcox_score",
          signature( experiment = "matrix", sets="CMAPCollection"),
          function( experiment, sets, adjust.ties=FALSE, keep.scores=FALSE, ... ) {

            ## subset objects to shared genes
            matched.features <- match( row.names(experiment), featureNames( sets ))
            matched.sets <- sets[na.omit(matched.features),]
            
            ## extract scores for each gene set
            sets.up     <- mclapply( seq(ncol(matched.sets)),
                                    function( x ) which(members( matched.sets )[ ,x ] == 1 ))
            
            sets.down     <- mclapply( seq(ncol(matched.sets)),
                                      function( x ) which(members( matched.sets )[ ,x ] == -1 ))

            if( keep.scores == TRUE ){
            ## store per-gene scores as data-column:gene-set list-of-list
            gene.scores <- featureScores( experiment, sets)
            } else {
              gene.scores <- NA
            }
            
            ## transform experiment to ranks
            experiment <- apply(experiment, 2, rank)

            ## perform wilcox rank sum test
            stats <- apply( experiment, 2, function( m ) {
              sapply( seq_along( sets.up ), function( n ) {
                res <- signedRankSumTest(m, sets.up[[n]], sets.down[[n]],
                                         adjust.ties=adjust.ties, input.is.ranks=TRUE)
                list(p.value=res["p"],
                     U=res["U"],
                     effect=res["z"]
                     )
              })
            })
            
            ## store results
            results <- mclapply( seq_along( stats ), function(n) {
              
              if( ! all(is.na( gene.scores ))) { 
                geneScores <- I(gene.scores[[n]])
                } else {
                  geneScores <- NA
                }
                
              res <- CMAPResults(
                                 data=data.frame(
                                   set = sampleNames(sets),
                                   trend = ifelse( unlist(stats[[n]]["effect",]) >= 0, "correlated", "anticorrelated"),
                                   pval  = unlist(stats[[n]]["p.value",]),
                                   padj  = p.adjust( unlist( stats[[n]]["p.value",]), method="BH"),
                                   effect = unlist(stats[[n]]["effect",]),
                                   nSet = Matrix::colSums( abs( members (sets)[,sampleNames(sets)] ) ),
                                   nFound = Matrix::colSums( abs( members (matched.sets) ) ),
                                   geneScores = geneScores,
                                   pData(sets))
                                 ,
                                 docs ="\n Results from a two-tailed Wilcox-Rank Sum test \n p-values were adjusted using the 'p.adjust' function with method 'BH'."
                                 )
              
              varMetadata(res)$labelDescription <- c("SetName",
                                                     "Direction",
                                                     "Wilcox Rank Sum test p-value",
                                                     "Adjusted p-value (BH)",
                                                     "z-score",
                                                     "Number of genes annotated in the query set",
                                                     "Number of genes set members found in the dataset",
                                                     "Per-gene raw scores",                                                          
                                                     colnames(pData(sets)))
              res
            })
            names(results) <- colnames(experiment)
            ## return single or list of CMAPResults object(s)
            if( length( results) == 1 ) return( results[[1]] )
            return( results)
          }
          )

setMethod(
          "wilcox_score",
          signature( experiment="numeric", sets="CMAPCollection"),
          function( experiment, sets, ...) {
            wilcox_score(
                         matrix(experiment, ncol=1, dimnames=list(names(experiment), "Profile")),
                         sets
                         )
            
          }
          )

setMethod(
          "wilcox_score",
          signature( experiment="eSet", sets = "CMAPCollection"),
          function( experiment, sets, element = "z", ...) {
            stopifnot( element %in% assayDataElementNames( experiment) )
            wilcox_score( as( assayDataElement(experiment, element), "matrix"), sets, ...)
          }
          )


setMethod(
          "wilcox_score",
          signature( experiment="numeric",sets = "GeneSet"),
          function( experiment,sets,... ) {
            wilcox_score(
                         matrix(experiment, ncol=1, dimnames=list(names(experiment), "Profile")),
                         as( sets, "CMAPCollection"),
                         ...)
          }
          ) 

setMethod(
          "wilcox_score",
          signature( experiment="matrix", sets = "GeneSet"),
          function(experiment, sets,... ) {
            wilcox_score( experiment, as( sets, "CMAPCollection"),...)
          }
          ) 

setMethod(
          "wilcox_score",
          signature( experiment="eSet", sets = "GeneSet"),
          function( experiment, sets, element="z", ... ) {
            stopifnot( element %in% assayDataElementNames( experiment) )
            wilcox_score( as( assayDataElement(experiment, element), "matrix"), as(sets, "CMAPCollection"), ...)
          }
          ) 

setMethod(
          "wilcox_score",
          signature( experiment="numeric", sets = "GeneSetCollection"),
          function( experiment, sets, ... ) {
            wilcox_score(
                         matrix(experiment, ncol=1, dimnames=list(names(experiment), "Profile")),
                         as( sets, "CMAPCollection"),
                         ...)
          }
          ) 

setMethod(
          "wilcox_score",
          signature( experiment="matrix", sets = "GeneSetCollection"),
          function( experiment, sets, ... ) {
            wilcox_score( experiment, as( sets, "CMAPCollection"),...)
          }
          ) 

##-------- queries are gene set(s)

setMethod(
          "wilcox_score",
          signature( experiment = "CMAPCollection", sets="eSet"),
          function( experiment, sets, element="z", adjust.ties=FALSE, keep.scores=FALSE, ... ) {

            stopifnot( element %in% assayDataElementNames(sets) )
 
            ## store raw per-gene expression scores
            if( keep.scores == TRUE ){
            gene.scores =   featureScores(experiment, sets, element=element )
            } else {
              gene.scores <- NA
            }
            score.matrix <- as( assayDataElement( sets, element), "matrix")
            
            ## subset gene sets
            matched.features <- match( featureNames( sets ), featureNames(experiment))
            matched.sets <- experiment[na.omit( matched.features ),]

            ## extract scores for each gene set
            sets.up     <- mclapply( seq(ncol(matched.sets)),
                                 function( x ) which(members( matched.sets )[ ,x ] == 1 ))
            
            sets.down     <- mclapply( seq(ncol(matched.sets)),
                                 function( x ) which(members( matched.sets )[ ,x ] == -1 ))

            ## transform scores to ranks
            score.matrix <- apply(score.matrix, 2, rank)

            ## perform wilcox rank sum test
            stats <- mclapply(seq_along( sets.up ), function( n ) {
               sapply( seq(ncol(score.matrix)), function(m) {
                res <- signedRankSumTest(score.matrix[,m], sets.up[[n]], sets.down[[n]],
                                         adjust.ties=adjust.ties, input.is.ranks=TRUE)
                list(p.value=res["p"],
                     U=res["U"],
                     effect=res["z"]
                  )
              })
            })
            names(stats) <- sampleNames( experiment )

            ## store results
            results <- mclapply( seq_along( stats ), function(n) {
              
              if( ! is.na( gene.scores)) {
                geneScores <- I(gene.scores[[n]])
              } else {
                geneScores <- NA
              }
                res <- CMAPResults(
                                 data=data.frame(
                                   set = sampleNames(sets),
                                   trend = ifelse( unlist(stats[[n]]["effect",]) >= 0, "correlated", "anticorrelated"),
                                   pval = unlist(stats[[n]]["p.value",]),
                                   padj = p.adjust( unlist( stats[[n]]["p.value",]), method="BH"),
                                   effect = unlist(stats[[n]]["effect",]),
                                   geneScores = geneScores,
                                   pData(sets))
                                 ,
                                 docs ="\n Results from a two-tailed Wilcox-Rank Sum test \n p-values were adjusted using the 'p.adjust' function with method 'BH'."
                                 )
              
              varMetadata(res)$labelDescription <- c("SetName",
                                                     "Direction",
                                                     "Wilcox Rank Sum test p-value",
                                                     "Adjusted p-value (BH)",
                                                     "z-score",
                                                     "Per-gene raw scores",
                                                     colnames(pData(sets)))
              res
            })
            names(results) <- sampleNames(experiment)
            ## return single or list of CMAPResults object(s)
            if( length( results) == 1 ) return( results[[1]] )
            return( results)
          }
          )

setMethod(
          "wilcox_score",
          signature( experiment="CMAPCollection", sets="numeric"),
          function( experiment, sets, ...) {
            wilcox_score(
                         experiment,
                         ExpressionSet(
                                       matrix(sets, ncol=1, dimnames=list(names(sets), "Profile"))
                                       ),
                         element="exprs",
                         ...)
          }
          )

setMethod(
          "wilcox_score",
          signature( experiment="CMAPCollection", sets = "matrix"),
          function( experiment, sets, ...) {
            wilcox_score( experiment, ExpressionSet(sets), element="exprs", ...)
          }
          )


setMethod(
          "wilcox_score",
          signature( experiment="GeneSet",sets = "numeric"),
          function( experiment,sets,  ... ) {
            wilcox_score(
                         as( experiment, "CMAPCollection"),
                         ExpressionSet(
                                       matrix(sets, ncol=1, dimnames=list(names(sets), "Profile"))
                                       ),
                         element="exprs",
                         ...)
          }
          ) 

setMethod(
          "wilcox_score",
          signature( experiment="GeneSet", sets = "matrix"),
          function(experiment, sets, ... ) {
            wilcox_score( as( experiment, "CMAPCollection"), ExpressionSet(sets), element="exprs",...)
          }
          ) 

setMethod(
          "wilcox_score",
          signature( experiment="GeneSet", sets = "eSet"),
          function( experiment, sets,element="z",... ) {
            wilcox_score( as(experiment, "CMAPCollection"), sets,...)
          }
          ) 

setMethod(
          "wilcox_score",
          signature( experiment="GeneSetCollection", sets = "numeric"),
          function( experiment, sets, ... ) {
            wilcox_score(
                         as( experiment, "CMAPCollection"),
                         ExpressionSet(
                                       matrix(sets, ncol=1, dimnames=list(names(sets), "Profile"))
                                       ),
                         element="exprs",
                         ...)
          }
          ) 

setMethod(
          "wilcox_score",
          signature( experiment="GeneSetCollection", sets = "matrix"),
          function( experiment, sets, ... ) {
            wilcox_score( as( experiment, "CMAPCollection"),
                         ExpressionSet(sets),
                         element="exprs",...)
          }
          ) 

setMethod(
          "wilcox_score",
          signature( experiment="GeneSetCollection", sets = "eSet"),
          function( experiment, sets,element="z",... ) {
            wilcox_score( as( experiment, "CMAPCollection"),
                         sets,
                         ...)
          }
          ) 
