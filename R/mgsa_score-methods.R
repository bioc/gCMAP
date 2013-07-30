if( ! isClass("MgsaMcmcResults")){
	setClass(
        "MgsaMcmcResults"
        )
}

setMethod(
  "mgsa_score",
  signature( query = "GeneSet", sets = "CMAPCollection", universe = "character" ),
  function( query, sets, universe, keep.scores=FALSE,
            alpha=seq(0.01,0.3, length.out=10), 
            beta=seq(0.1,0.8, length.out=10),
            p=seq( min(0.1, 1/length(sets)), 
                   min(0.3, 20/length(sets)), length.out=10),
            steps=1e6, restarts=5, 
            threads=getOption("mc.cores", default=0)) {
    
    if(!suppressWarnings(require("mgsa", quietly=TRUE, character.only=TRUE))){
      stop("To use this method, please install the Bioconductor package 'mgsa'.")
    }
    
    if( length( geneIds(query )) == 0){
      stop("Query sets is empty.")
    }
    collection.ids <- geneIds( sets )
    res <- mgsa::mgsa( geneIds( query), collection.ids, 
                       population=universe, alpha=alpha,
                       beta=beta, p=p, steps=steps,
                       restarts=restarts, threads=threads)   
    res <- as( res, "CMAPResults")
    res@data <- data.frame( res@data, pData( sets ), stringsAsFactors=FALSE)
    
    ## store per-gene scores as data-column:gene-set list-of-list
    if( keep.scores == TRUE) {
      res@data$geneScore <- featureScores( as(query, "CMAPCollection"), sets)
    }
    return( res )
  }
)

setMethod(
  "mgsa_score",
  signature( query = "GeneSet", sets = "NChannelSet", universe = "character" ),
  function( query, sets, universe, element, lower=NULL, higher=NULL, 
            keep.scores=FALSE, min.set.size=5, ...) {
    
    if(all( is.null( c(lower, higher )))){
      stop("Please provide at least one of the 'higher' and 'lower' parameters.")
    }
    
    ## induce CMAPCollection from NChannelSet
    induced.sets <- induceCMAPCollection( sets, element=element, 
                                          lower=lower, higher=higher)
    if( ncol(induced.sets) == 0){
      stop("None of the genes in the reference dataset passed the score cutoff to induce gene sets.")
    }
    
    if( !is.null( min.set.size )){
      induced.sets <- minSetSize(induced.sets, min.members = min.set.size)
      if( ncol(induced.sets) == 0){
        stop(sprintf("None of the induced gene sets had more than %s members.", min.set.size))
      }
    }
    
    results <- mgsa_score( query, induced.sets, universe=universe, keep.scores=FALSE, ...)
    
    ## store per-gene scores as data-column:gene-set list-of-list
    query  <- query[intersect(geneIds( query), universe),]
    sets   <- sets[ intersect(featureNames( sets ),  universe), sampleNames( induced.sets)]
    
    ## store per-gene scores as data-column:gene-set list-of-list
    if( keep.scores == TRUE) {
      results@data$geneScores <- featureScores( query, sets, simplify=FALSE)
    }
    return( results )
  }
)

setMethod(
  "mgsa_score",
  signature( query = "SignedGeneSet", sets = "CMAPCollection", universe = "character" ),
  function( query, sets, universe, ...) {
    mgsa_score(as(query,"GeneSet"), sets, universe=universe, ...)
  }
)

setMethod(
  "mgsa_score",
  signature( query = "SignedGeneSet", sets = "NChannelSet", universe = "character" ),
  function( query, sets, universe, ...) {
    mgsa_score(as(query,"GeneSet"), sets, universe=universe, ...)
  }
)


setMethod(
  "mgsa_score",
  signature( query = "GeneSetCollection", sets = "CMAPCollection", universe = "character" ),
  function( query, sets, universe, ... ) {
    res <- lapply( query, function( x ){
      if( length( geneIds( x )) == 0){
        NULL
      } else {
        mgsa_score(x, sets, universe, ...)
      }
    })
    names( res ) <- sapply( query, setName)
    return( res )
  })

setMethod(
  "mgsa_score",
  signature( query = "GeneSetCollection", sets = "NChannelSet", universe = "character" ),
  function( query, sets, universe, ... ) {
    res <- lapply( query, function( x ){
      if( length( geneIds( x )) == 0){
        NULL
      } else {
        mgsa_score(x, sets, universe, ...)
      }
    })
    names( res ) <- sapply( query, setName)
    return( res )
  })

setMethod(
  "mgsa_score",
  signature( query = "GeneSetCollection", sets = "GeneSetCollection", universe = "character" ),
  function( query, sets, universe, ... ) {
    mgsa_score(query, as(sets, "CMAPCollection"), universe=universe, ...)
  }
)

setMethod(
  "mgsa_score",
  signature( query = "CMAPCollection", sets = "GeneSetCollection", universe = "character" ),
  function( query, sets, universe, ... ) {
    mgsa_score(as(query, "GeneSetCollection"), 
               as(sets, "CMAPCollection"), universe=universe, ...)
  }
)

setMethod(
  "mgsa_score",
  signature( query = "GeneSet", sets = "GeneSetCollection", universe = "character" ), ## sets from cmap
  function( query, sets, universe, ...) {
    mgsa_score(query, as(sets, "CMAPCollection"), universe, ...)
  }
)


setMethod(
  "mgsa_score",
  signature( query = "GeneSet", sets = "GeneSet", universe = "character" ),
  function( query, sets, universe, ... ) {
    mgsa_score(query, as(sets, "CMAPCollection"), universe, ...)
  }
)

setMethod(
  "mgsa_score",
  signature( query = "CMAPCollection", sets = "CMAPCollection", universe = "character" ),
  function( query, sets, universe, ... ) {
    if( ncol( query ) >1 ){
      mgsa_score(as(query, "GeneSetCollection"), sets, universe)
    } else { 
      mgsa_score(as(query, "GeneSet"), sets, universe, ...)
    }
  }
)