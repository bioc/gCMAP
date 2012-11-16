setClass("CMAPCollection",
         contains = "eSet",
         prototype = prototype(
           new("VersionedBiobase",
               versions=c(classVersion("eSet"), CMAPCollection="1.0.0")),
           phenoData = new("AnnotatedDataFrame",
             data=data.frame(),
             varMetadata=data.frame(
               labelDescription=character(0),
               channelDescription=factor())))
         )

setClass("CMAPResults",
         contains="AnnotatedDataFrame",
         representation(
           docs = "character",
           errors = "list"
         ),
         prototype = prototype(
           docs = "",
           errors = list()
           )         
         )

setClass(
         "SignedGeneSet",
         contains = "GeneColorSet",
         validity = function( object ) {
           if (! all(levels( object@geneColor ) %in% c( "down", "up" ) ) )
             return( "Levels for 'geneColor'/'geneSign' factor must be 'down' and/or 'up'." )
           if ( any( is.na( object@geneColor ) ) )
             return( "All 'geneColor'/'geneSign' entries must be 'down' or 'up'." )
           return( TRUE )
         }
         )

setClassUnion(
              "matrix_or_big.matrix",
              c("matrix", "big.matrix")
              )
