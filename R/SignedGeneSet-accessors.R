## Inherit from GeneColorSet despite awkward names, to avoid defining lots of
## methods. Force "down" and "up" levels for geneColor, and strip phenotype and
## phenotypeColor information, since it is unused.


setMethod(
          "initialize",
          signature( "SignedGeneSet" ),
          function(
                   .Object,
                   .Template = .Object,
                   ...,
                   geneIds = .Template@geneIds,
                   geneColor = .Template@geneColor,
                   geneSign,
                   phenotype,
                   phenotypeColor
                   )
          {
            ## Handle one extra argument (geneSign), and trap unused slots
            ## (phenotype and phenotypeColor). Handle templating with same class
            ## or superclass, following GSEABase initialize code.
            if ( !missing( geneSign ) )
              geneColor <- geneSign
            if ( length( geneColor ) == 0 || all( geneColor == "" ) )
                geneColor <- rep( "up", length( geneIds ) )
            if ( !is.factor( geneColor ) )
              geneColor <- factor( geneColor, levels = c( "down", "up" ) )            
            callNextMethod(
                           .Object,
                           .Template,
                           ...,
                           geneIds = geneIds,
                           geneColor = geneColor,
                           phenotype = "",
                           phenotypeColor = factor( character( length( geneIds ) ) )
                           )
          }
          )




## Show method

setMethod(
          "show",
          signature( "SignedGeneSet" ),
          function( object ) {
            ## Based on method for GeneColorSet
            cat("setName:", setName(object), "\n")
            cat("geneIds:", paste(selectSome(geneIds(object), maxToShow = 4), 
                                  collapse = ", "), paste("(total: ", length(geneIds(object)), 
                                                          ")\n", sep = ""), sep = " ")
            show(geneIdType(object))
            show(collectionType(object))
            cat(
                "geneSign: ",
                paste( selectSome( as.character( geneSign( object ) ), maxToShow = 4 ), collapse = ", " ),
                "\n  levels: ",
                paste( levels( geneSign( object ) ), collapse = ", " ),
                "\n",
                sep = ""
                )
            cat( "details: use 'details(object)'\n" )
          }
          )




## Constructor

setMethod( "SignedGeneSet", "missing", function( ... ) new( "SignedGeneSet", ... ) )

setMethod(
          "SignedGeneSet",
          "ANY",
          function( type, ..., geneSign, phenotype ) {
            if ( missing( geneSign ) )
              new( "SignedGeneSet", GeneColorSet( type, ..., phenotype = "" ) )
            else
              new( "SignedGeneSet", GeneColorSet( type, ..., phenotype = "" ), geneSign = geneSign )
          }
          )

setAs("CMAPCollection", "SignedGeneSet",
      function( from ) {
        as(as( from, "GeneSet"), "SignedGeneSet")
      })


## Additional ID accessors

setMethod( "upIds", "SignedGeneSet", function( object ) geneIds( object )[ geneSign( object ) == "up" ] )

setMethod( "downIds", "SignedGeneSet", function( object ) geneIds( object )[ geneSign( object ) == "down" ] )




## geneSign accessors

setMethod( "geneSign", "SignedGeneSet", function( obj ) geneColor( obj ) )

setMethod(
          "geneSign<-",
           c( "SignedGeneSet", "factor" ),
          function( object, value ) {
            geneColor( object ) <- value
            return( object )
          }
          )

setMethod(
          "geneSign<-",
           c( "SignedGeneSet", "character" ),
          function( object, value ) {
            geneColor( object ) <- factor( value, levels = c( "down", "up" ) )
            return( object )
          }
          )
