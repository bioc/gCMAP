setMethod(
          "mapIdentifiers",
          signature( what = "SignedGeneSet", to="GeneIdentifierType", from="GeneIdentifierType"),
          function (what, to, from, ... ) 
          {
            if( to@type == from@type) return( what )

            ## split SignedGenSet into two separate GeneSets for each sign
            split.ids <- split(geneIds( what), geneSign(what))
            split.sets <- lapply( split.ids, function( ids ){
              GeneSet( ids, geneIdType=from)
            })
            split.sets.new <- lapply( split.sets, mapIdentifiers, to)
            new.ids <- do.call( c, lapply( split.sets.new, geneIds))
            names( new.ids ) <- ifelse( grepl( "down", names( new.ids)), "down","up")

            ## remove duplicates ( same id & same sign)
            id.sign <- cbind( as.character( new.ids ), names( new.ids))
            id.sign <- id.sign[!duplicated(id.sign, MARGIN=1),]
            new.ids <- id.sign[,1]
            names(new.ids) <- id.sign[,2]

            ## remove genes appearing both up- and down-regulated
            inconsistent <- sum(duplicated( new.ids))
            if( inconsistent != 0){
              message(paste("Removed", inconsistent, "gene(s) that appeared both up- and downregulated."))
            }
            new.ids <- new.ids[ !duplicated( new.ids) ]

            ## return new SignedGeneSet
            SignedGeneSet(new.ids,
                          geneSign=names(new.ids),
                          setName=setName(what),
                          geneIdType=geneIdType( split.sets.new[[1]])
                          )
          }
          )

setMethod(
          "mapIdentifiers",
          signature( what = "SignedGeneSet", to="GeneIdentifierType", from="missing"),
          function (what, to, ...) 
          {
            if ( annotation(geneIdType(what)) == "" )
              {
                stop("SignedGeneSet does not include annotation information.\nPlease provide 'from' parameter.")
              }
            from = geneIdType(what)
            mapIdentifiers(what, to, from)
          }
          )

setMethod(
          "mapIdentifiers",
          signature( what = "SignedGeneSet", to="missing", from="ANY"),
          function(what, ...)
          {
            stop("Please specify the 'to' parameter.")
          }
          )
