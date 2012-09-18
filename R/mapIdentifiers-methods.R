setMethod(
          "mapIdentifiers",
          signature( what = "SignedGeneSet", to="GeneIdentifierType", from="GeneIdentifierType"),
          function (what, to, from, ... ) 
          {
            if (to@type != "EntrezId")
              {
                stop("Currently, only conversion to EntrezIds is supported for SignedGeneSets")
              }

            ## map identifiers
            if (grepl("IGIS", annotation(geneIdType(what)), ignore.case = TRUE)) { ## Genentech IGIS gene annotation
              ids <- sub("GeneID:", "", geneIds(what))
              geneIds(what) <- ids
              return(what)

            } else { ## other platform annotation
              ids <- .lookupAnno(geneIds(what), annotation(from))
              message(sprintf("Successfully mapped gene set members to %s Entrez id(s).", length(unique(ids))))
              
              
              ## remove multi-mapping probes with inconsistent signs
              sign <- as.character(geneSign(what))
              id.sign <- unique(cbind(ids, sign))
              consistent <- id.sign[,1] %in%  names(which(table(id.sign[,1]) == 1))
              message(sprintf("Discarded %s Entrez id(s) because of sign inconsistency.",
                              length(unique(ids))-length(id.sign[consistent,1])))
            }

            ## return new SignedGeneSet
            SignedGeneSet(id.sign[consistent,1],
                          geneSign=id.sign[consistent,2],
                          setName=setName(what),
                          geneIdType=to)
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
