setMethod(
          "incidence",
          signature( x = "GeneSetCollection" ),
          function (x, ...)  {
            args <- c(x, ...)
            ## extract / create geneSign information
            signList <- lapply(args, function (x) {
              if ( inherits(x, "SignedGeneSet") ) {
                geneSign(x)
              } else {
                factor(rep("up", length(geneIds(x) ) ) )
              }
            }
                               )
            .signed.incidence(lapply(args, geneIds), lapply(args, setName), signList)
          }
          )

setMethod(
          "incidence",
          signature( x = "SignedGeneSet" ),
          function (x, ...)
          {
            args <- c(x, ...)
            .signed.incidence(lapply(args, geneIds), lapply(args, setName), lapply(args, geneSign))
          }
          )

setMethod(
          "incidence",
          signature( x = "CMAPCollection" ),
          function (x, ...)
          {
            Matrix::t( members(x) )
          }
          )

setMethod(
          "incidence",
          signature( x = "list" ),
          function (x, ...)  {
            .incidence(x, names( x ))
          }
          )

.signed.incidence <- function (gidList, gnmList, signList) {
  uids   <- unique(unlist(gidList)) ## all genes
  setids <- unlist(gnmList) ## all sets
  signList <- mclapply( signList, ## all signs
                     function( x ) {
                       ifelse(x == "up", 1, -1) } )
  
  ## assign each geneId its index in the uids vector
  isIn <- lapply( seq( length( gidList ) ),
                 function(n) {
                   isIn <- na.omit(match(gidList[[n]], uids))}
                 )
  
  ## create sparseMatrix with GeneSets in rows and genes in columns
  i <- unlist(lapply(seq( length( isIn ) ), function(n) rep(n, length(isIn[[n]]))))
  j <- unlist(isIn)
  x <- unlist(lapply( seq( length (signList ) ), function(n) signList[[n]] ) )
  sparseMatrix(i=i, j=j, x=x, dims=list(length(setids), length(uids)), dimnames = list(setids, uids) )
}

.incidence <- function (gidList, gnmList)  {
  uids <- unique(unlist(gidList))
  isIn <- mclapply(gidList, function(g, u) match(u, g, nomatch = 0),
                 uids)
  Matrix::t(Matrix(as.integer(unlist(isIn) > 0), ncol = length(gidList),
           dimnames = list(uids, gnmList)))
}
