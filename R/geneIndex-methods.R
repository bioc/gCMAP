setMethod(
  "geneIndex",
  signature( gene.sets = "GeneSet", gene.ids = "character"),
  function( gene.sets, gene.ids, remove.empty=TRUE ) {
    unlist(.map2indices(list(geneIds(gene.sets)), gene.ids, remove.empty=remove.empty))
  }
)

setMethod(
  "geneIndex",
  signature( gene.sets = "GeneSetCollection", gene.ids = "character"),
  function( gene.sets, gene.ids, remove.empty=TRUE ) {
    .map2indices(geneIds(gene.sets), gene.ids, remove.empty=remove.empty)
  }
)

setMethod(
  "geneIndex",
  signature( gene.sets = "CMAPCollection", gene.ids = "character"),
  function( gene.sets, gene.ids, remove.empty=TRUE ) {
    .map2indices(geneIds(gene.sets), gene.ids, remove.empty=remove.empty)
  }
)

.map2indices <- function(...){
  if(exists("symbols2indices", mode = "function")){ ## old limma version
    symbols2indices(...)
  } else if(exists("id2indices", mode = "function")){ ## new limma version
    id2indices(...)
  }
}