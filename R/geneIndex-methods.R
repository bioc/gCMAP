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
  if( packageVersion("limma") >= "3.20.0"){
    ids2indices(...)
  } else {
    symbols2indices(...)
  }
}
