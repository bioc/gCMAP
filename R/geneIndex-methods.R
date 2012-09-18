setMethod(
  "geneIndex",
  signature( gene.sets = "GeneSet", gene.ids = "character"),
  function( gene.sets, gene.ids, remove.empty=TRUE ) {
    unlist(symbols2indices(list(geneIds(gene.sets)), gene.ids, remove.empty=remove.empty))
  }
  )

setMethod(
  "geneIndex",
  signature( gene.sets = "GeneSetCollection", gene.ids = "character"),
  function( gene.sets, gene.ids, remove.empty=TRUE ) {
    symbols2indices(geneIds(gene.sets), gene.ids, remove.empty=remove.empty)
  }
  )

setMethod(
  "geneIndex",
  signature( gene.sets = "CMAPCollection", gene.ids = "character"),
  function( gene.sets, gene.ids, remove.empty=TRUE ) {
    symbols2indices(geneIds(gene.sets), gene.ids, remove.empty=remove.empty)
  }
  )
