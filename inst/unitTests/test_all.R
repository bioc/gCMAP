test_zScores <- function() {
  checkEqualsNumeric(zScores(1),0)
  checkEqualsNumeric(zScores(0),zScores(.Machine$double.xmin))
}

test_camera <- function() {
  data(gCMAPData)
  gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)
  
  ## random score matrix
  y <- matrix(rnorm(1000*6),1000,6, dimnames=list(featureNames(gCMAPData), 1:6))
  predictor <- c( rep("Control", 3), rep("Case", 3))
  res <- camera_score(y, gene.set.collection, predictor = predictor)

  checkTrue( inherits(res, "CMAPResults"))
}

test_SignedGeneSet <- function() {
  gene.ids <- letters[1:10]
  gene.signs <- rep(c("up","down"), each=5)
  s <- SignedGeneSet(gene.ids, geneSign=gene.signs, setName="set1")

  checkTrue( inherits(s, "GeneSet"))
  checkEquals( as.character( geneSign(s)), rep(c("up","down"), each=5) )
}

test_connectivity <- function() {
  data(gCMAPData)
  
  ## induce CMAPCollection from z-scores
  sets <- induceCMAPCollection(gCMAPData, "z", lower=-3, higher=3)
  res <- connectivity_score(gCMAPData[,1], sets, element="z")
 
  checkTrue( inherits(res, "CMAPResults"))
}

test_fisher <- function() {
  data(gCMAPData)
  gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)
  universe = featureNames(gCMAPData)
  res <-  fisher_score(gene.set.collection[,1], gene.set.collection, universe = universe)
  checkTrue( inherits(res, "CMAPResults"))
  checkEqualsNumeric( effect( res )[1], -Inf)
}

test_geneIndex <- function() {
  data(gCMAPData)
  gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)
  gene.ids <- geneIds(gene.set.collection[,2]) ## geneIds of the second set
  res <- geneIndex(gene.set.collection, gene.ids)
  
  checkTrue( inherits(res, "list"))
  checkEqualsNumeric( length(res),3)
  checkEquals( res[[3]], c(4,10,20,32,43,50))
}

test_generateNChannelSet <- function() {
  data("sample.ExpressionSet") ## from Biobase
  es.list <- list( sample.ExpressionSet[,1:4], sample.ExpressionSet[,5:8], sample.ExpressionSet[,9:12])
  names(es.list) <- paste( "Instance", 1:3, sep=".")
  de <- generate_gCMAP_NChannelSet(
                                   es.list,
                                   1:3,
                                   platform.annotation = annotation(es.list[[1]]),
                                   control_perturb_col="type",
                                   control="Control",
                                   perturb="Case") 

  checkTrue( inherits(de, "NChannelSet"))
  checkEquals( channelNames(de), c("exprs", "log_fc", "p", "z"))
}

test_gsealm_jg <- function() {
  data(gCMAPData)
  gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)
  profile <- assayDataElement(gCMAPData, "z")[,1]
  gsealm_jg_score(profile, gene.set.collection)
  res <- gsealm_jg_score(assayDataElement(gCMAPData, "z")[,1], gene.set.collection)
  
  checkTrue( inherits(res, "CMAPResults"))
  checkEqualsNumeric( padj( res )[1], 0)
}

test_mapNmerge <- function() {
  if( require( "hgu95av2.db" )) {
    data(sample.ExpressionSet) ## from Biobase
    entrez <- mapNmerge(sample.ExpressionSet)
    
    checkTrue( inherits( entrez, "ExpressionSet"))
    checkEqualsNumeric( nrow( entrez), 290, tolerance=0.1 )
    checkEqualsNumeric( ncol( entrez), 26, tolerance=0.1 )
  } else {
    TRUE
  }
}
test_minSetSize <- function() {
  data(gCMAPData)
  gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)
  res <- minSetSize(gene.set.collection, min.members=100)
  
  checkTrue( inherits( res, "CMAPCollection"))
  checkEqualsNumeric( ncol( res ), 1)
  checkEqualsNumeric( nrow( res ), 1000)
}

test_mroast <- function() {
  set.seed(123)
  data(gCMAPData)
  gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)
  y <- matrix(rnorm(1000*6),1000,6, dimnames=list(featureNames(gCMAPData), 1:6))
  predictor <- c( rep("Control", 3), rep("Case", 3))
  res <- mroast_score(y, gene.set.collection, predictor = predictor)
  
  checkTrue( inherits( res, "CMAPResults"))
}

test_romer <- function() {
  set.seed(123)
  data(gCMAPData)
  gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)
  y <- matrix(rnorm(1000*6),1000,6, dimnames=list(featureNames(gCMAPData), 1:6))
  predictor <- c( rep("Control", 3), rep("Case", 3))
  res <- romer_score(y, gene.set.collection, predictor = predictor)
  
  checkTrue( inherits( res, "CMAPResults"))
}

test_wilcox <- function() {
  data(gCMAPData)
  gene.set.collection <- induceCMAPCollection(gCMAPData, "z", higher=2, lower=-2)
  
  profile <- assayDataElement(gCMAPData[,1], "z")
  res <- wilcox_score(profile, gene.set.collection)
  checkTrue( inherits( res, "CMAPResults"))
  checkEqualsNumeric( pval( res )[1], 0, tolerance=1e-100)
}

test_signedRankSum <- function() {
  genes.up   <- c(1:10)
  genes.down <- c(21:30)
  
  set.seed(123)
  scores <- matrix(rnorm(200), ncol=2)
  
  ## the first gene set receives increased /
  ## decreased scores in the first experiment
  scores[genes.up,1]   <- scores[genes.up  ,1] + 1 
  scores[genes.down,1] <- scores[genes.down,1] - 1
  
  ## significantly greater
  res <- signedRankSumTest( statistics = scores[,1],
                    index.up   = genes.up,
                    index.down = genes.down)
  checkEqualsNumeric( res["p"], 7.169693e-08, tolerance=1e-6)

  ## not signficant
  res <- signedRankSumTest( statistics = scores[,2],
                    index.up   = genes.up,
                    index.down = genes.down)
  checkTrue( res["p"] > 0.05)
}

test_featureScores <- function() {  
  data(gCMAPData)
  sets <- induceCMAPCollection(gCMAPData, "z", higher=-2, lower=2)
  res <- featureScores(sets, gCMAPData)
  checkTrue( length(res) == 3, "one matrix per data column is retrieved")
  res <- res[[1]][,1] ## vector for first data column 
  checkTrue( all( res == assayDataElement(gCMAPData, "z")[names(res),1]), "correct scores were extracted")
}
