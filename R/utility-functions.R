eSetOnDisk <- function(eset, out.file=NULL) {
  if( is.null(out.file)) {
   stop("Please provide path/basename of the output file(s).")
  }
  
  storageMode(assayData( eset )) <- "list" ## unlock
  
  for ( element.name in assayDataElementNames( eset ) ) {
    element <- assayDataElement( eset, element.name )
    bm <- BigMatrix(element,
                      backingfile = paste( out.file, element.name, sep="_"),
                      )
    
    assayData( eset )[[element.name]] <- bm
  }
  
  storageMode( assayData( eset ) ) <- "lockedEnvironment" ## lock
  assign( basename(out.file), eset)
  save( list=basename(out.file), file=paste(out.file, "rdata", sep="."))
  eset
}

eSetProbeToGene <- function(eset, get="ENTREZID", channel=NULL, genes=NULL, na.rm=TRUE) {
  gene.ids <- .probeAnno(eset, get=get)

  if ( all(is.na(gene.ids))) {
    stop("None of the feature ids could be mapped !")

  } else {
    message(
          paste(
                length(
                       gene.ids[!is.na(gene.ids)]
                       ),
                "of",
                dim(eset)[1],
                "feature ids were successfully mapped to",
                length(unique(na.omit(gene.ids))),
                paste(get,"s", sep="")
                )
          )
  }      

  ## subset eSet if requested
  if ( !is.null( genes ) ) {
    matched.features <- which(gene.ids %in% genes)
    eset <- eset[ matched.features,]
    gene.ids <- gene.ids[matched.features] 
  }
  
  ## create empty eSet for output
  annotation.types <- list(ENTREZID="EntrezId", SYMBOL="Symbol") ## match GSEABase nomeclature
  gene.eset <- new("NChannelSet") ## eset[,FALSE] ## metaData only
  pData(gene.eset) <- pData(eset)
  featureData(gene.eset) <- new("AnnotatedDataFrame")
  annotation(gene.eset) <- annotation.types[[get]]
  
  ## for probes mapped to the same gene identifier, apply 'fun'
  if( is.null(channel) )
    {
      channels <- assayDataElementNames(eset)
    } else {
      channels <- intersect(channel, assayDataElementNames(eset))
      message(paste("Retaining channels:", paste(channels, collapse=", ")))
    }
  
  for (element in channels ) {
    assayDataElement(gene.eset,element) <- matrix(
                                                  apply(assayDataElement(eset, element),
                                                        2,
                                                        function(column) {
                                                          ave(column,
                                                              gene.ids,
                                                              FUN=function(x)mean(x, na.rm=na.rm))
                                                        }
                                                        ),
                                                  ncol=ncol(eset),
                                                  dimnames=list(
                                                    gene.ids,
                                                    sampleNames(eset)),
                                                  )
    
  }
  gene.eset[! (is.na(featureNames(gene.eset)) | duplicated(featureNames(gene.eset))),]
}


.lookupAnno <- function(gene.ids, platform, get="ENTREZID") {
  unlist(
         lookUp(
                gene.ids,
                platform,
                what=get,
                load=TRUE
                )
         )
}

.probeAnno <- function(eset, get="ENTREZID") {

  if ( length(annotation(eset)) == 0 ) {
    stop("eSet does not contain annotation information.")
  }

  ## IGIS annotation
  if ( grepl("igis", annotation(eset), ignore.case=TRUE ) )
    {
      if( get != "ENTREZID" ) {
        stop("Currently IGIS ids can only be converted to ENTREZIDs.")
      } else {
        ## for IGIS ids, simply drop the prefix      
        return( sub("GeneID:", "", featureNames(eset)) )
      }
      
    } else {
      ## for other platforms, use bioconductor annotation package
      .lookupAnno(featureNames(eset), annotation(eset), get=get)
    }
}

addFeatureData <- function(eset) {
  entrez.table <- toTable(get(paste(annotation(eset), "ENTREZID", sep = "")))
  symbol.table <- toTable(get(paste(annotation(eset), "SYMBOL", sep = "")))
  conversion.table <- merge(entrez.table, symbol.table)
  colnames(conversion.table) <- c("Annotation", "Entrez", "Symbol")
  featureData(eset) <- as(conversion.table, "AnnotatedDataFrame")
  eset
}

foldChange <- function(vst, perturbations, controls) {
  ## calculates the log2 fold changes for each instance by subtracting
  ## the rowMeans of perturbation and control samples
  ## perturbations = sample names of perturbation columns
  ## controls = sample names of control columns

  mod_fc <- rowMeans( subset(vst, select=perturbations) ) - rowMeans( subset(vst, select=controls) )
  if( length(mod_fc) == 0) {
    return(as.numeric(rep(NA, dim(vst)[1])))
  } else {
    return(mod_fc)
  }
}
##----------------------------------------
## standard t-test for array data
##----------------------------------------
pairwise_compare <- function( eset,
                             control_perturb_col = "cmap", 
                             control="control",
                             perturb="perturbation") {
  
  ## Immediately fail for cases that have no ExpressionSet
  if(!is( eset, "ExpressionSet" ))	
    {
      stop(paste("Expected ExpressionSet found: ", class(eset)))
    }
  
  ## Check column 'control_perturb_col' exists
  if(!control_perturb_col %in% varLabels(eset))
    {
      stop(paste("Column label <", control_perturb_col, "> not found "))
    }
  
  ## Check control and perturb exist
  ncontrol = length( which(pData(eset)[,control_perturb_col] == control)) 
  nperturb = length( which(pData(eset)[,control_perturb_col] == perturb)) 
    
  if(ncontrol == 0 | nperturb == 0)
    {
      stop(paste("Expected at least one control (found ", ncontrol, ") ",
                 "and a at least one perturb (found ", nperturb, ")"))
    }

  if(ncontrol + nperturb != length(sampleNames(eset)))
  {
    warning(paste("ExpressionSet includes ", 
                  length(sampleNames(eset)) - ncontrol - nperturb, 
                  " non-control/non-perturb samples, these will be ingored"))
    eset <- eset[,pData(eset)[,control_perturb_col] %in% c(control, perturb)]
  }
  
  ## Compare control and perturb
  if(ncontrol == 1 & nperturb == 1)
    {
      AveExpr = rowMeans(exprs(eset), na.rm=TRUE)
      log_fc = exprs(eset)[,which(pData(eset)[,control_perturb_col] == perturb)] - 
        exprs(eset)[,which(pData(eset)[,control_perturb_col] == control)]
      p = rep(NA_real_, length(featureNames(eset)))
      z = rep(NA_real_, length(featureNames(eset)))
    }
  else
    {
      comparison <- factor(pData(eset)[,control_perturb_col], levels=c(perturb, control), ordered=TRUE)
      t_res= rowttests(exprs(eset), fac=comparison)
      AveExpr = rowMeans(exprs(eset), na.rm=TRUE)
      log_fc = t_res$dm  ## perturb - control
      p = t_res$p.value
      z =  zScores(p, log_fc)
    }
  
  data.frame(exprs=AveExpr, log_fc, p, z, row.names = featureNames(eset))  
}

##----------------------------------------
## limma's moderated t-test for array data
##----------------------------------------
pairwise_compare_limma <- function( eset,
                                   control_perturb_col = "cmap", 
                                   control="control",
                                   perturb="perturbation",
                                   limma.index=2)
{
  ## Immediately fail for cases that have no ExpressionSet
  if( ! is( eset, "ExpressionSet" ) )	
    {
      stop( paste( "Expected ExpressionSet found: ", class( eset ) ) )
    }
  
  ## Check column, control_perturb_col, exists
  if( !control_perturb_col %in% varLabels( eset ) )
    {
      stop( paste( "Column label <", control_perturb_col, "> not found " ) )
    }
  
  ## Check control and perturb exist
  ncontrol = length( which( pData( eset )[,control_perturb_col] == control ) ) 
  nperturb = length( which( pData( eset )[,control_perturb_col] == perturb ) ) 
  
  if( ncontrol + nperturb != length( sampleNames( eset ) ) )
    {
      warning( paste( "ExpressionSet includes ", 
                     length( sampleNames( eset ) ) - ncontrol - nperturb, 
                     " non-control/non-perturb samples, these will be ingored" ) )
    }
  
  if( ncontrol == 0 | nperturb == 0 )
    {
      stop( paste( "Expected at least one control ( found ", ncontrol, " ) ",
                  "and a at least one perturb ( found ", nperturb, " )" ) )
      eset <- eset[,pData( eset )[,control_perturb_col] %in% c( control, perturb )]
    }

  ## Differential expression analysis
  if( ncontrol == 1 & nperturb == 1 ) ## no replicates, no statistical testing
    {
      AveExpr = rowMeans( exprs( eset ), na.rm=TRUE )
      log_fc = exprs( eset )[,which( pData( eset )[,control_perturb_col] == perturb )] - 
        exprs( eset )[,which( pData( eset )[,control_perturb_col] == control )]
      p = rep( NA_real_, length( featureNames( eset ) ) )
      z = rep( NA_real_, length( featureNames( eset ) ) )
      
    } else {
      design <- model.matrix( ~factor( pData( eset )[[control_perturb_col]], level=c( control, perturb ) ) )
      fit <- eBayes( lmFit( eset, design ) )

      ## Check that the residual degress of freedom and stdev.unscaled do not vary.
      if ( !all( fit$df.residual == fit$df.residual[1] ) )
        stop( "Code assumes that 'df.residual' is not gene specific. Was there missing data?" )
      
      if ( !all( fit$stdev.unscaled[,limma.index] == fit$stdev.unscaled[1,limma.index] ) )
        stop( "Code assumes that 'stdev.unscaled' is the same for all genes. Was there missing data?" )
      
      ## Extract relevant columns from limma results
      topT <- topTable( fit, coef=limma.index, sort.by="none", number=Inf )
      
      AveExpr = topT$AveExpr
      log_fc   = topT$logFC
      p        = topT$P.Value
      z         = zScores( p, log_fc )      
    }
  ## return results as a data.frame
  data.frame(exprs=AveExpr,
             log_fc,
             p,
             z,
             row.names = featureNames( eset )
             ) 
}


##----------------------------------------
## DESeq's nbinomTest for count data
##----------------------------------------
pairwise_DESeq <- function( cds,
                           vst,
                           control_perturb_col = "condition", 
                           control="control",
                           perturb="perturbation",
                           try.hard=FALSE)
{
  if( is.element("DESeq", installed.packages()[,1])){
    require( DESeq )
  } else {
    stop("To run this function, please install the Bioconductor package 'DESeq'.")
  }
  
  ## Immediately fail for cases that have no CountDataSet
  if(!is( cds, "CountDataSet" ))	
    {
      stop( paste( "Expected CountDataSet, found: ", class( cds )))
    }
  
  ## Check column, control_perturb_col, exists
  if( ! control_perturb_col %in% varLabels( cds ) )
    {
      stop( paste( "Column label <", control_perturb_col, "> not found "))
    }
  
  ## Check control and perturb exist
  ncontrol = length( which( pData( cds )[, control_perturb_col ] == control )) 
  nperturb = length( which( pData( cds )[, control_perturb_col ] == perturb )) 
  
  if( ncontrol + nperturb != length( sampleNames( cds ) ) )
    {
      warning( paste( "CountDataset includes ", 
                    length( sampleNames( cds ) ) - ncontrol - nperturb, 
                    " non-control/non-perturb samples, these will be ingored"))
      cds <- cds[, pData( cds )[, control_perturb_col ] %in% c( control, perturb )]
    }
  
  if( ncontrol == 0 | nperturb == 0 )
    {
      stop( paste( "Expected at least one control (found ", ncontrol, ") ",
                 "and a at least one perturb (found ", nperturb, ")"
                  )
           )
    }

  ## Compare control and perturb
  design <- model.matrix( ~factor( pData( cds )[[ control_perturb_col ]],
                                 level = c( control, perturb ) ) )
  res <- try( .DESeq_nbinom(cds,
                           control,
                           perturb,
                           try.hard,
                           control_perturb_col=control_perturb_col)
             )
  
  if( class(res) != "try-error" ) {
    ## compute moderated fold-changes from variance-transformed counts (vst)
    res$mod_logFC <- .moderated_logFC(vst, cds, control_perturb_col, control, perturb)
  }

  p        = res$pval
  z        = zScores( p, res$log2FoldChange )

  p[ is.na( p ) ] <- 1 ## replace NA with p-value 1
  z[ is.na( z ) ] <- 0 ## replace NA with z-score 0

  data.frame(
             exprs  = res$baseMean,
             log_fc = res$log2FoldChange,
             mod_fc = res$mod_logFC,
             p = p,
             z = z,
             row.names = featureNames( cds )
             )
}

##----------------------------------------
## calulation of moderated log2 fold
## changes from variance-transformed data
##----------------------------------------
.moderated_logFC <- function(vst, cds, control_perturb_col, control, perturb) {

  ## identify indices of control and perturbation columns
    perturbations <- sampleNames(cds)[which(pData(cds)[,control_perturb_col] == perturb)]
    controls <- sampleNames(cds)[which(pData(cds)[,control_perturb_col] == control)]

    ## calculate fold changes
    mod_fc <- .foldChange( vst, perturbations, controls )  
    mod_fc <- rowMeans( subset( vst, select=perturbations) ) - rowMeans( subset( vst, select=controls ) )
    if( length(mod_fc) == 0) {
      mod_fc <- as.numeric(rep(NA, dim( vst )[ 1 ] ) )
    }
    return( mod_fc )
  }

.foldChange <- function(vst, perturbations, controls) {
  ## perturbations = sample names of perturbation columns
  ## controls = sample names of control columns
  mod_fc <- rowMeans( subset(vst, select=perturbations) ) - rowMeans( subset(vst, select=controls) )
  if( length(mod_fc) == 0) {
    return(as.numeric(rep(NA, dim(vst)[1])))
  } else {
    return(mod_fc)
  }
}

##----------------------------------------
## nbinomTest of count data
##----------------------------------------
.DESeq_nbinom <- function(cds,
                         control="control",
                         perturb="perturbation",
                         try.hard=FALSE,
                         control_perturb_col="cmap",
                         ...) {
  if( is.element("DESeq", installed.packages()[,1])){
    require( DESeq )
  } else {
    stop("To use this function, please install the Bioconductor package 'DESeq'.")
  }
  
  conditions(cds) <- pData(cds)[, control_perturb_col]
  cds <- DESeq::estimateSizeFactors( cds )
  
  ## can we use replicates to estimate dispersions ?
  if ( max( table( conditions(cds ) )[ c( control, perturb )] ) > 1 ) {
    
    ## yay, replicates for at least one condition !
    cds <- try( DESeq::estimateDispersions( cds ) )
  } else {
    ## nay, no replicates
    cds <- try( DESeq::estimateDispersions( cds, method = "blind", sharingMode = "fit-only", ... ) )
  }
  
  if( class( cds ) == "try-error" ) {
    if( try.hard == TRUE ) {
      warning( "DESeq analysis with parametric dispersion estimate failed. Trying local fit instead." )
      
      ## try local fit instead
      cds <- try( .DESeq_nbinom( cds, try.hard=FALSE, fitType="local" ) )
    } else {
      stop( "Error: DESeq dispersion estimate failed." )
    }
  }
  
  res <- DESeq::nbinomTest( cds, control, perturb )
  res
}

generate_gCMAP_NChannelSet <- function(
                                       data.list,
                                       uids=1:length(data.list),
                                       sample.annotation=NULL,
                                       platform.annotation="",
                                       control_perturb_col="cmap",
                                       control="control",
                                       perturb="perturbation",
                                       limma=TRUE,
                                       limma.index=2,
                                       big.matrix=NULL,
                                       center.z="peak",
                                       center.log_fc="peak"
                                       )
{

  ## Check that sample.annotation is complete or NULL
  if ( !is.null(sample.annotation) )
    if ( !is.data.frame(sample.annotation) || nrow(sample.annotation) != length(data.list) )
      stop("sample.annotation should be a data frame with one row per element of data.list")
  
  ## Check that all data elements are of the same class
  data.classes = unique( unlist( lapply( data.list, class) ) )
  
  if(length(data.classes) > 1) {
    stop("All data.list elements must be of the same class (ExpressionSet or CountDataSet).")
  }

  if(! data.classes %in% c("ExpressionSet", "CountDataSet")) {
    stop("'Data.list' must be a list of either ExpressionSet or CountDataSet objects")
  }
  
  ## Check featureNames of expression sets are identical
  fnames = lapply(data.list, featureNames)
  test = unique(unlist(lapply(fnames, all.equal, fnames[[1]], check.attributes=FALSE)))
  if(length(test) > 1 | test[1] != TRUE) {
    stop("featureNames for are not consistent for all items of data.list")
  }    
  feature.number <- nrow(data.list[[1]])

  ## Check unique ids are really unique
  uids <- as.character( uids )
  if( any( duplicated( uids ) ) )  {
    stop("uids must be a vector of unique ids")
  }

  ## start processing
  if( data.classes == "ExpressionSet") {
    .process_arrays(
                    data.list = data.list,
                    uids = uids,
                    sample.annotation = sample.annotation,
                    platform.annotation = platform.annotation,
                    control_perturb_col = control_perturb_col,
                    control = control,
                    perturb = perturb,
                    limma = limma,
                    limma.index = limma.index,
                    big.matrix = big.matrix,
                    center.z = center.z,
                    center.log_fc = center.log_fc
                    )
    
  } else if (data.classes == "CountDataSet") {
    if( is.element("DESeq", installed.packages()[,1])){
      require( DESeq )
    } else {
      stop("To analyze RNAseq count data, please install the Bioconductor package 'DESeq'.")
    }
    data.list <- lapply( data.list, function(x){
      if( ! identical(varLabels (x ), c("sizeFactor", "condition") )){
        newCountDataSet( counts( x), conditions=pData( x)[,control_perturb_col])
      } else {
        return( x )
      }
    })
    .process_counts(
                    data.list = data.list,
                    uids = uids,
                    sample.annotation = sample.annotation,
                    platform.annotation = platform.annotation,
                    control_perturb_col = "condition",
                    control = control,
                    perturb = perturb,
                    big.matrix = big.matrix,
                    center.z = center.z,
                    center.log_fc = center.log_fc
                    )
  }
}

.process_arrays <- function(data.list,
                            uids,
                            sample.annotation,
                            platform.annotation,
                            control_perturb_col,
                            control,
                            perturb,
                            limma,
                            limma.index,
                            big.matrix,
                            center.z,
                            center.log_fc
                            ) 
{
  stopifnot(center.z %in% c("none", "mean", "median", "peak") )
  stopifnot(center.log_fc %in% c("none", "mean", "median", "peak") )

  if ( limma == TRUE ) {
    res <- lapply( data.list, function( x ) 
                  try( pairwise_compare_limma( 
                                              x,
                                              control_perturb_col = control_perturb_col,
                                              control=control, perturb=perturb ) ,
                      silent=TRUE ) 
                  ) 

  } else { ## standard t-tests
    res <- lapply( data.list, function( x ) 
                  try( pairwise_compare( 
                                        x,
                                        control_perturb_col = control_perturb_col,
                                        control=control, perturb=perturb ) ,
                      silent=TRUE ) 
                  ) 
  }
  
  ## remove failed instances
  bad.instances = sapply( res, is, "try-error" ) 
  
  if( any( bad.instances ) ) {
    warning( sprintf( "The following instances could not be analyzed: %s", paste( names( res[ bad.instances ] ) , collapse=", " ) ) ) 
    res = res[ ! bad.instances ]
    uids = uids[ ! bad.instances ]
    if( ! is.null( sample.annotation ) ) {
      sample.annotation <- sample.annotation[ which( ! bad.instances ) ,]
    }
  }
  
  ## collect results from all instances
  AveExpr = sapply( res, "[[", "exprs" ) 
  z = sapply( res, "[[", "z" ) 
  p = sapply( res, "[[", "p" ) 
  log_fc = sapply( res, "[[", "log_fc" ) 
  
  dimnames( AveExpr ) <- dimnames( z ) <- dimnames( p ) <- dimnames( log_fc ) <- list( row.names( res[[ 1 ]] ) , uids ) 

  ## create in-memory NChannelSet
  assay.data <- assayDataNew( exprs=AveExpr, z=z, p=p, log_fc=log_fc ) 

  fdata = data.frame( 
    probeid=row.names( res[[1]] ) ,
    row.names=row.names( res[[1]] ) 
    ) 
  
  pdata = data.frame( 
    UID = uids,
    row.names = uids
    ) 
  
  if ( !is.null( sample.annotation ) ) {
    pdata = cbind( pdata, sample.annotation ) 
  }

  vdata <- data.frame( 
                      labelDescription=colnames( pdata ) , 
                      channel=factor( rep( "_ALL_", ncol( pdata ) ) , levels=c( "_ALL_", colnames( pdata ) ) ) ) 
  
  ncs <- new( "NChannelSet",
             assayData = assay.data,
             featureData = new( "AnnotatedDataFrame", data = fdata ) ,
             phenoData = new( "AnnotatedDataFrame", data = pdata, varMetadata = vdata ) ,
             annotation = platform.annotation
             ) 
  
  ## center z-scores and log_fc channels
  ncs <- center_eSet( ncs, "z", center=center.z)
  ncs <- center_eSet( ncs, "log_fc", center=center.log_fc)
  
  ## create NChannelSet on disk
  if( ! is.null( big.matrix ) ) { ## big.matrix = path to BigMatrix file on disk
    eSetOnDisk( ncs, out.file=big.matrix ) 
  }

  return(ncs)
}

.process_counts <- function(data.list,
                            uids,
                            sample.annotation,
                            platform.annotation,
                            control_perturb_col,
                            control,
                            perturb,
                            big.matrix,  
                            center.z,
                            center.log_fc
                            )
{
  if( is.element("DESeq", installed.packages()[,1])){
    require( DESeq )
  } else {
    stop("To analyze RNAseq count data, please install the Bioconductor package 'DESeq'.")
  }
  
  stopifnot(center.z %in% c("none", "mean", "median", "peak") )
  stopifnot(center.log_fc %in% c("none", "mean", "median", "peak") )

  vst <- .vst_transform( data.list ) ## variance-stabilizing transformation
  res <- lapply( data.list, function( x )
                try(
                    pairwise_DESeq( x,
                                   vst,
                                   control_perturb_col = control_perturb_col,
                                   control=control, perturb=perturb,
                                   try.hard=FALSE),
                    silent=TRUE)
                )      
  
  ## remove failed instances
  bad.instances = sapply(res, is, "try-error")
  
  if ( any( bad.instances ) ) {
    warning( sprintf( "The following instances could not be analyzed: %s",
                     paste( names( res[ bad.instances ] ),
                           collapse=", ")
                     )
            )

    res  = res [ ! bad.instances ]
    uids = uids[ ! bad.instances ]
    if( ! is.null( sample.annotation ) ) {
      sample.annotation <- sample.annotation[ which( ! bad.instances ), ]
    }
  }

  ## collect results from all instances
  AveExpr = sapply( res, "[[", "exprs" )
  z = sapply( res, "[[", "z" )
  p = sapply( res, "[[", "p" )
  log_fc = sapply( res, "[[", "log_fc" )
  mod_fc = sapply( res, "[[", "mod_fc" )

  dimnames( AveExpr ) <- dimnames( z ) <- dimnames( p ) <- dimnames( mod_fc ) <- dimnames( log_fc ) <- list( row.names( res[[1]] ), uids)

  
  ## create in-memory NChannelSet
  assay.data <- assayDataNew( exprs=AveExpr, z=z, p=p, log_fc=log_fc, mod_fc=mod_fc)

  fdata = data.frame(
    probeid   = row.names( res[[ 1 ]] ),
    row.names = row.names( res[[ 1 ]] )
    )
  
  pdata = data.frame(
    UID = uids,
    row.names = uids
    )
  
  if ( !is.null( sample.annotation ) ) {
    pdata = cbind( pdata, sample.annotation)
  }
  
  vdata <- data.frame(
                      labelDescription = colnames( pdata ), 
                      channel=factor( rep( "_ALL_", ncol( pdata )),
                        levels=c( "_ALL_", colnames( pdata )
                          )
                        )
                      )
  
  ## in-memory NChannelSet
  ncs <- new( "NChannelSet",
             assayData   = assay.data,
             featureData = new( "AnnotatedDataFrame", data = fdata ),
             phenoData   = new( "AnnotatedDataFrame", data = pdata, varMetadata=vdata),
             annotation  = platform.annotation
             )

  ## center z-scores and log_fc channels
  ncs <- center_eSet( ncs, "z", center=center.z)
  ncs <- center_eSet( ncs, "log_fc", center=center.log_fc)
  ncs <- center_eSet( ncs, "mod_fc", center=center.log_fc)
  
  ## create NChannelSet on disk
  if(  ! is.null( big.matrix) ) { ## big.matrix = path to BigMatrix file on disk
    eSetOnDisk(  ncs, out.file=big.matrix )
  }
  return( ncs )
}

.vst_transform <- function( cds.list ) {
  if( is.element("DESeq", installed.packages()[,1])){
    require( DESeq )
  } else {
    stop("To perform variance stabilizing normalization, please install the Bioconductor package 'DESeq'.")
  }
  
  ## collect all counts in a single data matrix
  list.of.counts <- lapply( cds.list, counts )
  counts.matrix <- do.call( cbind, list.of.counts ) ## contains duplicates of control samples
  counts.matrix <- counts.matrix[, ! duplicated( colnames( counts.matrix ) ) ]
  
  ## estimate SizeFactors and dispersions
  cds <- DESeq::newCountDataSet( counts.matrix, colnames( counts.matrix ) )
  cds <- DESeq::estimateSizeFactors( cds )
  try( cds.fit <- DESeq::estimateDispersions( cds, method = "blind", sharingMode =  "fit-only") )

  if( is(cds.fit,"try-error")) {
    warning(print(".vst_transform: Parametric dispersion estimate for variance stabilizing tranformation failed. Trying local fit instead."))
    try( cds.fit <- DESeq::estimateDispersions( cds, method = "blind", sharingMode =  "fit-only", fitType="local") )
    
    if( is(cds.fit,"try-error")) {
      stop(".vst_transform: Dispersion estimate failed using parametric or local fit.")
    }
  }
  
  ## variance-stabilizing transformation
  vst <- DESeq::getVarianceStabilizedData( cds.fit )
}

mapNmerge <- function(eset, translation.fun = NULL, get="ENTREZID", verbose=FALSE, summary.fun=function(x)mean(x, na.rm=TRUE)) {

  if( !is.null( translation.fun )) {
    featureNames(eset) <- translation.fun(featureNames(eset))
    annotation(eset) <- get
    gene.ids <- featureNames(eset)

  } else {
    if ( length(annotation(eset)) == 0 ) {
      stop("eSet does not contain annotation information and 'translation.fun' has not been specified.")
    }

    ## load annotation package an retrieve ENTREIDs
    gene.ids <- unlist(
                       lookUp(
                              featureNames(eset),
                              annotation(eset),
                              what=get,
                              load=TRUE
                              )
                       )

    if ( all(is.na(gene.ids))) {
      stop("None of the feature ids could be mapped to EntrezIds.")
    }
    
    if( verbose == TRUE) {
      message(
            paste(
                  length(
                         gene.ids[!is.na(gene.ids)]
                         ),
                  "of",
                  dim(eset)[1],
                  "feature ids were successfully mapped to",
                  length(unique(na.omit(gene.ids))),
                  paste(get,"s", sep="")
                  )
            )
    }
  }

  if( any( duplicated( gene.ids ) ) & is.null( summary.fun )) {
    stop("Multiple features mapped to the same identifier, but no 'summary.fun' was specified.")

  } else if (is.null( summary.fun )) {
    new.eset <- eset
    
  } else {

    ## create new output eSet
    new.eset <- eset[,FALSE] ## metaData only
    annotation(new.eset) <- get

    ## apply 'summary.fun' to multiple probes mapped to the same EntrezId
    for (element in assayDataElementNames(eset) ) {
      element.sum <- matrix(
                            apply(assayDataElement(eset, element),
                                  2,
                                  function(column) {
                                    ave(column,
                                        gene.ids,
                                        FUN=summary.fun)
                                  }
                                  ),
                            dimnames=list(
                              gene.ids,
                              sampleNames(eset)),
                            ncol=dim(eset)[2]
                            )

      if( class(eset) == "CountDataSet") { ## integers allowed only
        element.sum <- round(element.sum,0)
      }

      element.sum <- as.matrix(element.sum[! (is.na(row.names(element.sum)) | duplicated(row.names(element.sum))) ,])
      assayDataElement(new.eset, element) <- element.sum
    }
    featureData(new.eset) <- as(data.frame(row.names=row.names(element.sum)), "AnnotatedDataFrame")
    phenoData(new.eset) <- phenoData(eset)
  }
  return( new.eset )
}

memorize <- function (object,
                      names=assayDataElementNames( object ),
                      ...) {
    if (any(duplicated(names)))
      stop("Channel 'names' must be unique")

    channelNames <- assayDataElementNames(object)
    badNames <- !names %in% channelNames

    if (any(badNames))
      stop("Channel 'names' must be channels")

    dropChannels <- channelNames[!channelNames %in% names]

    assayData <- .assayDataSubsetElements(assayData(object), names)
    metadata <- varMetadata(object)[["channel"]]
    okMetadata <- !metadata %in% dropChannels
    phenoData <- phenoData(object)[, okMetadata]
    varMetadata(phenoData)[["channel"]] <- factor(metadata[okMetadata],
                                                  levels = unique(c(names, "_ALL_")))
    initialize(object, assayData = assayData, phenoData = phenoData,
               featureData = featureData(object), experimentData = experimentData(object),
               annotation = annotation(object), protocolData = protocolData(object)
               ,...)
  }

.assayDataSubsetElements <- function (object, elts) {
    if (any(duplicated(elts)))
      stop("'AssayData' element names must be unique")
    
    storageMode <- if (is(object, "list")) {
      "list"
    } else if (environmentIsLocked(object)) {
      "lockedEnvironment" 
    } else {
      "environment"
    }

    names <- if (storageMode(object) == "list") {
      names(object)
    } else {
      ls(object)
    }
    
    if (!all(elts %in% names))
      stop("'AssayData' missing elements: '", paste(elts[!elts %in%
                                                         names], collapse = "', '", sep = ""), "'")
    switch(storageMode,
           lockedEnvironment = {
             assayData <- new.env(parent = emptyenv())
             for (nm in elts) {
               assayData[[nm]] <- as( object[[nm]], "matrix" )
             }
             lockEnvironment(assayData, bindings = TRUE)
             assayData
           },           
           environment = {
             assayData <- new.env(parent = emptyenv())
             for (nm in elts) assayData[[nm]] <- as( object[[nm]], "matrix" )
             assayData
           }, list = {
             object[elts]
           })
  }

pickChannels <- function (object, names, ...) {
    if (any(duplicated(names)))
      stop("Channel 'names' must be unique")

    channelNames <- assayDataElementNames(object)
    badNames <- !names %in% channelNames

    if (any(badNames))
      stop("Channel 'names' must be channels")

    dropChannels <- channelNames[!channelNames %in% names]

    assayData <- .assayDataSubsetElements(assayData(object), names)
    metadata <- varMetadata(object)[["channel"]]
    okMetadata <- !metadata %in% dropChannels
    phenoData <- phenoData(object)[, okMetadata]
    varMetadata(phenoData)[["channel"]] <- factor(metadata[okMetadata],
                                                  levels = unique(c(names, "_ALL_")))
    initialize(object, assayData = assayData, phenoData = phenoData,
               featureData = featureData(object), experimentData = experimentData(object),
               annotation = annotation(object), protocolData = protocolData(object)
               ,...)
  }

.assayDataSubsetElements <- function (object, elts) {
    if (any(duplicated(elts)))
      stop("'AssayData' element names must be unique")
    
    storageMode <- if (is(object, "list")) {
      "list"
    } else if (environmentIsLocked(object)) {
      "lockedEnvironment" 
    } else {
      "environment"
    }

    names <- if (storageMode(object) == "list") {
      names(object)
    } else {
      ls(object)
    }
    
    if (!all(elts %in% names))
      stop("'AssayData' missing elements: '", paste(elts[!elts %in%
                                                         names], collapse = "', '", sep = ""), "'")
    switch(storageMode,
           lockedEnvironment = {
             assayData <- new.env(parent = emptyenv())
             for (nm in elts) {
               assayData[[nm]] <- as( object[[nm]], "matrix" )
             }
             lockEnvironment(assayData, bindings = TRUE)
             assayData
           },           
           environment = {
             assayData <- new.env(parent = emptyenv())
             for (nm in elts) assayData[[nm]] <- as( object[[nm]], "matrix" )
             assayData
           }, list = {
             object[elts]
           })
  }


eset_instances <- function(instance.matrix, 
                           eset, 
                           control_perturb_col="cmap", 
                           control="control", 
                           perturb="perturbation") {
  stopifnot(row.names( instance.matrix) %in% sampleNames( eset ))
  stopifnot(all( instance.matrix %in% c(-1,0,1)))
  
  apply( instance.matrix, 2, function( x ){
    control.samples <- sampleNames( eset )[x == -1]
    perturbation.samples <- sampleNames( eset )[x == 1]
    instance.eset <- eset[, c(control.samples, perturbation.samples) ]
    pData( instance.eset )[,control_perturb_col] <- c( rep( control, 
                                                           length(control.samples )
                                                           ), 
                                                      rep( perturb, 
                                                          length(perturbation.samples )
                                                          )
                                                      )
    return( instance.eset )
  })
}

splitPerturbations <- function( eset, 
                                control="none",
                                controlled.factors="none", ## all, none or vector
                                factor.of.interest="Compound",
                                ignore.factors=NULL,
                                cmap.column="cmap",
                                prefix=NULL){
  if( inherits( eset, "CountDataSet")){
    ignore.factors <- c(ignore.factors, "sizeFactor")
  }
  
  ## get annotation information from the phenoData slot
  pd <- pData( eset )
  
  ## remove columns matching the 'ignore.factors' parameter
  if( !is.null(ignore.factors)){
    ignore.factors <- sapply(ignore.factors, function( x ){
      if( is.null(prefix)){
        m <- grep( x, colnames(pd), value=TRUE )        
      } else {
        m <- grep( x, grep( prefix, colnames(pd), value=TRUE), value=TRUE )
      }
      if( length(m) == 0){
        return(NA)
      } else {
        return( m )
      }
    })
    
    if( all(is.na(ignore.factors))){
      stop("None of the 'ignore.factors' could be found in the phenoData.")
    } else {
      message( sprintf("The following factor(s) will be ignored: %s",
                       paste(ignore.factors, collapse=", ")))
      pd <- pd[,!colnames(pd) %in% ignore.factors]
    }
  }
  
  ## identify experimental factors by the "Factor" prefix
  if( is.null( prefix )){
    factor.all <- colnames(pd)
  } else {
    factor.all <- grep( prefix, colnames(pd), value=TRUE)
  }
  
  ## match / extract the user-specified factors of interest
  factor.of.interest <- grep( factor.of.interest, factor.all, value=TRUE)
  if( length( factor.of.interest ) == 0){
    stop( "Your factor of interest could not be found in the dataset.")
  } else if( length( factor.of.interest) > 1) {
    stop( "The 'factor.of.interest' parameter matches multiple annotation columns.")
  }
  message(sprintf("Using %s as factor of interest.", factor.of.interest))
  
  other.factors <- setdiff( factor.all, factor.of.interest)
  perturbations <- unique( setdiff( pd[,factor.of.interest], control))
  
  if( length( perturbations) == 0){
    stop( "No perturbations could be found in the 'factor.of.interest' column.")
  }
  
  if( !identical( controlled.factors,"all") & !identical( controlled.factors,"none")){
    controlled.factors <- sapply(controlled.factors, function( x) { 
      grep( x, factor.all,value=TRUE)
    })
  } else if( controlled.factors == "all"){
    controlled.factors <- other.factors
  } else if( controlled.factors == "none"){
    controlled.factors <- controlled.factors
  }
  
  ## identify control and perturbation instances
  if( length( other.factors) == 0){
    control.instances <- pd[pd[, factor.of.interest] == control, factor.of.interest, drop=FALSE ]
  } else {
    control.instances <- pd[pd[, factor.of.interest] == control, other.factors,drop=FALSE ]
  }
  
  if( nrow( control.instances) == 0){
    stop( sprintf("No control instances could be found matching %s in the %s column.", control, factor.of.interest))
  }
  
  perturb.instances <- pd[pd[, factor.of.interest] != control, factor.all, drop=FALSE ]
  
  ## iterate over perturbations and compile single-factor experiments
  all.instances <- lapply( perturbations, function( x ) {
    ## unique experimental conditions
    if( length(other.factors) == 0){
      perturb.unique <- unique(  perturb.instances[ which( perturb.instances[,factor.of.interest] == x),factor.of.interest, drop=FALSE])
      perturb.all <- perturb.instances[which( perturb.instances[,factor.of.interest] == x), factor.of.interest, drop=FALSE]
    } else {
    perturb.unique <- unique(  perturb.instances[ which( perturb.instances[,factor.of.interest] == x),other.factors, drop=FALSE])
    perturb.all <- perturb.instances[which( perturb.instances[,factor.of.interest] == x), other.factors, drop=FALSE]
    }
    
    ## assign samples to unique experimental conditions
    matched.samples <- lapply( 1:nrow( perturb.unique), 
                               function( n ){
                                 sapply( 1:nrow( perturb.all),
                                         function( m ){
                                           all( perturb.all[m,] == perturb.unique[n,])
                                         })
                               })
    
    matched.samples <- lapply( matched.samples, function(x){ 
      row.names(perturb.all)[x]
    })
    
    if( any( sapply( matched.samples, length) != 1)){
      matched.samples <- I(matched.samples)
    } else {
      matched.samples <- unlist(matched.samples)
    }
    
    ## identify suitable controls
    if( identical( controlled.factors,"none")){
      matched.controls <- lapply( 1:nrow( perturb.unique),function(n){
        row.names( control.instances )
      })
    } else {
      matched.controls <- lapply( 1:nrow( perturb.unique), 
                                  function( n ){
                                    sapply( 1:nrow( control.instances),
                                            function( m ){
                                              all( control.instances[m, controlled.factors] == perturb.unique[n,controlled.factors])
                                            })
                                  })
      matched.controls <- lapply( matched.controls, function(x){ 
        row.names(control.instances)[x]
      })
    }
    
    if( all( is.na(unlist(matched.controls)))){
      warning(sprintf("No control samples found for %s", x), call.=FALSE)
      return( NA)
    }
    
    if( any( sapply( matched.controls, length) != 1)){
      matched.controls <- I(matched.controls)
    } else {
      matched.controls <- unlist(matched.controls)
    }
    
    ## return results as a data.frame
    data.frame(perturbagen=x, 
               perturb.unique, 
               data.frame( matched.samples, stringsAsFactors=FALSE), 
               data.frame( matched.controls,stringsAsFactors=FALSE))
  })
  names( all.instances) <- perturbations
  
  instances <- unlist(
    lapply( all.instances[!is.na(all.instances)], function( x ){
      lapply( 1:nrow( x ), function( n) {
        selected.sampes <- c( x[n,"matched.samples"], 
                              x[n,"matched.controls"])
        selected.sampes <- unlist( selected.sampes)
      })
    }), recursive=FALSE)
  
  eset.list <- lapply( instances, function( x ){
    instance.eset <- eset[, x]
    pData(instance.eset)[,cmap.column] <- ifelse( pData( instance.eset)[,factor.of.interest] == control, "control", "perturbation")
    return( instance.eset)
  })

  return( eset.list)
}

annotate_eset_list <- function(eset.list, cmap.column="cmap", 
                               perturbation="perturbation") {
  common.varLabels <- Reduce('intersect', lapply( eset.list, varLabels))
  common.varLabels <- setdiff( common.varLabels, cmap.column )
  
  sample.anno <-
    t(sapply( eset.list, function( x ){
      perturb <- pData( x )
      perturb <- perturb[ perturb[,cmap.column] == perturbation, 
                          common.varLabels, drop=FALSE]
      y <- apply( perturb, 2, function(x) {
        paste( unique( x ), collapse=", ")
      })
      as.character( y )
    }))
  
  sample.anno <- data.frame( apply( sample.anno, 2, as.character))
  colnames(sample.anno)  <- common.varLabels
  return( sample.anno )
}


mergeCMAPs <- function(x, y){
  ## basic checks
  if( ! inherits(x, "eSet")){
    stop( "Object 'x' is not an eSet.")
  }
  
  if( ! inherits(y, "eSet")){
    stop( "Object 'y' is not an eSet.")
  }
     
  if( ! class(x ) == class( y )){
    stop( "Objects 'x' and 'y' are not of the same class.")
  }

  if( annotation( x ) != annotation( y )){
    stop( "Objects 'x' and 'y' have different 'annotation' slots.")
  }

  if( any( assayDataElementNames( x ) != assayDataElementNames( y))){
    stop( "Objects 'x' and 'y' have different AssayDataElementNames.")
  }

  if( length( intersect( sampleNames( x ), sampleNames( y ))) != 0 ){
    stop( sprintf( "Objects 'x' and 'y' share %s sampleNames. sampleNames must be unique.", length( intersect( sampleNames( x ), sampleNames( y )))))
  }

  if(  any( varLabels( x ) !=  varLabels( y ))){
    stop( "Objects 'x' and 'y' have different pData columns.")
  }
  
  common.features <- intersect( featureNames( x ), featureNames( y ))
  message(sprintf("eSets 'x' and 'y' share %s common features.", length(common.features)))
  
  ## create empty eSet for output
  merged.eset <- new( class( x ) ) 
  annotation(merged.eset) <- annotation( x )
  
  ## merge data for each channel separately
  for (element in assayDataElementNames( x ) ) {
    merged.channel <- merge(
                            assayDataElement( x, element)[,],
                            assayDataElement( y, element)[,],
                            by.x=0, by.y=0, all=TRUE
                            )
    row.names( merged.channel) <- merged.channel$Row.names
    merged.channel$Row.names <- NULL
    assayDataElement( merged.eset, element) <- as.matrix( merged.channel )
  }

  ## merge pData tables => requires the same column names in each eSet !
  pData(merged.eset) <- rbind( pData( x ), pData( y ))
  featureData(merged.eset) <- AnnotatedDataFrame(
                                                 data.frame(
                                                            probeid=row.names( merged.channel),
                                                            row.names=row.names(merged.channel)
                                                            )
                                                 )
  return( merged.eset)
}

center.function <- function(x, type) {
  switch(type,
         mean = mean(x, na.rm=TRUE),
         median = median(x, na.rm=TRUE),
         peak = {
           d <- density( x, adjust=2, na.rm=TRUE  )
           d$x[ d$y == max(d$y)][1]
         })
}

center_eSet <- function( eset,
                                channel,
                                center="peak"){
  
  stopifnot( inherits( eset, "eSet"))
  stopifnot( center %in% c("none", "mean", "median", "peak") )
  stopifnot( channel %in% assayDataElementNames( eset ))  
  
  center.function <- function(x, type) {
    switch(type,
           mean = mean(x, na.rm=TRUE),
           median = median(x, na.rm=TRUE),
           peak = {
             d <- density( x, adjust=2, na.rm=TRUE  )
             d$x[ d$y == max(d$y)][1]
           })
  }
  
  dat <- assayDataElement( eset, channel)
  
  if( center != "none"){
    dat.shift <- apply( dat, 2, function(x){
      y <- try( center.function( x, center), silent=TRUE)
      if( inherits( y, "try-error")){
        return( NA )
      } else {
        return( y )
      }
    })
    dat <- sweep( dat, 2, dat.shift)
    assayDataElement( eset, channel) <- dat
  } else {
    dat.shift <- rep(NA, ncol( eset ))
  }
  pData( eset )[,paste(channel, "shift", sep=".")] <- dat.shift
  varMetadata( eset )[paste(channel, "shift", sep="."),
                      "labelDescription"] <- sprintf("center of the uncorrected %s distribution", channel)
  return( eset)
}


reactome2cmap <- function(species, annotation.package){ 
  if( is.element("reactome.db", installed.packages()[,1])){
    require( "reactome.db",character.only = TRUE )
  } else {
    stop("To run this function, please install the Bioconductor package 'reactome.db'.")
  }

  if( is.element(annotation.package, installed.packages()[,1])){
    require( annotation.package,character.only = TRUE )
  } else {
    stop(sprintf("The specified annotation package %s is not installed on this system.", annotation.package))
  }
  
  pathways <- as.list(reactomePATHID2EXTID)
  
  ## retrieve names
  pathway.names <- unlist(AnnotationDbi::mget(names(pathways), reactomePATHID2NAME))
  pathway.names <- pathway.names[ match(names( pathways),
                                        names( pathway.names )) ]
  
  ## remove categories with duplicated or missing names
  filtered.names <- duplicated( names( pathway.names)) | is.na(pathway.names)
  pathways <- pathways[ ! filtered.names ]
  pathway.names <- pathway.names[ ! filtered.names]
  
  selected.species <- grepl( paste("^", species, sep=""), pathway.names)
  
  pheno.data <- as(
    data.frame(name=pathway.names[ selected.species ],
               row.names=names(pathways[ selected.species ])
    ),
    "AnnotatedDataFrame")
  
  i.matrix <- Matrix::t( incidence( pathways[ selected.species ] ) )
  reactome <- CMAPCollection( i.matrix,
                                 phenoData=pheno.data,
                                 annotation=annotation.package,
                                 signed=rep( FALSE, ncol(i.matrix)) )
  return( reactome)
}

KEGG2cmap <- function( species, annotation.package ){
  if( is.element("KEGG.db", installed.packages()[,1])){
    require( "KEGG.db",character.only = TRUE )
  } else {
    stop("To run this function, please install the Bioconductor package 'KEGG.db'.")
  }

  if( is.element(annotation.package, installed.packages()[,1])){
    require( annotation.package,character.only = TRUE )
  } else {
    stop(sprintf("The specified annotation package %s is not installed on this system.", annotation.package))
  }
  
  ## retrieve entrez ids of pathw ay members
  pathways <- as.list(KEGGPATHID2EXTID)
  
  ## retrieve names
  pathway.names <- unlist(AnnotationDbi::mget(sub("^...", "",names(pathways)), KEGGPATHID2NAME))
  
  ## species-specific CMAPCollections
  selected.species <- grepl( paste("^", species, sep=""), names( pathways ))
  pheno.data <- as(
    data.frame(name=pathway.names[ selected.species ],
               row.names=names(pathways[ selected.species ])
    ),
    "AnnotatedDataFrame")
  
  i.matrix <- Matrix::t( incidence( pathways[ selected.species ] ) )
  kegg <- CMAPCollection( i.matrix,
                                 phenoData=pheno.data,
                                 annotation=annotation.package,
                                 signed=rep( FALSE, ncol(i.matrix)) )
  return( kegg )
}

wiki2cmap <- function( species, annotation.package ){  
  if( is.element(annotation.package, installed.packages()[,1])){
    require( annotation.package,character.only = TRUE )
  } else {
    stop(sprintf("The specified annotation package %s is not installed on this system.", annotation.package))
  }
  
  ## Download wikipathways in plain text format
  url <- sprintf( "http://www.wikipathways.org//wpi/batchDownload.php?species=%s&fileType=txt&tag=Curation:AnalysisCollection", URLencode( species))
  download.file(url, file.path(tempdir(),"wiki.zip" ))
  unzip( file.path(tempdir(),"wiki.zip" ), 
         exdir=file.path( tempdir(), "wiki")
  )
  pathways <- lapply( list.files( file.path( tempdir(), "wiki"), 
                                  full.names=TRUE, pattern=".txt"), 
                      read.delim)
  names( pathways ) <- sub(".txt", "",  list.files( file.path( tempdir(), "wiki"), 
                                                    pattern=".txt" ))
  
  pathways <- lapply( pathways, function(x) {
    entrez <- as.character( x[which(x$Database == "Entrez Gene"), "Identifier"])
    ensembl.ids <- grep( "Ensembl", x$Database)
    if( length( ensembl.ids ) > 0){
      ensembl.ids <- as.character( x[ensembl.ids, "Identifier"])
      ## construct environment name    
      ## query environment for Entrez Ids
      s2e <- get( paste( sub( ".db$", "", annotation.package), "ENSEMBL2EG", sep=""))
      entrez.ids <- AnnotationDbi::mget(ensembl.ids, s2e, ifnotfound=NA)
      entrez.ids <- sapply( entrez.ids, "[[", 1) ## for multi-matches, use first Entrez Id
    } else {
      entrez.ids <- NULL
    }
    return( na.omit(unique(c(entrez, entrez.ids))))
  })

  wiki.sets <- as( pathways, "CMAPCollection")
  return( wiki.sets )
}

go2cmap <- function( annotation.package="org.Hs.eg.db", ontology="BP", evidence=NULL){
  GOBPANCESTOR <- GOMFANCESTOR <- GOCCANCESTOR <- NULL
  
  if( is.element("GO.db", installed.packages()[,1])){
    require( "GO.db",character.only = TRUE )
    message("Successfully loaded GO.db package.")
  } else {
    stop("To use this function, please install the Bioconductor package 'GO.db'.")
  }
  
  if( is.element(annotation, installed.packages()[,1])){
    require( annotation,character.only = TRUE )
    message(paste("Successfully loaded",annotation, "package."))
  } else {
    stop(sprintf("The specified annotation page % was not found.", annotation))
  }
  
  consolidate.gomap = function(gomap) {
    colnames(gomap)[2] = "go1"
    gomap = gomap[gomap$go1!="",]
    
    ## build GO ancestors
    cat("consolidate ancestor terms...\n")
    allgo = unique(gomap$go1)
    ac = c(as.list(GOBPANCESTOR), as.list(GOMFANCESTOR), as.list(GOCCANCESTOR))
    ac = ac[allgo]
    ac = ac[sapply(ac, length)>0]
    acdf = data.frame(go1=rep(names(ac), sapply(ac, length)), go=unlist(ac), stringsAsFactors=FALSE)
    acdf = rbind(acdf, data.frame(go1=allgo, go=allgo, stringsAsFactors=FALSE))
    acdf = acdf[acdf$go!="all",]
    
    ## consolidate gomap
    gomap = merge(gomap, acdf, by="go1", all.x=TRUE)
    gomap$go1= NULL
    unique(gomap)
  }
  
  ## build gene/terms association maps
  build.gomaps <- function(anno, ontology=NULL, evidence=NULL) {
    ggo <- as.list(get(anno))
    
    if( !is.null( evidence )) {
      ggo <- mclapply( ggo, function( gene ) {
        if( inherits( gene, "list" ) ) {
          gene[ sapply(gene, function( category ) { category[["Evidence"]] %in% evidence }) ]
        } else {
          NA
        }
      })
    }
    
    
    if( !is.null( ontology )) {
      ggo <- mclapply( ggo, function( gene ) {
        if( inherits( gene, "list" ) ) {
          gene[ sapply(gene, function( category ) { category[["Ontology"]] %in% ontology }) ]
        } else {
          NA
        }
      })
    }
    
    ggo <- ggo[ ! sapply( ggo, function(x) inherits(x, "logical")) ]
    ggo <- mclapply(ggo, function(g) names(g))
    gomap <- data.frame(gene=rep(names(ggo), sapply(ggo, length)), go=unlist(ggo), stringsAsFactors=FALSE)
    gomap <- consolidate.gomap(gomap)
    return(gomap)
  }
  
  anno <- paste( sub(".db$", "", annotation),"GO", sep="")
  go  <-build.gomaps( anno=anno, ontology=ontology )
  gene.sets <- split(go$gene, go$go)
  
  ## retrieve go-term names
  go.names <- as.data.frame(Term(names( gene.sets )))
  colnames(go.names) <- "Name"
  
  ## create CMAP
  go.cmap <- CMAPCollection( Matrix::t(incidence(gene.sets)), phenoData=as(go.names, "AnnotatedDataFrame"))
  signed(go.cmap) <- rep(FALSE, ncol( go.cmap))
  
  ## add additional annotations
  experimentData( go.cmap)@title <- paste("GO", ontology, "ontology")
  experimentData( go.cmap)@abstract <- sprintf( "%s categories from the GO %s ontology", ncol( go.cmap), ontology)
  return( go.cmap )
}
