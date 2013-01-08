setGeneric("CMAPCollection",
           function(assayData,
                    phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                    featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                    annotation=character(),
                    protocolData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                    ...)
           standardGeneric("CMAPCollection"),
           signature="assayData")


setGeneric("induceCMAPCollection",
           def = function( eset, ... ) standardGeneric( "induceCMAPCollection" )
           )


setGeneric("CMAPResults", function(data, varMetadata, ...)
           standardGeneric("CMAPResults")
           )

setGeneric("setSizes",
           def = function( object) standardGeneric( "setSizes" )
)

setGeneric("members",
           def = function( object) standardGeneric( "members" )
           )

setGeneric(
           "signed",
           def = function( object) standardGeneric( "signed" )
           )

setGeneric("signed<-",
           def=function(x, value) standardGeneric("signed<-")
           )


setGeneric("mergeCollections",
           def = function( x, y) standardGeneric( "mergeCollections" )
           )

setGeneric("cmapTable",
           def=function(object, ...) standardGeneric("cmapTable")
           )

setGeneric("set",
           def=function(object) standardGeneric("set")
           )

setGeneric("set<-",
           def=function(x, value) standardGeneric("set<-")
           )

setGeneric("trend",
           def=function(object) standardGeneric("trend")
           )

setGeneric("trend<-",
           def=function(x, value) standardGeneric("trend<-")
           )

setGeneric("pval",
           def=function(object) standardGeneric("pval")
           )

setGeneric("pval<-",
           def=function(x, value) standardGeneric("pval<-")
           )

setGeneric("padj",
           def=function(object) standardGeneric("padj")
           )

setGeneric("padj<-",
           def=function(x, value) standardGeneric("padj<-")
           )

setGeneric("zscores",
           def=function(object) standardGeneric("zscores")
           )

setGeneric("effect",
           def=function(object) standardGeneric("effect")
           )

setGeneric("effect<-",
           def=function(x, value) standardGeneric("effect<-")
           )

setGeneric("nSet",
           def=function(object) standardGeneric("nSet")
           )

setGeneric("nSet<-",
           def=function(x, value) standardGeneric("nSet<-")
           )

setGeneric("nFound",
           def=function(object) standardGeneric("nFound")
           )

setGeneric("nFound<-",
           def=function(x, value ) standardGeneric("nFound<-")
           )

setGeneric("docs",
           def=function(object) standardGeneric("docs")
           )

setGeneric("docs<-",
           def=function(x, value) standardGeneric("docs<-")
           )

setGeneric("errors",
           def=function(object) standardGeneric("errors")
           )

setGeneric("errors<-",
           def=function(x, value) standardGeneric("errors<-")
           )

setGeneric("labels",
           def=function(object) standardGeneric("labels")
           )

setGeneric("labels<-",
           def=function(x, value) standardGeneric("labels<-")
           )

setGeneric("geneScores",
           def=function(object) standardGeneric("geneScores")
           )

setGeneric("plot",
           def=function(x, y, ...) standardGeneric("plot")
           )

setGeneric("SignedGeneSet",
           def = function( type, ... ) standardGeneric( "SignedGeneSet" )
           )

setGeneric("upIds",
           def = function( object, ... ) standardGeneric( "upIds" )
           )

setGeneric("downIds",
           def = function( object, ... ) standardGeneric( "downIds" )
           )

setGeneric("geneSign",
           def = function( obj ) standardGeneric( "geneSign" )
           )

setGeneric("geneSign<-",
           def = function( object, value ) standardGeneric( "geneSign<-" )
           )

setGeneric("camera_score",
           def = function( experiment, sets, ... ) standardGeneric( "camera_score" )
           )

setGeneric("connectivity_score",
           def = function( experiment, query, ... ) standardGeneric( "connectivity_score" )
           )

setGeneric("featureScores",
           def = function( query, dat, ... ) standardGeneric( "featureScores" )
           )

setGeneric("fisher_score",
           def = function( query, sets, universe, ... ) standardGeneric( "fisher_score" )
           )

setGeneric("geneIndex",
           def = function( gene.sets, gene.ids, ...) standardGeneric( "geneIndex" )
           )

setGeneric("gsealm_score",
           def = function( query, set, ... ) standardGeneric( "gsealm_score" )
           )

setGeneric("gsealm_jg_score",
           def = function( query, sets, ... ) standardGeneric( "gsealm_jg_score" )
           )

setGeneric("minSetSize",
           def = function( sets, ... ) standardGeneric( "minSetSize" )
           )

setGeneric("mroast_score",
           def = function( experiment, sets, ... ) standardGeneric( "mroast_score" )
           )

setGeneric("romer_score",
           def = function( experiment, sets, ... ) standardGeneric( "romer_score" )
           )

setGeneric("wilcox_score",
           def = function( experiment, sets, ... ) standardGeneric( "wilcox_score" )
           )

