# Privates

setGeneric(".getProteins", 
           function(object) {
             standardGeneric(".getProteins")
           }
)

setGeneric(".getTaxId", 
           function(object,...) {
             standardGeneric(".getTaxId")
           }
)

setGeneric(".getTaxName", 
           function(object,...) {
             standardGeneric(".getTaxName")
           }
)

setGeneric(".getInteractions", 
           function(object,...) {
             standardGeneric(".getInteractions")
           }
)

setGeneric(".getReferenceInteractions", 
           function(object,...) {
             standardGeneric(".getReferenceInteractions")
           }
)

setGeneric(".addReferenceSpecies", 
           function(object,taxName,taxId,proteins,irefindex.file) {
             standardGeneric(".addReferenceSpecies")
           }
)


setGeneric(".removeReferenceSpecies", 
           function(object,species) {
             standardGeneric(".removeReferenceSpecies")
           }
)

setGeneric(".getAllTaxIds", 
           function(object) {
             standardGeneric(".getAllTaxIds")
           }
)

setGeneric(".getReferenceSpecies", 
           function(object,species,info=FALSE,detailed=FALSE) {
             standardGeneric(".getReferenceSpecies")
           }
)

setGeneric(".getSourceDBs", 
           function(object, irefindex) {
             standardGeneric(".getSourceDBs")
           }
)

setGeneric(".countSpecies", 
           function(object) {
             standardGeneric(".countSpecies")
           }
)

setGeneric(".getHomologs", 
           function(object,taxId,...) {
             standardGeneric(".getHomologs")
           }
)

setGeneric(".addHomologs", 
           function(object,taxId,homologs) {
             standardGeneric(".addHomologs")
           }
)

setGeneric(".removeHomologs", 
           function(object,species) {
             standardGeneric(".removeHomologs")
           }
)

# Publics

setGeneric("addReference", 
           function(path2ppi,taxName,taxId,proteins,irefindex, homologs) {
             standardGeneric("addReference")
           }
)

setGeneric("showReferences", 
           function(path2ppi,species=NA,returnValue=NA) {
             standardGeneric("showReferences")
           }
)

setGeneric("plot.Path2PPI", 
           function (x,type="ppi",multiple.edges=FALSE,scores=FALSE,
                     species.colors=c(),vertices.opacity=0.8,
                     use.identifiers=FALSE,protein.labels=NA,show.legend=TRUE,
                     vertices.coordinates=NA,return.coordinates=FALSE,
                     tkplot=FALSE,...) {
             standardGeneric("plot.Path2PPI")
           }
)

setGeneric("getHybridNetwork", 
           function(path2ppi, igraph=FALSE) {
             standardGeneric("getHybridNetwork")
           }
)

setGeneric("getPPI", 
           function(path2ppi, raw=FALSE, igraph=FALSE) {
             standardGeneric("getPPI")
           }
)

setGeneric("showInteraction", 
           function(path2ppi,interaction, mode="default", verbose=TRUE) {
             standardGeneric("showInteraction")
           }
)

setGeneric("predictPPI", 
           function(path2ppi,mode="both",h.thresh=1e-5,h.range=c(1e-100,1e-20),
                    i.thresh=0.7,consider.complexes=FALSE,max.complex.size=5,
                    decline.self.interaction.ref=FALSE,
                    decline.self.interaction.tar=TRUE,verbose=TRUE) {
             standardGeneric("predictPPI")
           }
)

setGeneric("removeReference", 
           function(path2ppi,species) {
             standardGeneric("removeReference")
           }
)
