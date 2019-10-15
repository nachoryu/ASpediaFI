setOldClass("igraph")
setClass("ASpediaFI",
         representation = representation(samples = "data.frame",
                                         events = "list",
                                         psi = "SummarizedExperiment",
                                         gtf = "GRanges",
                                         network = "igraph",
                                         gene.table = "data.frame",
                                         as.table = "data.frame",
                                         pathway.table = "data.frame"))

#Constructor
ASpediaFI <- function(sample.names, bam.files, conditions){

    sample.info <- data.frame(name = sample.names, path = bam.files,
                          condition = conditions,
                          stringsAsFactors = FALSE)
    new("ASpediaFI", samples = sample.info, events = list(),
        psi = SummarizedExperiment(), gtf = GRanges(),
        network = make_empty_graph(),
        gene.table = data.frame(), as.table = data.frame(),
        pathway.table = data.frame())

}

#Accessors
setGeneric("samples", function(object) {
    standardGeneric("samples")
})
setGeneric("samples<-", function(object, value) {
    standardGeneric("samples<-")
})
setMethod("samples", "ASpediaFI", function(object){
    return(object@samples)
})
setReplaceMethod("samples", "ASpediaFI", function(object, value) {

    #Make changes to object
    outFI <- object
    outFI@samples <- value

    return(outFI)
})

#Accessor for AS event annotations
setGeneric("events", function(object,...){
    standardGeneric("events")
})
setGeneric("events<-", function(object, value) {
    standardGeneric("events<-")
})
setMethod("events", "ASpediaFI", function(object){
    return(object@events)
})
setReplaceMethod("events", "ASpediaFI", function(object, value) {

    #Make changes to object
    outFI <- object
    outFI@events <- value

    return(outFI)
})

#Accessor for PSI values
setGeneric("psi", function(object,...){
    standardGeneric("psi")
})
setGeneric("psi<-", function(object, value) {
    standardGeneric("psi<-")
})
setMethod("psi", "ASpediaFI", function(object){
    return(object@psi)
})
setReplaceMethod("psi", "ASpediaFI", function(object, value) {

    #Make changes to object
    outFI <- object
    outFI@psi <- value

    return(outFI)
})

#Accessor for DRaWR network
setGeneric("network", function(object,...){
    standardGeneric("network")
})
setMethod("network", "ASpediaFI", function(object){
    return(object@network)
})

#Accessor for gene table
setGeneric("gene.table", function(object,...){
    standardGeneric("gene.table")
})
setMethod("gene.table", "ASpediaFI", function(object){
    return(object@gene.table)
})

#Accessor for as.table
setGeneric("as.table", function(object,...){
    standardGeneric("as.table")
})
setMethod("as.table", "ASpediaFI", function(object){
    return(object@as.table)
})

#Accessor for pathway.table
setGeneric("pathway.table", function(object,...){
    standardGeneric("pathway.table")
})
setMethod("pathway.table", "ASpediaFI", function(object){
    return(object@pathway.table)
})




