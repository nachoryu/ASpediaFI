#' @rdname ASpediaFI-accessor
#' @export
setGeneric("samples", function(object, ...) {
    standardGeneric("samples")
})

#' @rdname ASpediaFI-accessor
#' @export
setMethod("samples", "ASpediaFI", function(object) {
    return(object@samples)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("samples<-", function(object, value) {
    standardGeneric("samples<-")
})

#' @rdname ASpediaFI-accessor
#' @export
setReplaceMethod("samples", "ASpediaFI", function(object, value) {
    outFI <- object
    outFI@samples <- value
    return(outFI)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("events", function(object, ...) {
    standardGeneric("events")
})

#' @rdname ASpediaFI-accessor
#' @export
setMethod("events", "ASpediaFI", function(object) {
    return(object@events)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("events<-", function(object, value) {
    standardGeneric("events<-")
})

#' @rdname ASpediaFI-accessor
#' @export
setReplaceMethod("events", "ASpediaFI", function(object, value) {
    outFI <- object
    outFI@events <- value
    return(outFI)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("psi", function(object, ...) {
    standardGeneric("psi")
})

#' @rdname ASpediaFI-accessor
#' @export
setMethod("psi", "ASpediaFI", function(object) {
    return(object@psi)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("psi<-", function(object, value) {
    standardGeneric("psi<-")
})

#' @rdname ASpediaFI-accessor
#' @export
setReplaceMethod("psi", "ASpediaFI", function(object, value) {
    outFI <- object
    outFI@psi <- value
    return(outFI)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("gtf", function(object, ...) {
    standardGeneric("gtf")
})

#' @rdname ASpediaFI-accessor
#' @export
setMethod("gtf", "ASpediaFI", function(object) {
    return(object@gtf)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("gtf<-", function(object, value) {
    standardGeneric("gtf<-")
})

#' @rdname ASpediaFI-accessor
#' @export
setReplaceMethod("gtf", "ASpediaFI", function(object, value) {
    outFI <- object
    outFI@gtf <- value
    return(outFI)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("network", function(object, ...) {
    standardGeneric("network")
})

#' @rdname ASpediaFI-accessor
#' @export
setMethod("network", "ASpediaFI", function(object) {
    return(object@network)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("network<-", function(object, value) {
    standardGeneric("network<-")
})

#' @rdname ASpediaFI-accessor
#' @export
setReplaceMethod("network", "ASpediaFI", function(object, value) {
    outFI <- object
    outFI@network <- value
    return(outFI)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("gene.table", function(object, ...) {
    standardGeneric("gene.table")
})

#' @rdname ASpediaFI-accessor
#' @export
setMethod("gene.table", "ASpediaFI", function(object) {
    return(object@gene.table)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("gene.table<-", function(object, value) {
    standardGeneric("gene.table<-")
})

#' @rdname ASpediaFI-accessor
#' @export
setReplaceMethod("gene.table", "ASpediaFI", function(object, value) {
    outFI <- object
    outFI@gene.table <- value
    return(outFI)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("as.table", function(object, ...) {
    standardGeneric("as.table")
})

#' @rdname ASpediaFI-accessor
#' @export
setMethod("as.table", "ASpediaFI", function(object) {
    return(object@as.table)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("as.table<-", function(object, value) {
    standardGeneric("as.table<-")
})

#' @rdname ASpediaFI-accessor
#' @export
setReplaceMethod("as.table", "ASpediaFI", function(object, value) {
    outFI <- object
    outFI@as.table <- value
    return(outFI)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("pathway.table", function(object, ...) {
    standardGeneric("pathway.table")
})

#' @rdname ASpediaFI-accessor
#' @export
setMethod("pathway.table", "ASpediaFI", function(object) {
    return(object@pathway.table)
})

#' @rdname ASpediaFI-accessor
#' @export
setGeneric("pathway.table<-", function(object, value) {
    standardGeneric("pathway.table<-")
})

#' @rdname ASpediaFI-accessor
#' @export
setReplaceMethod("pathway.table", "ASpediaFI", function(object, value) {
    outFI <- object
    outFI@pathway.table <- value
    return(outFI)
})
