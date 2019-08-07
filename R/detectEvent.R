#'Detection of AS events
#'
#'This function detects AS events from an input GTF.
#'
#'@param file Input GTF or GFF3 file.
#'@param Ncor The number of cores for parallel processing
#'@return A list of AS events
#'@details This function wraps around \code{makeTxDbFromGFF} in the
#'\code{GenomicFeatures} package and \code{Splicingfinder} in the
#'\code{IVAS} package.
#'@references Han, A. et al. (2017). Genome wide discovery of genetic variants
#'affecting alternative splicing patterns in human using bioinformatics method.
#'\emph{Genes & Genomics}, 39.
#'@importFrom GenomicFeatures makeTxDbFromGFF
#'@export
detectEvent <- function(file, Ncor){
    gtfdb <- makeTxDbFromGFF(file)
    saveDb(gtfdb, file = "gtf.db")
    aslist.raw <- detectEventFromTxDb(gtfdb, Ncor)
    gtfdb <- loadDb("gtf.db")
    aslist.filtered <- filterEvent(aslist.raw, gtfdb)
    aslist.annotated <- annotateEvent(aslist.filtered)
}
