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
#'@importFrom GenomicFeatures makeTxDbFromGFF exonsBy intronsByTranscript
#'@importFrom biomaRt select
#'@export
detectEvent <- function(file, Ncor){
    gtfdb <- makeTxDbFromGFF(file)
    total.exon.range <- exonsBy(gtfdb, by = "tx")
    total.intron.range <- intronsByTranscript(gtfdb)
    txTable <- select(gtfdb, keys = names(total.exon.range),
                      columns = c("TXCHROM", "TXNAME", "GENEID",
                                  "TXSTART", "TXEND", "TXSTRAND"),
                      keytype = "TXID")
    aslist.raw <- detectEventFromTxDb(gtfdb, Ncor)
    aslist.filtered <- filterEvent(aslist.raw, total.exon.range,
                                   total.intron.range, txTable)
    aslist.annotated <- annotateEvent(aslist.filtered)
}
