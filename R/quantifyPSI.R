#' AS event quantification
#'
#' Compute PSI values of AS events. This function borrows code
#' from the \code{IMAS} package.
#'
#' @param object Object of class ASpediaFI
#' @param read.type a type of RNA-seq reads ('single' or 'paired')
#' @param read.length read length
#' @param insert.size insert size
#' @param min.reads a minimum number of reads mapped to a given exon
#' @param num.cores the number of cores for parallel processing
#' @export
#' @references Han, S. et al. (2017). IMAS: Integrative
#' analysis of Multi-omics data for Alternative Splicing.
#' R package version 1.8.0.
#' @return ASpediaFI object with PSI values
#' @examples
#' bamWT <- system.file('extdata/GSM3167290.subset.bam', package = 'ASpediaFI')
#' GSE114922.ASpediaFI <- ASpediaFI(
#'     sample.names = 'GSM3167290',
#'     bam.files = bamWT, conditions = 'WT'
#' )
#' \dontrun{
#' GSE114922.ASpediaFI <- quantifyPSI(GSE114922.ASpediaFI,
#'     read.type = 'paired',
#'     read.length = 100, insert.size = 300,
#'     min.reads = 3, num.cores = 1
#' )
#' }
quantifyPSI <- function(object, read.type = "paired", read.length,
                            insert.size, min.reads, num.cores = 1) {
    outFI <- object
    psi(outFI) <- quantify(events(outFI), samples(outFI), read.type,
                            read.length, insert.size, min.reads, num.cores)
    return(outFI)
}
