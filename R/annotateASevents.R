#' AS event annotation
#'
#' Detect and annotate AS events from GTF. This function borrows code
#' from the \code{IVAS} package.
#'
#' @param object Object of class ASpediaFI
#' @param gtf.file an input GTF file
#' @param num.cores the number of cores for parallel processing
#' @export
#' @importFrom rtracklayer import.gff
#' @references Han, S. et al. (2017). Genome wide
#' discovery of genetic variants affecting alternative splicing
#' patterns in human using bioinformatics method. \emph{Genes & Genomics}, 39.
#' @return ASpediaFI object with a list of AS event annotations
#' @examples
#' fi <- new("ASpediaFI")
#' gtf <- system.file("extdata/GRCh38.subset.gtf", package = "ASpediaFI")
#' fi <- annotateASevents(fi, gtf.file = gtf, num.cores = 1)
#' sapply(events(fi), length)
#' head(events(fi)$SE)
annotateASevents <- function(object, gtf.file, num.cores = 1){
    outFI <- object
    as.list <- detect(gtf.file, num.cores)
    outFI@gtf <- import.gff(gtf.file)
    outFI@events <- annotate(as.list, outFI@gtf)
    return(outFI)
}
