setOldClass("igraph")
.ASpediaFI <- setClass("ASpediaFI", slots = c(samples = "data.frame",
                        events = "list", psi = "SummarizedExperiment",
                        gtf = "GRanges", network = "igraph",
                        gene.table = "data.frame", as.table = "data.frame",
                        pathway.table = "data.frame"))

#' @title ASpediaFI class
#'
#' @description ASpediaFI class is a wrapper of ASpediaFI functionalities and
#' a container of inputs and outputs.
#'
#' @param sample.names a character vector of sample names (or IDs)
#' @param bam.files a character vector of paths to RNA-Seq BAM files
#' @param conditions a vector of sample conditions (e.g. mutation status)
#'
#' @section Slots:
#' \describe{
#' \item{}{
#'     \code{samples}: a data frame containing information about samples. The
#'     first three columns should be names, BAM file paths, and conditions.
#' }
#' \item{}{
#'     \code{events}: a list of AS events extracted from a GTF file.
#' }
#' \item{}{
#'     \code{gtf}: a \code{GRanges} object containing genomic features
#'     extracted from a GTF file.
#' }
#' \item{}{
#'     \code{psi}: a \code{SummarizedExperiment} object containing AS event
#'     quantification
#' }
#' \item{}{
#'     \code{network}: an \code{igraph} object containing a query-specific
#'     subnetwork as a result of DRaWR.
#' }
#' \item{}{
#'     \code{gene.table, as.table, pathway.table}: data
#'     frames containing gene nodes, AS event nodes, and pathway nodes.
#' }
#' }
#'
#' @section Accessors:
#'
#' In the following, 'x' represents a \code{ASpediaFI} object:
#' \describe{
#' \item{}{
#'     \code{samples(x), samples(x) <- value}: get or set sample information.
#'     value must be a data frame containing sample information.
#' }
#' \item{}{
#'     \code{events(x), events(x) <- value}: get or set AS event annotations.
#'     value must be a list of annotations.
#' }
#' \item{}{
#'     \code{gtf(x), gtf(x) <- value}: get or set a \code{GRanges} object
#'     containing GTF. value must be a \code{GRanges} object.
#' }
#' \item{}{
#'     \code{psi(x), psi(x) <- value}: get or set PSI values.
#'     value must be a \code{SummarizedExperiment} object.
#' }
#' \item{}{
#'     \code{network(x), network(x) <- value}: get or set final subnetwork.
#'     value must be an \code{igraph} object.
#' }
#' \item{}{
#'     \code{gene.table(x), gene.table(x) <- value}: get or set gene node
#'     tables. value must be a data frame containing information about gene
#'     nodes.
#' }
#' \item{}{
#'     \code{as.table(x), as.table(x) <- value}: get or set AS node tables.
#'     value must be a data frame containing information about AS nodes.
#' }
#' \item{}{
#'     \code{pathway.table(x), pathway.table(x) <- value}: get or set pathway
#'     node tables. value must be a data frame containing information about
#'     pathway nodes.
#' }
#' }
#' @examples
#' bamWT <- system.file('extdata/GSM3167290.subset.bam', package = 'ASpediaFI')
#' GSE114922.ASpediaFI <- ASpediaFI(
#'     sample.names = 'GSM3167290',
#'     bam.files = bamWT, conditions = 'WT'
#' )
#' @name ASpediaFI-class
#' @export ASpediaFI
#' @exportClass ASpediaFI
#' @return ASpediaFI object
#' @import methods
#' @import igraph
#' @importFrom GenomicRanges GRanges
#' @import SummarizedExperiment
ASpediaFI <- function(sample.names, bam.files, conditions) {
    sample.info <- data.frame(name = sample.names, path = bam.files,
                            condition = conditions, stringsAsFactors = FALSE)
    new("ASpediaFI", samples = sample.info, events = list(),
            psi = SummarizedExperiment(), gtf = GRanges(),
            network = make_empty_graph(), gene.table = data.frame(),
            as.table = data.frame(), pathway.table = data.frame())
}
