#' AS event and pathway visualization
#'
#' Visualize AS event or pathway. If an AS event node is given, the
#' function modified from the \code{plotTranscripts} function in the
#' \code{maser} package is used to visualize the event. If a
#' pathway node is given, a subnetwork pertaining to the pathway is
#' visualized.
#'
#' @param object Object of class ASpediaFI
#' @param node the name of AS event or pathway
#' @param zoom a logical to determine if genomic coordinates are zoomed
#' (for AS event visualization)
#' @param n the number of genes and AS events to be shown (for pathway
#' visualization)
#' @export
#' @return a plot demonstrating AS event or pathway
#' @references Veiga, D. (2019). maser: Mapping
#' Alternative Splicing Events to pRoteins. R package version 1.2.0.
#' https://github.com/DiogoVeiga/maser
#' @examples
#' \dontrun{
#' # Visualize AS event
#' visualize(GSE114922.ASpediaFI,
#'     node = as.table(GSE114922.ASpediaFI)$EventID[1],
#'     zoom = FALSE
#' )
#'
#' # Visualize pathway
#' visualize(GSE114922.ASpediaFI, node = 'HALLMARK_HEME_METABOLISM', n = 10)
#' }
visualize <- function(object, node, zoom = NULL, n = NULL) {
    if (node %in% as.table(object)$EventID) {
        visualizeEvent(node, gtf(object), psi(object), zoom)
    }
    if (node %in% pathway.table(object)$Pathway) {
        visualizeNetwork(node, network(object), n)
    }
}
