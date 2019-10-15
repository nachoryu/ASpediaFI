#' Export network to GML format
#'
#' Export a subnetwork pertaining to the given pathway to GML format
#' which can be used in Cytoscape. If no pathway is given, the entire
#' final subnetwork is exported.
#'
#' @param object Object of class ASpediaFI
#' @param node the name of pathway. If NULL, the entire subnetwork is exported.
#' @param file the file name to export the network
#' @export
#' @return a GML file containing a subnetwork
#' @import igraph
#' @examples
#' library(igraph)
#' fi <- new("ASpediaFI", network = make_empty_graph(n = 0))
#' exportNetwork(fi, node = NULL, file = "empty.gml")
exportNetwork <- function(object, node = NULL, file){
    if(is.null(node)){
        write_graph(object@network, file, "gml")
    }
    else{
        neighbor.genes <- neighbors(object@network, node, "all")
        AS.vs <- V(object@network)[V(object@network)$type == "AS"]
        neighbor.AS <- unique(unlist(sapply(neighbor.genes, function(x){
            tmp <- neighbors(object@network, x, "all")
            tmp <- tmp$name[tmp$type == "AS"]
            tmp
        })))
        neighbor.AS <- AS.vs[AS.vs$name %in% neighbor.AS]
        subnet <- suppressWarnings(subgraph(object@network,
                                            c(neighbor.genes,
                                              neighbor.AS)))
        write_graph(subnet, file, "gml")
    }
}
