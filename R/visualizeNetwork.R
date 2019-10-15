#'Network visualization
#'
#'This function visualizes a heterogeneous network pertaining to specific
#'pathway.
#'
#'@param pathway.id name of pathway
#'@param net an igraph object containing a heterogeneous network produced
#'by DRaWR.
#'@param n number of genes and AS events to show
#'@keywords internal
#'@return plot of a pathway-specific network
#'@import igraph
#'@importFrom graphics par
#'@importFrom limma strsplit2
#'@noRd
visualizeNetwork <- function(pathway.id, net, n = 10){

    neighbor.genes <- neighbors(net, pathway.id, "all")
    neighbor.genes <- neighbor.genes$name[order(neighbor.genes$prob,
                                                decreasing = TRUE)][1:n]
    AS.vs <- V(net)[V(net)$type == "AS"]
    neighbor.AS <- unique(unlist(sapply(neighbor.genes, function(x){
        tmp <- neighbors(net, x, "all")
        tmp <- tmp$name[tmp$type == "AS"]
        tmp
    })))
    neighbor.AS <- AS.vs[AS.vs$name %in% neighbor.AS]
    neighbor.AS <- neighbor.AS$name[order(neighbor.AS$prob,
                                          decreasing = TRUE)][1:n]
    subnet <- suppressWarnings(subgraph(net, c(neighbor.genes,
                                               neighbor.AS)))

    par(mar = c(0, 0, 0, 0))
    V(subnet)$label <- ifelse(V(subnet)$type == "gene", V(subnet)$name,
                              apply(strsplit2(V(subnet)$name,
                                              split = ":")[,1:2], 1,
                                    function(x) paste(x, collapse = "\n")))
    V(subnet)$color <- ifelse(V(subnet)$type == "gene", "white", "lightblue")
    V(subnet)$shape <- ifelse(V(subnet)$type == "gene", "circle", "square")
    E(subnet)$lty <- ifelse(E(subnet)$type == "PPI", 1, 2)
    E(subnet)$weight <- 1

    plot(subnet, vertex.label.color = "black", vertex.size = 25,
         layout = layout_with_lgl(subnet, area = vcount(subnet)^3),
         edge.color = "black",
         edge.width = 2)

}
