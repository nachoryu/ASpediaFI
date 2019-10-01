setOldClass("igraph")

#' A Reference Class for functional interaction analysis of AS events
#'
#' @field samples a data frame containing the names, bam file paths, and
#' conditions of samples
#' @field events a list of annotated AS events extracted from GTF
#' @field psi a \code{SummarizedExperiment} object containing PSI of AS events
#' @field gtf a \code{GRanges} object containing GTF
#' @field network a heterogeneous network consisting of genes, AS events, and
#' pathways
#' @field gene.table a data frame containing gene nodes
#' @field as.table a data frame containing AS nodes
#' @field pathway.table a data frame containing Pathway nodes
#' @import igraph
#' @import GenomicRanges
#' @import methods
#' @importFrom mGSZ geneSetsList
#' @importFrom rtracklayer import.gff
#' @export ASpediaFI
#' @exportClass ASpediaFI
#' @examples
#' ###Instantiate ASpediaFI reference class
#' GSE114922.ASpediaFI <- ASpediaFI()
#'
#' ###Detect and annotate AS events from a subset of GTF file
#' #gtf <- system.file("extdata/GRCh38.subset.gtf", package = "ASpediaFI")
#'
#' ###Compute PSI values of AS events
#' #bamMT <- system.file("extdata/GSM3167287.subset.bam", package = "ASpediaFI")
#' #bamWT <- system.file("extdata/GSM3167290.subset.bam", package = "ASpediaFI")
#' #GSE114922.ASpediaFI$samples <- data.frame(name = c("GSM3167287", "GSM3167290"),
#' #                                         path = c(bamMT, bamWT),
#' #                                         condition = c("MT", "WT"))
#' #GSE114922.ASpediaFI$quantifyPSI(read.type = "paired", read.length = 100,
#' #                               insert.size = 300, min.reads = 3, num.cores = 3)
ASpediaFI <- setRefClass(
    "ASpediaFI",
    fields = list(samples = "data.frame",
                  events = "list",
                  psi = "SummarizedExperiment",
                  gtf = "GRanges",
                  network = "igraph",
                  gene.table = "data.frame",
                  as.table = "data.frame",
                  pathway.table = "data.frame"),
    methods = list(
        annotateASevents = function(gtf.file, num.cores = 1){
            "Detect and annotate AS events from GTF. This function borrows code
             from the \\code{IVAS} package.
             \\subsection{Arguments}{\\itemize{
                \\item{\\code{gtf.file} an input GTF file}
                \\item{\\code{num.cores} the number of cores for parallel
                       processing}}}
             \\subsection{References}{Han, S. et al. (2017). Genome wide
             discovery of genetic variants affecting alternative splicing
             patterns in human using bioinformatics method.
             \\emph{Genes & Genomics}, 39.}
            "
            as.list <- detect(gtf.file, num.cores)
            .self$gtf <- import.gff(gtf.file)
            .self$events <- annotate(as.list, .self$gtf)
        },
        quantifyPSI = function(read.type = "paired", read.length, insert.size,
                               min.reads = 3, num.cores = 1){
            "Compute PSI values of AS events. This function borrows code from
             the \\code{IMAS} package.
             \\subsection{Arguments}{\\itemize{
                 \\item{\\code{read.type} a type of RNA-seq reads ('single' or
                        'paired')}
                 \\item{\\code{read.length} read length}
                 \\item{\\code{insert.size} insert size}
                 \\item{\\code{min.reads} a minimum number of reads mapped to
                        a given exon}
                 \\item{\\code{num.cores} the number of cores for parallel
                        processing}}}
             \\subsection{References}{Han, S. et al. (2017). IMAS: Integrative
             analysis of Multi-omics data for Alternative Splicing.
             R package version 1.8.0.}
            "
            .self$psi <- quantify(.self$events, .self$samples, read.type,
                                  read.length, insert.size, min.reads, num.cores)
        },
        analyzeFI = function(query, expr, ppi = NULL, pathways = NULL,
                             restart = 0.7, num.folds = 5, num.feats = 100,
                             low.expr = 1, low.var = NULL, prop.na = 0.05,
                             prop.extreme = 1, cor.threshold = 0.3){
            "Analyze functional interactions of AS events using Discriminative
             Random Walk with Restart (DRaWR). It runs a DRaWR on a
             heterogeneous network containing genes, AS events, and pathways.
             It then performs GSEA on gene sets related to query genes.
             \\subsection{Arguments}{\\itemize{
                 \\item{\\code{query} a character vector or a data frame
                        containing query genes}
                 \\item{\\code{expr} a \\code{SummarizedExperiment} object or
                        matrix containing gene expression profiles (FPKM)}
                 \\item{\\code{ppi} an \\code{igraph} object containing known
                        interactions between genes. If NULL, an \\code{igraph}
                        object containing human gene-gene interactions will be
                        used.}
                 \\item{\\code{pathways} a GMT file or a named list of pathway
                        gene sets. If NULL, a combined list of HALLMARK, KEGG,
                        and REACTOME pathway gene sets will be used. }
                 \\item{\\code{restart} a restart probability}
                 \\item{\\code{num.folds} the number of folds for
                        cross-validation}
                 \\item{\\code{num.feats} the number of feature nodes to be
                        retained in the final subnetwork}
                 \\item{\\code{low.expr} Genes with mean expression below
                        low.expr are excluded. AS events for corresponding
                        genes are also excluded.}
                 \\item{\\code{low.var} AS events with variance below low.var
                        are excluded. If NULL, top 10,000 variable events are
                        used for analysis.}
                 \\item{\\code{prop.na} AS events with the higher proportion of
                        missing values than prop.na are excluded.}
                 \\item{\\code{prop.extreme} AS events with the higher
                        proportion of extreme values (0 or 1) than prop.extreme
                        are excluded.}
                 \\item{\\code{cor.threshold} a pair of AS event and gene with
                        Spearman's correlation greather than cor.threshold are
                        connected in a heterogeneous network.}}}
             \\subsection{References}{Blatti, C. et al. (2016). Characterizing
             gene sets using discrminative random walks with restart on
             heterogeneous biological networks. \\emph{Bioinformatics}, 32.}
            "

            if(is.null(ppi)){
                ppi <- ppi.human
            }
            if(is.null(pathways)){
                pathways <- pathways.human
            }

            if(!is(pathways, "list")){
                pathways <- geneSetsList(pathways)
            }
            res <- analyze(query = query, psi = .self$psi, expr = expr,
                           ppi = ppi, pathways = pathways, restart = restart,
                           num.folds = num.folds, num.feats = num.feats,
                           low.expr = low.expr, low.var = low.var,
                           prop.na = prop.na, prop.extreme = prop.extreme,
                           cor.threshold = cor.threshold)
            .self$network <- res$network
            .self$gene.table <- res$gene.table
            .self$as.table <- res$as.table
            .self$pathway.table <- res$pathway.table
        },
        visualize = function(node, zoom = NULL, n = NULL){
            "Visualize AS event or pathway. If an AS event node is given, the
             function modified from the \\code{plotTranscripts} function in the
             \\code{maser} package is used to visualize the event. If a
             pathway node is given, a subnetwork pertaining to the pathway is
             visualized.
             \\subsection{Arguments}{\\itemize{
                 \\item{\\code{node} the name of AS event or pathway}
                 \\item{\\code{zoom} a logical to determine if genomic
                        coordinates are zoomed (for AS event visualization)}
                 \\item{\\code{n} the number of genes and AS events to be shown
                        (for pathway visualization)}}}
             \\subsection{References}{Veiga, D. (2019). maser: Mapping
             Alternative Splicing Events to pRoteins. R package version 1.2.0.
             https://github.com/DiogoVeiga/maser}
            "
            if(node %in% .self$as.table$EventID){
                visualizeEvent(node, .self$gtf, .self$psi, zoom)
            }
            if(node %in% .self$pathway.table$Pathway){
                visualizeNetwork(node, .self$network, n)
            }
        },
        exportNetwork = function(node = NULL, file){
            "Export a subnetwork pertaining to the given pathway to GML format
             which can be used in Cytoscape. If no pathway is given, the entire
             final subnetwork is exported.
            \\subsection{Arguments}{\\itemize{
            \\item{\\code{node} the name of pathway. If NULL, the entire
                   subnetwork is exported.}
            \\item{\\code{file} the file name to export the network}}}
            "
            if(is.null(node)){
                write_graph(.self$network, file, "gml")
            }
            else{
                neighbor.genes <- neighbors(.self$network, node, "all")
                AS.vs <- V(.self$network)[V(.self$network)$type == "AS"]
                neighbor.AS <- unique(unlist(sapply(neighbor.genes, function(x){
                    tmp <- neighbors(.self$network, x, "all")
                    tmp <- tmp$name[tmp$type == "AS"]
                    tmp
                })))
                neighbor.AS <- AS.vs[AS.vs$name %in% neighbor.AS]
                subnet <- suppressWarnings(subgraph(.self$network,
                                                    c(neighbor.genes,
                                                      neighbor.AS)))
                write_graph(subnet, file, "gml")
            }
        }
    )
)


