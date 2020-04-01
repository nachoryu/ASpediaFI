#' Functional interaction analysis of AS events
#'
#' This function analyzes comprehensvie functional interactions of AS events
#' using Discriminative Random Walk with Restart (DRaWR). It runs a DRaWR on a
#' heterogeneous network containing genes, AS events, and gene sets
#' (or pathways). It then performs GSEA on gene sets related to query genes.
#'
#' @param query a character vector or data frame containing query genes
#' @param expr a SummarizedExperiment object or matrix containing gene
#' expression profiles (FPKM)
#' @param psi a SummarizedExperiment object containing PSI of AS events
#' and an optional data frame describing samples
#' @param pathways a named list of pathway gene sets. If NULL, a combined
#' list of HALLMARK, KEGG, and REACTOME pathway gene sets will be used.
#' @param ppi an \code{igraph} object containing known interactions between
#' genes. If NULL, an \code{igraph} object containing human gene-gene
#' interactions will be used.
#' @param restart a restart probability
#' @param num.folds number of folds for cross-validation
#' @param num.feats number of feature nodes to be kept in the final network
#' @param low.expr Genes with mean expression below low.expr are excluded. AS
#' events for corresponding genes are also excluded.
#' @param low.var AS events with variance below low.var are excluded. If NULL,
#' top 10,000 variable events are used for analysis.
#' @param prop.na AS events with the higher proportion of missing values than
#' prop.na are excluded.
#' @param prop.extreme AS events with the higher proportion of extreme values
#' (0 or 1) than prop.extreme are excluded.
#' @param cor.threshold a pair of AS event and gene with Spearman's correlation
#' greather than cor.threshold are connected in a heterogeneous network.
#' @return A matrix containing PSI of AS events
#' @details This function wraps around the \code{DRaWR} function in the
#' \code{DRaWR} package with data processing steps.
#' @references Blatti, C. et al. (2016). Characterizing gene sets using
#' discrminative random walks with restart on heterogeneous biological networks
#' \emph{Bioinformatics}, 32.
#' @keywords internal
#' @import igraph
#' @importFrom limma strsplit2
#' @importFrom reshape2 melt
#' @importFrom fgsea fgsea
#' @importFrom stats t.test var cor
#' @importFrom e1071 impute
#' @importFrom methods is
#' @noRd
analyze <- function(query, psi, expr, pathways = NULL, ppi = NULL,
                        restart = 0.7,
                        num.folds = 5, num.feats = 100, low.expr = 1,
                        low.var = NULL, prop.na = 0.05, prop.extreme = 1,
                        cor.threshold = 0.3) {
    psi.mat <- assays(psi)[[1]]

    if (is(expr, "SummarizedExperiment")) {
        exp.mat <- assays(expr)[[1]]
    } else {
        exp.mat <- expr
    }

    if (!all.equal(colnames(psi.mat), colnames(exp.mat))) {
        stop("Sample names do not match.")
    }

    # Exclude genes with low expression
    exp.mat <- rbind(exp.mat[rowMeans(exp.mat) > low.expr, ])

    # Exclude AS events with high proportions of missing values
    psi.mat <- rbind(psi.mat[apply(psi.mat, 1,
                        function(x) sum(is.na(x))/ncol(psi.mat)) < prop.na, ])
    psi.mat <- t(impute(t(psi.mat), what = "median"))

    # Exclude AS events with high proportions of extreme values
    psi.mat <- rbind(psi.mat[apply(psi.mat, 1,
            function(x) sum(x == 0 | x == 1)/ncol(psi.mat)) < prop.extreme, ])

    # Exclude AS events on genes with low expression
    psi.mat <- rbind(psi.mat[strsplit2(rownames(psi.mat),
                        split = ":")[, 1] %in% rownames(exp.mat), ])

    # Exclude AS events with low variance
    psi.var <- apply(psi.mat, 1, var)
    if (is.null(low.var)) {
        psi.mat <- rbind(psi.mat[names(sort(psi.var,
            decreasing = TRUE))[seq_len(min(nrow(psi.mat), 10000))], ])
    } else {
        psi.mat <- rbind(psi.mat[psi.var >= low.var, ])
    }

    # Network component 1: gene-gene
    exp.mat <- rbind(exp.mat[intersect(rownames(exp.mat),
                        unique(names(V(ppi)))), ])
    universe.genes <- rownames(exp.mat)
    ppi.sub <- suppressWarnings(subgraph(ppi, universe.genes))

    gg <- data.frame(as_edgelist(ppi.sub), stringsAsFactors = FALSE)
    coexpr <- cor(t(exp.mat), method = "pearson")
    gg$X3 <- apply(gg, 1, function(x) coexpr[x[1], x[2]])
    gg <- gg[!is.na(gg$X3), ]
    gg$X3 <- abs(gg$X3)
    gg$X4 <- "PPI"
    rm(coexpr)

    # Network component 2: gene-AS
    exprpsi <- t(rbind(log2(exp.mat + 1), psi.mat))
    corr <- cor(exprpsi, method = "spearman")
    corr.ep <- corr[(nrow(exp.mat) + 1):nrow(corr), seq_len(nrow(exp.mat))]
    corr.melt <- melt(corr.ep)
    corr.melt <- data.frame(corr.melt, stringsAsFactors = FALSE)
    corr.melt <- corr.melt[!is.na(corr.melt$value), ]
    corr.melt$value <- abs(corr.melt$value)

    gas <- corr.melt[corr.melt$value > cor.threshold, ]
    gas$X4 <- "AS"
    rm(corr, corr.ep, corr.melt, exprpsi)

    # Network component 3: gene-pathway
    num.neighbors <- vapply(pathways, function(x) sum(universe.genes %in% x),
                                numeric(1))
    pathways <- pathways[num.neighbors > 10]
    gpw <- NULL
    for (gid in universe.genes) {
        pwid <- names(pathways)[vapply(pathways, function(x) gid %in% x,
                                        logical(1))]
        if (length(pwid) > 0) {
            gpw <- rbind(gpw, cbind(pwid, gid))
        }
    }
    gpw <- data.frame(gpw, stringsAsFactors = FALSE)
    gpw$X3 <- 1
    gpw$X4 <- "Pathway"

    # Combine into a heterogeneous network
    colnames(gg) <- colnames(gas) <- colnames(gpw) <- c("X1", "X2", "X3", "X4")
    hetnet <- data.frame(rbind(gg, gas, gpw), stringsAsFactors = FALSE)

    # Prepare inputs for DRaWR
    universe <- data.frame(gene = universe.genes, weight = 1,
                            stringsAsFactors = FALSE)

    if (is.null(dim(query))) {
        query <- cbind(as.character(query), rep(1, length(query)))
    }
    query <- data.frame(query, stringsAsFactors = FALSE)
    colnames(query) <- c("gene", "weight")
    query <- query[query$gene %in% universe.genes, ]

    # Run DRaWR
    drawr.result <- drawr(query.genes = query, universe = universe,
                            network = hetnet, restart = restart,
                            num.feats = num.feats, num.folds = num.folds)

    # Gene table
    gene.table <- drawr.result$genes
    gene.table$GeneSymbol <- gene.table$node
    gene.table$StatP <- as.numeric(gene.table$prob)
    gene.table$PermPvalue <- as.numeric(gene.table$pval)
    gene.table <- gene.table[, c("GeneSymbol", "StatP", "PermPvalue")]
    gene.nodes <- gene.table$GeneSymbol

    # AS event table
    as.table <- drawr.result$features[drawr.result$features$type == "AS", ]
    as.table$EventID <- as.table$node
    as.table$StatP <- as.numeric(as.table$prob)
    as.table$GeneSymbol <- strsplit2(as.table$EventID, split = ":")[, 1]
    as.table$EventType <- strsplit2(as.table$EventID, split = ":")[, 2]
    as.table$Rank <- seq_len(nrow(as.table))
    as.table <- as.table[, c("EventID", "GeneSymbol", "EventType", "Rank",
                                "StatP")]
    rownames(as.table) <- seq_len(nrow(as.table))
    as.nodes <- as.table$EventID

    # Pathway table
    pathway.table <- drawr.result$features[drawr.result$features$type ==
                                                "Pathway", ]
    pathway.table$prob <- as.numeric(pathway.table$prob)
    pathway.table$Pathway <- pathway.table$node
    pathway.table$StatP <- pathway.table$prob
    pathway.table$Rank <- seq_len(nrow(pathway.table))
    pathway.table <- pathway.table[, c("Pathway", "Rank", "StatP")]
    pathway.nodes <- pathway.table$Pathway

    # Final network
    gg <- gg[gg$X1 %in% gene.nodes & gg$X2 %in% gene.nodes, ]
    gas <- gas[gas$X1 %in% as.nodes & gas$X2 %in% gene.nodes, ]
    gpw <- gpw[gpw$X1 %in% pathway.nodes & gpw$X2 %in% gene.nodes, ]

    edges <- rbind(gg, gas, gpw)
    colnames(edges) <- c("src", "target", "weight", "type")
    typetable <- aggregate(weight ~ type, data = edges, FUN = sum)
    rownames(typetable) <- as.character(typetable[, 1])
    edges$weight <- edges$weight/typetable[as.character(edges$type), "weight"]
    vertices <- data.frame(name = c(gene.nodes, as.nodes, pathway.nodes),
                                type = c(rep("gene", length(gene.nodes)),
                                            rep("AS", length(as.nodes)),
                                            rep("Pathway",
                                                length(pathway.nodes))),
                                prob = c(gene.table$StatP,
                                            as.table$StatP,
                                            pathway.table$StatP))

    drawr.net <- graph_from_data_frame(edges, directed = FALSE,
                                        vertices = vertices)

    # Gene set analysis for pathway nodes
    if ("condition" %in% colnames(colData(psi))) {
        condition <- colData(psi)$condition
        tstat <- apply(log2(exp.mat + 1), 1,
                        function(x) t.test(x ~ factor(condition))$statistic)
        gsea <- suppressWarnings(fgsea(pathways, tstat, 10000))
        gsea <- data.frame(gsea, stringsAsFactors = FALSE)
        rownames(gsea) <- gsea$pathway
        pathway.table$Pvalue <- gsea[pathway.nodes, "pval"]
        pathway.table$Adj.Pvalue <- gsea[pathway.nodes, "padj"]
        pathway.table$EnrichmentScore <- gsea[pathway.nodes, "ES"]
        pathway.table$NormalizedEnrichmentScore <- gsea[pathway.nodes, "NES"]
    }

    pathway.table$Size <- vapply(pathways, length, numeric(1))[pathway.nodes]
    pathway.table$Count <- vapply(pathway.nodes,
                            function(x) sum(universe.genes %in% pathways[[x]]),
                                numeric(1))
    pathway.table$AvgRank <- vapply(pathway.nodes,
            function(x) mean(which(drawr.result$genes$node %in% pathways[[x]])),
                numeric(1))
    pathway.table$NumEvents <- vapply(pathway.nodes,
            function(x) sum(unique(gas$X1[gas$X2 %in% pathways[[x]]]) %in%
                                                                    as.nodes),
                numeric(1))
    rownames(pathway.table) <- seq_len(nrow(pathway.table))
    return(list(network = drawr.net, gene.table = gene.table,
                as.table = as.table, pathway.table = pathway.table))
}
