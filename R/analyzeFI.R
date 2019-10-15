#' Functional interaction analysis of AS events
#'
#' Analyze functional interactions of AS events using Discriminative
#' Random Walk with Restart (DRaWR). It runs a DRaWR on a
#' heterogeneous network containing genes, AS events, and pathways.
#' It then performs GSEA on gene sets related to query genes.
#'
#' @param object Object of class ASpediaFI
#' @param query a character vector or a data frame containing query genes
#' @param expr a \code{SummarizedExperiment} object or matrix containing gene
#' expression profiles (FPKM)
#' @param ppi an \code{igraph} object containing known interactions between
#' genes. If NULL, an \code{igraph} object containing human gene-gene
#' interactions will be used.
#' @param pathways a GMT file or a named list of pathway gene sets. If NULL,
#' a combined list of HALLMARK, KEGG, and REACTOME pathway gene sets will be
#' used.
#' @param restart a restart probability
#' @param num.folds the number of folds for cross-validation
#' @param num.feats the number of feature nodes to be retained in the final
#' subnetwork
#' @param low.expr Genes with mean expression below low.expr are excluded.
#' AS events for corresponding genes are also excluded.
#' @param low.var AS events with variance below low.var are excluded. If NULL,
#' top 10,000 variable events are used for analysis.
#' @param prop.na AS events with the higher proportion of missing values than
#' prop.na are excluded.
#' @param prop.extreme AS events with the higher proportion of extreme values
#' (0 or 1) than prop.extreme are excluded.
#' @param cor.threshold a pair of AS event and gene with Spearman's correlation
#' greather than cor.threshold are connected in a heterogeneous network.
#' @export
#' @references Blatti, C. et al. (2016). Characterizing
#' gene sets using discrminative random walks with restart on
#' heterogeneous biological networks. \emph{Bioinformatics}, 32.
#' @importFrom mGSZ geneSetsList
#' @return ASpediaFI object with results of functional interaction analysis
#' @examples
#' library(limma)
#' data(GSE114922.fpkm)
#' data(GSE114922.psi)
#' design <- cbind(WT = 1, MvsW = colData(GSE114922.psi)$condition == "MUT")
#' fit <- lmFit(log2(GSE114922.fpkm + 1), design = design)
#' fit <- eBayes(fit, trend = TRUE)
#' tt <- topTable(fit, number = Inf, coef = "MvsW")
#' query <- rownames(tt[tt$logFC > 1 &tt$P.Value < 0.1,])
#' head(query)
#' #
#' #GSE114922.ASpediaFI <- analyzeFI(GSE114922.ASpediaFI, query, GSE114922.fpkm)
analyzeFI <- function(object, query, expr, ppi = NULL,
                      pathways = NULL, restart = 0.7, num.folds = 5,
                      num.feats = 100, low.expr = 1, low.var = NULL,
                      prop.na = 0.05, prop.extreme = 1,
                      cor.threshold = 0.3){

    outFI <- object

    if(is.null(ppi)){
        ppi <- ppi.human
    }
    if(is.null(pathways)){
        pathways <- pathways.human
    }

    if(!is(pathways, "list")){
        pathways <- geneSetsList(pathways)
    }

    res <- analyze(query = query, psi = outFI@psi, expr = expr,
                   ppi = ppi, pathways = pathways, restart = restart,
                   num.folds = num.folds, num.feats = num.feats,
                   low.expr = low.expr, low.var = low.var,
                   prop.na = prop.na, prop.extreme = prop.extreme,
                   cor.threshold = cor.threshold)

    outFI@network <- res$network
    outFI@gene.table <- res$gene.table
    outFI@as.table <- res$as.table
    outFI@pathway.table <- res$pathway.table

    return(outFI)

}

