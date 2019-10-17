#' Example dataset containing PSI values
#'
#' A SummarizedExperiment containing PSI values of 5,000 AS events.
#' We downloaded
#' RNA-Seq reads of 82 MDS patients from GEO database (GSE114922),
#' aligned with STAR, and computed PSI values using rMATS. AS events with a lot
#' of missing values or extreme values, or those on genes with low expression
#' were filtered out and 5,000 most variable AS events were selected.
#'
#' @name GSE114922.psi
#' @docType data
#' @references Pellagatti, A. et al. (2018). Impact of spliceosome mutations on
#' RNA splicing in myelodysplasia: dysregulated genes/pathways and clinical
#' associations. \emph{Blood}, 132.
#' @examples
#' data(GSE114922.psi)
"GSE114922.psi"

#' Example gene expression dataset
#'
#' A matrix containing gene expression values. We downloaded RNA-Seq reads of
#' 82 MDS patients from GEO database (GSE114922), aligned with STAR, and
#' obtained FPKM values using RSEM.
#'
#' @name GSE114922.fpkm
#' @docType data
#' @references Pellagatti, A. et al. (2018). Impact of spliceosome mutations on
#' RNA splicing in myelodysplasia: dysregulated genes/pathways and clinical
#' associations. \emph{Blood}, 132.
#' @examples
#' data(GSE114922.fpkm)
"GSE114922.fpkm"
