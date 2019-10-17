#' Map transcripts to A5SS events
#'
#' This function maps transcripts to A5SS events
#'
#' @param eventGr A GRangesList object of AS event
#' @param gtf_exons A GRanges object containing exons in gtf
#' @return a list of transcripts
#' @details This function is borrowed from the \code{maser} package.
#' @references Veiga, D. (2019). maser: Mapping Alternative Splicing Events
#' to pRoteins. R package version 1.2.0. https://github.com/DiogoVeiga/maser
#' @keywords internal
#' @importFrom GenomicRanges findOverlaps GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors subjectHits
#' @noRd
mapTranscriptsA5SSevent <- function(eventGr, gtf_exons) {
    ovl.e1 <- findOverlaps(eventGr$exon_short, gtf_exons, type = "equal")
    ovl.e2 <- findOverlaps(eventGr$exon_long, gtf_exons, type = "equal")
    ovl.e3 <- findOverlaps(eventGr$exon_flanking, gtf_exons, type = "equal")
    mytx.ids.e1 <- gtf_exons$transcript_id[subjectHits(ovl.e1)]
    mytx.ids.e2 <- gtf_exons$transcript_id[subjectHits(ovl.e2)]
    mytx.ids.e3 <- gtf_exons$transcript_id[subjectHits(ovl.e3)]
    intron.short <- GRanges(seqnames = seqnames(eventGr$exon_short),
                        ranges = IRanges(start = end(eventGr$exon_short) + 1,
                        end = start(eventGr$exon_flanking) - 1),
                        strand = strand(eventGr$exon_short))
    intron.long <- GRanges(seqnames = seqnames(eventGr$exon_long),
                        ranges = IRanges(start = end(eventGr$exon_long) + 1,
                        end = start(eventGr$exon_flanking) - 1),
                        strand = strand(eventGr$exon_long))
    ovl.intron.short <- findOverlaps(intron.short, gtf_exons, type = "any")
    mytx.ids.intron.short <-
        gtf_exons$transcript_id[subjectHits(ovl.intron.short)]
    ovl.intron.long <- findOverlaps(intron.long, gtf_exons, type = "any")
    mytx.ids.intron.long <-
        gtf_exons$transcript_id[subjectHits(ovl.intron.long)]
    mytx.ids.short <- intersect(mytx.ids.e1, mytx.ids.e3)
    mytx.ids.long <- intersect(mytx.ids.e2, mytx.ids.e3)
    mytx.ids.short <- setdiff(mytx.ids.short, mytx.ids.intron.short)
    mytx.ids.long <- setdiff(mytx.ids.long, mytx.ids.intron.long)
    tx_ids <- list()
    tx_ids[["txn_short"]] <- mytx.ids.short
    tx_ids[["txn_long"]] <- mytx.ids.long
    return(tx_ids)
}
