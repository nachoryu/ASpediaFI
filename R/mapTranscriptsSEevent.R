#' Map transcripts to SE events
#'
#' This function maps transcripts to SE events
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
mapTranscriptsSEevent <- function(eventGr, gtf_exons) {
    ovl.e1 <- findOverlaps(eventGr$exon_upstream, gtf_exons, type = "equal")
    ovl.e2 <- findOverlaps(eventGr$exon_target, gtf_exons, type = "equal")
    ovl.e3 <- findOverlaps(eventGr$exon_downstream, gtf_exons, type = "equal")
    mytx.ids.e1 <- gtf_exons$transcript_id[subjectHits(ovl.e1)]
    mytx.ids.e2 <- gtf_exons$transcript_id[subjectHits(ovl.e2)]
    mytx.ids.e3 <- gtf_exons$transcript_id[subjectHits(ovl.e3)]
    intron.skipping <- GRanges(seqnames = seqnames(eventGr$exon_target),
                        ranges = IRanges(start = end(eventGr$exon_upstream) + 1,
                        end = start(eventGr$exon_downstream) - 1),
                        strand = strand(eventGr$exon_target))
    intron.inclusion <- GRanges(seqnames = seqnames(eventGr$exon_target),
                    ranges = IRanges(start = c(end(eventGr$exon_upstream) + 1,
                    end(eventGr$exon_target) + 1),
                    end = c(start(eventGr$exon_target) - 1,
                    start(eventGr$exon_downstream) - 1)),
                    strand = strand(eventGr$exon_target))
    ovl.intron.inclusion <-
        findOverlaps(intron.inclusion, gtf_exons, type = "any")
    mytx.ids.intron.inclusion <-
        gtf_exons$transcript_id[subjectHits(ovl.intron.inclusion)]
    ovl.intron.skipping <-
        findOverlaps(intron.skipping, gtf_exons, type = "any")
    mytx.ids.intron.skipping <-
        gtf_exons$transcript_id[subjectHits(ovl.intron.skipping)]
    mytx.ids.3exons <- intersect(mytx.ids.e1, mytx.ids.e3)
    mytx.ids.3exons <- intersect(mytx.ids.3exons, mytx.ids.e2)
    mytx.ids.2exons <- intersect(mytx.ids.e1, mytx.ids.e3)
    mytx.ids.3exons <- setdiff(mytx.ids.3exons, mytx.ids.intron.inclusion)
    mytx.ids.2exons <- setdiff(mytx.ids.2exons, mytx.ids.intron.skipping)
    tx_ids <- list()
    tx_ids[["txn_3exons"]] <- mytx.ids.3exons
    tx_ids[["txn_2exons"]] <- mytx.ids.2exons
    return(tx_ids)
}
