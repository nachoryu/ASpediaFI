#'Map transcripts to RI events
#'
#'This function maps transcripts to RI events
#'
#'@param eventGr A GRangesList object of AS event
#'@param gtf_exons A GRanges object containing exons in gtf
#'@return a list of transcripts
#'@details This function is borrowed from the \code{maser} package.
#'@references Veiga, D. (2019). maser: Mapping Alternative Splicing Events
#'to pRoteins. R package version 1.2.0. https://github.com/DiogoVeiga/maser
#'@keywords internal
#'@importFrom GenomicRanges findOverlaps GRanges
#'@importFrom IRanges IRanges
#'@importFrom S4Vectors subjectHits
mapTranscriptsRIevent <- function (eventGr, gtf_exons){
    ovl.e1 <- findOverlaps(eventGr$exon_upstream,
                           gtf_exons, type = "equal")
    ovl.e2 <- findOverlaps(eventGr$exon_ir, gtf_exons,
                           type = "equal")
    ovl.e3 <- findOverlaps(eventGr$exon_downstream,
                           gtf_exons, type = "equal")
    mytx.ids.e1 <- gtf_exons$transcript_id[subjectHits(ovl.e1)]
    mytx.ids.e2 <- gtf_exons$transcript_id[subjectHits(ovl.e2)]
    mytx.ids.e3 <- gtf_exons$transcript_id[subjectHits(ovl.e3)]
    intron <- GRanges(seqnames = seqnames(eventGr$exon_ir),
                      ranges = IRanges(start = end(eventGr$exon_upstream) +
                                           1, end = start(eventGr$exon_downstream) - 1),
                      strand = strand(eventGr$exon_ir))
    ovl.intron <- findOverlaps(intron, gtf_exons,
                               type = "any")
    mytx.ids.intron <- gtf_exons$transcript_id[subjectHits(ovl.intron)]
    tx.ids.nonRetention <- intersect(mytx.ids.e1, mytx.ids.e3)
    tx.ids.nonRetention <- setdiff(tx.ids.nonRetention, mytx.ids.intron)
    tx.ids.Retention <- mytx.ids.e2
    tx_ids <- list()
    tx_ids[["txn_nonRetention"]] <- tx.ids.nonRetention
    tx_ids[["txn_retention"]] <- tx.ids.Retention
    return(tx_ids)
}
