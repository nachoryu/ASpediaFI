#'Map transcripts to MXE events
#'
#'This function maps transcripts to MXE events
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
#'@noRd
mapTranscriptsMXEevent <- function (eventGr, gtf_exons){
    ovl.e1 <- findOverlaps(eventGr$exon_1, gtf_exons,
                           type = "equal")
    ovl.e2 <- findOverlaps(eventGr$exon_2, gtf_exons,
                           type = "equal")
    ovl.e3 <- findOverlaps(eventGr$exon_upstream,
                           gtf_exons, type = "equal")
    ovl.e4 <- findOverlaps(eventGr$exon_downstream,
                           gtf_exons, type = "equal")
    mytx.ids.e1 <- gtf_exons$transcript_id[subjectHits(ovl.e1)]
    mytx.ids.e2 <- gtf_exons$transcript_id[subjectHits(ovl.e2)]
    mytx.ids.e3 <- gtf_exons$transcript_id[subjectHits(ovl.e3)]
    mytx.ids.e4 <- gtf_exons$transcript_id[subjectHits(ovl.e4)]
    intron.mxe.exon1 <- GRanges(seqnames = seqnames(eventGr$exon_1),
                                ranges = IRanges(start = c(end(eventGr$exon_upstream) +
                                                               1, end(eventGr$exon_1) + 1),
                                                 end = c(start(eventGr$exon_1) -
                                                             1,
                                                         start(eventGr$exon_downstream) - 1)),
                                strand = strand(eventGr$exon_1))
    intron.mxe.exon2 <- GRanges(seqnames = seqnames(eventGr$exon_2),
                                ranges = IRanges(start = c(end(eventGr$exon_upstream) +
                                                               1, end(eventGr$exon_2) + 1),
                                                 end = c(start(eventGr$exon_2) -
                                                             1,
                                                         start(eventGr$exon_downstream) - 1)),
                                strand = strand(eventGr$exon_2))
    ovl.mxe.exon1 <- findOverlaps(intron.mxe.exon1,
                                  gtf_exons, type = "any")
    mytx.ids.intron1 <- gtf_exons$transcript_id[subjectHits(ovl.mxe.exon1)]
    ovl.mxe.exon2 <- findOverlaps(intron.mxe.exon2,
                                  gtf_exons, type = "any")
    mytx.ids.intron2 <- gtf_exons$transcript_id[subjectHits(ovl.mxe.exon2)]

    mytx.ids.mxe.exon1 <- intersect(mytx.ids.e3, mytx.ids.e4)
    mytx.ids.mxe.exon1 <- intersect(mytx.ids.mxe.exon1, mytx.ids.e1)
    mytx.ids.mxe.exon2 <- intersect(mytx.ids.e3, mytx.ids.e4)
    mytx.ids.mxe.exon2 <- intersect(mytx.ids.mxe.exon2, mytx.ids.e2)

    mytx.ids.mxe.exon1 <- setdiff(mytx.ids.mxe.exon1, mytx.ids.intron1)
    mytx.ids.mxe.exon2 <- setdiff(mytx.ids.mxe.exon2, mytx.ids.intron2)
    tx_ids <- list()
    tx_ids[["txn_mxe_exon1"]] <- mytx.ids.mxe.exon1
    tx_ids[["txn_mxe_exon2"]] <- mytx.ids.mxe.exon2
    return(tx_ids)
}
