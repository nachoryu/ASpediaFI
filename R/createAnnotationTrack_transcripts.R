#'Create transcript tracks from AS event
#'
#'This function creates transcript tracks from AS event.
#'
#'@param eventGr A GRangesList object of AS event
#'@param gtf_exons A GRanges object containing exons in gtf
#'@param type AS event type
#'@return a list of GeneRegionTrack objects
#'@details This function is borrowed from the \code{maser} package.
#'@references Veiga, D. (2019). maser: Mapping Alternative Splicing Events
#'to pRoteins. R package version 1.2.0. https://github.com/DiogoVeiga/maser
#'@keywords internal
#'@noRd
createAnnotationTrack_transcripts <- function(eventGr, gtf_exons, type){
    if(type == "A3SS"){
        if (as.character(strand(eventGr[1])) == "+") {
            tx_ids <- mapTranscriptsA3SSevent(eventGr, gtf_exons)
        }
        else {
            tx_ids <- mapTranscriptsA5SSevent(eventGr, gtf_exons)
        }
        res.df <- createExonTable(gtf_exons, tx_ids$txn_short)
        skippingTrack <- createTxnTrack(res.df, "A3SS Short",
                                         "A3SS_Short")
        res.df <- createExonTable(gtf_exons, tx_ids$txn_long)
        inclusionTrack <- createTxnTrack(res.df, "A3SS Long",
                                        "A3SS_Long")
        txn_tracks <- list(inclusionTrack = inclusionTrack,
                           skippingTrack = skippingTrack)
        return(txn_tracks)
    }
    if(type == "A5SS"){
        if (as.character(strand(eventGr[1])) == "+") {
            tx_ids <- mapTranscriptsA5SSevent(eventGr, gtf_exons)
        }
        else {
            tx_ids <- mapTranscriptsA3SSevent(eventGr, gtf_exons)
        }
        res.df <- createExonTable(gtf_exons, tx_ids$txn_short)
        skippingTrack <- createTxnTrack(res.df, "A5SS Short",
                                         "A5SS_Short")
        res.df <- createExonTable(gtf_exons, tx_ids$txn_long)
        inclusionTrack <- createTxnTrack(res.df, "A5SS Long",
                                        "A5SS_Long")
        txn_tracks <- list(inclusionTrack = inclusionTrack,
                           skippingTrack = skippingTrack)
        return(txn_tracks)
    }
    if(type == "SE"){
        tx_ids <- mapTranscriptsSEevent(eventGr, gtf_exons)
        res.df <- createExonTable(gtf_exons, tx_ids$txn_3exons)
        inclusionTrack <- createTxnTrack(res.df, "Inclusion",
                                         "Inclusion")
        res.df <- createExonTable(gtf_exons, tx_ids$txn_2exons)
        skippingTrack <- createTxnTrack(res.df, "Skipping",
                                        "Skipping")
        txn_tracks <- list(inclusionTrack = inclusionTrack,
                           skippingTrack = skippingTrack)
        return(txn_tracks)
    }
    if(type == "MXE"){
        tx_ids <- mapTranscriptsMXEevent(eventGr, gtf_exons)
        res.df <- createExonTable(gtf_exons, tx_ids$txn_mxe_exon1)
        inclusionTrack <- createTxnTrack(res.df, "MXE Exon 1",
                                         "MXE_Exon1")
        res.df <- createExonTable(gtf_exons, tx_ids$txn_mxe_exon2)
        skippingTrack <- createTxnTrack(res.df, "MXE Exon 2",
                                        "MXE_Exon2")
        txn_tracks <- list(inclusionTrack = inclusionTrack,
                           skippingTrack = skippingTrack)
        return(txn_tracks)
    }
    if(type == "RI"){
        tx_ids <- mapTranscriptsRIevent(eventGr, gtf_exons)
        res.df <- createExonTable(gtf_exons, tx_ids$txn_retention)
        retention_Track <- createTxnTrack(res.df, "Retention",
                                          "Retention")
        res.df <- createExonTable(gtf_exons, tx_ids$txn_nonRetention)
        nonRetention_Track <- createTxnTrack(res.df, "Non-retention",
                                             "Non_Retention")
        txn_tracks <- list(inclusionTrack = retention_Track,
                           skippingTrack = nonRetention_Track)
    }
}
