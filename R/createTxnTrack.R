#'Create a transcript track
#'
#'This function creates a transcript track
#'
#'@param res.df A data frame containing an exon table
#'@param trackLabel Label of a transcript track
#'@param featureName Name of feature
#'@return a GeneRegionTrack object
#'@details This function is borrowed from the \code{maser} package.
#'@references Veiga, D. (2019). maser: Mapping Alternative Splicing Events
#'to pRoteins. R package version 1.2.0. https://github.com/DiogoVeiga/maser
#'@keywords internal
#'@importFrom Gviz GeneRegionTrack
createTxnTrack <- function(res.df, trackLabel, featureName){
    if (nrow(res.df) > 0) {
        res.df$feature <- featureName
        txnTrack <- GeneRegionTrack(range = res.df, name = trackLabel,
                                    transcriptAnnotation = "transcript",
                                    col = NULL, col.line = NULL)
    }
    else {
        txnTrack <- GeneRegionTrack(range = GRanges(), name = trackLabel,
                                    transcriptAnnotation = "transcript",
                                    col = NULL, col.line = NULL)
    }
    return(txnTrack)
}
