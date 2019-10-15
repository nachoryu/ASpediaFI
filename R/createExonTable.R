#'Create an exon table
#'
#'This function creates an exon table
#'
#'@param gtf_exons A GRanges object containing exons in gtf
#'@param ids transcript ids
#'@return a data frame containing an exon table
#'@details This function is borrowed from the \code{maser} package.
#'@references Veiga, D. (2019). maser: Mapping Alternative Splicing Events
#'to pRoteins. R package version 1.2.0. https://github.com/DiogoVeiga/maser
#'@keywords internal
#'@importFrom dplyr filter
#'@noRd
createExonTable <- function(gtf_exons, ids){
    transcript_id <- NULL
    res <- filter(as.data.frame(gtf_exons), transcript_id %in%
                      ids)
    res.df <- res[, c("seqnames", "start", "end",
                      "strand", "exon_id", "transcript_name")]
    colnames(res.df) <- c("chromosome", "start",
                          "end", "strand", "exon", "transcript")
    return(res.df)
}
