#' AS event visualization
#'
#' This function visualizes an AS event.
#'
#' @param event AS event ID in ASpedia format
#' @param gtf a GRanges object representing GTF
#' @param psi a SummarizedExperiment object containing PSI of AS events and
#' groups of samples
#' @param zoom a logical to determine if genomic coordinates are zoomed
#' @return a Gviz object
#' @details This function is modified from the \code{plotTranscripts} function
#' in the \code{maser} package.
#' @references Veiga, D. (2019). maser: Mapping Alternative Splicing Events
#' to pRoteins. R package version 1.2.0. https://github.com/DiogoVeiga/maser
#' @keywords internal
#' @import ggplot2
#' @import grid
#' @import SummarizedExperiment
#' @importFrom limma strsplit2
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom Gviz plotTracks DataTrack
#' @importFrom scales hue_pal
#' @noRd
visualizeEvent <- function(event, gtf, psi = NULL, zoom = NULL) {
    if (is.null(zoom)) {
        zoom <- FALSE
    }

    grlevent <- makeGRangesListFromEvent(event, gtf)
    eventGr <- grlevent$grl
    type <- grlevent$type
    gtf_exons <- gtf[gtf$type == "exon", ]

    eventTrack <- createAnnotationTrack_event(eventGr, type)
    txn_tracks <- createAnnotationTrack_transcripts(eventGr,
                                                        gtf_exons, type)

    trackList <- list(eventTrack, txn_tracks$inclusionTrack,
                        txn_tracks$skippingTrack)

    if (!is.null(psi) & "condition" %in% colnames(colData(psi)) &
            event %in% rownames(psi)) {
        event.gene <- strsplit2(event, split = ":")[1]
        event.type <- strsplit2(event, split = ":")[2]
        psidat <- data.frame(PSI = assays(psi)[[1]][event, ],
                                Condition = colData(psi)$condition)
        grid.newpage()
        pushViewport(viewport(height = 0.4, y = 1, just = "top"))
        bp <- ggplotGrob(ggplot(data = psidat,
                                    aes_string(x = "Condition", y = "PSI",
                                                fill = "Condition")) +
                                stat_boxplot(geom = "errorbar") +
                                geom_boxplot() + geom_jitter() +
                                ggtitle(paste(event.gene, event.type)) +
                                xlab(NULL) + theme_bw() +
                                theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
            legend.key.height = unit(0.1, "npc"),
                legend.key.width = unit(0.05, "npc"),
                legend.title.align = 0.5,
                plot.title = element_text(hjust = 0.5)))
        grid.draw(bp)
        popViewport(1)

        pushViewport(viewport(height = 0.6, y = 0, just = "bottom"))
        if (zoom) {
            plotTracks(trackList, Inclusion = "red", Skipping = "blue",
                        Retention = "red", Non_Retention = "blue",
                        MXE_Exon1 = "red", MXE_Exon2 = "blue",
                        A5SS_Short = "blue", A5SS_Long = "red",
                        A3SS_Short = "blue", A3SS_Long = "red",
                        from = start(range(unlist(eventGr))) - 500,
                        to = end(range(unlist(eventGr))) +
                500, add = TRUE)
        } else {
            plotTracks(trackList, Inclusion = "red", Skipping = "blue",
                        Retention = "red", Non_Retention = "blue",
                        MXE_Exon1 = "red", MXE_Exon2 = "blue",
                        A5SS_Short = "blue", A5SS_Long = "red",
                        A3SS_Short = "blue", A3SS_Long = "red", add = TRUE)
        }
    } else {
        if (zoom) {
            plotTracks(trackList, Inclusion = "red", Skipping = "blue",
                        Retention = "red", Non_Retention = "blue",
                        MXE_Exon1 = "red", MXE_Exon2 = "blue",
                        A5SS_Short = "blue", A5SS_Long = "red",
                        A3SS_Short = "blue", A3SS_Long = "red",
                        from = start(range(unlist(eventGr))) - 500,
                        to = end(range(unlist(eventGr))) +
                500)
        } else {
            plotTracks(trackList, Inclusion = "red", Skipping = "blue",
                        Retention = "red", Non_Retention = "blue",
                        MXE_Exon1 = "red", MXE_Exon2 = "blue",
                        A5SS_Short = "blue", A5SS_Long = "red",
                        A3SS_Short = "blue", A3SS_Long = "red")
        }
    }
}
