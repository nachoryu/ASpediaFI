#'Create an event track from AS event
#'
#'This function creates an event track from AS event.
#'
#'@param eventGr A GRangesList object of AS event
#'@param type AS event type
#'@return an AnnotationTrack object
#'@details This function is borrowed from the \code{maser} package.
#'@references Veiga, D. (2019). maser: Mapping Alternative Splicing Events
#'to pRoteins. R package version 1.2.0. https://github.com/DiogoVeiga/maser
#'@keywords internal
#'@importFrom Gviz AnnotationTrack feature
#'@noRd
createAnnotationTrack_event <- function(eventGr, type){
    if(type == "A3SS"){
        trackGr <- c(eventGr$exon_flanking, eventGr$exon_long,
                     eventGr$exon_flanking,
                     eventGr$exon_short)
        trackGr$group <- rep(c("A3SS Long", "A3SS Short"),
                             c(2, 2))
        trackGr$type <- rep("A3SS", 4)
        event_track <- AnnotationTrack(trackGr, name = "Event",
                                       groupAnnotation = "group",
                                       shape = "box", stacking = "squish",
                                       id = "A3SS", col = NULL,
                                       col.line = NULL)
        Gviz::feature(event_track) <- rep(c("A3SS_Long", "A3SS_Short"),
                                          c(2, 2))
        return(event_track)
    }
    if(type == "A5SS"){
        trackGr <- c(eventGr$exon_long, eventGr$exon_flanking,
                     eventGr$exon_short,
                     eventGr$exon_flanking)
        trackGr$group <- rep(c("A5SS Long", "A5SS Short"),
                             c(2, 2))
        trackGr$type <- rep("A5SS", 4)
        event_track <- AnnotationTrack(trackGr, name = "Event",
                                       groupAnnotation = "group", shape = "box",
                                       stacking = "squish", id = "A5SS",
                                       col = NULL, col.line = NULL)
        Gviz::feature(event_track) <- rep(c("A5SS_Long", "A5SS_Short"),
                                          c(2, 2))
        return(event_track)
    }
    if(type == "SE"){
        transcript_id <- NULL
        trackGr <- c(unlist(eventGr), unlist(eventGr[2:3]))
        trackGr$group <- rep(c("Inclusion", "Skipping"),
                             c(3, 2))
        trackGr$type <- rep("Exon skipping", 5)
        event_track <- AnnotationTrack(trackGr, name = "Event",
                                       groupAnnotation = "group", shape = "box",
                                       stacking = "squish",
                                       id = "Exon skipping",
                                       col = NULL, col.line = NULL)
        Gviz::feature(event_track) <- rep(c("Inclusion", "Skipping"),
                                          c(3, 2))
        return(event_track)
    }
    if(type == "MXE"){
        trackGr <- c(eventGr$exon_upstream, eventGr$exon_1,
                     eventGr$exon_downstream,
                     eventGr$exon_upstream, eventGr$exon_2,
                     eventGr$exon_downstream)
        trackGr$group <- rep(c("MXE_Exon1", "MXE_Exon2"),
                             c(3, 3))
        trackGr$type <- rep("Mutually Exclusive Exons", 3)
        event_track <- AnnotationTrack(trackGr, name = "Event",
                                       groupAnnotation = "group", shape = "box",
                                       stacking = "squish",
                                       id = "Mutually Exclusive Exons",
                                       col = NULL, col.line = NULL)
        Gviz::feature(event_track) <- rep(c("MXE_Exon1", "MXE_Exon2"),
                                          c(3, 3))
        return(event_track)
    }
    if(type == "RI"){
        transcript_id <- NULL
        trackGr <- c(eventGr$exon_ir, eventGr$exon_upstream,
                     eventGr$exon_downstream)
        trackGr$group <- rep(c("Retention", "Non-retention"),
                             c(1, 2))
        trackGr$type <- rep("Intron retention", 3)
        event_track <- AnnotationTrack(trackGr, name = "Event",
                                       groupAnnotation = "group", shape = "box",
                                       stacking = "squish",
                                       id = "Intron retention",
                                       col = NULL, col.line = NULL)
        Gviz::feature(event_track) <- rep(c("Retention", "Non_Retention"),
                                          c(1, 2))
        return(event_track)
    }
}
