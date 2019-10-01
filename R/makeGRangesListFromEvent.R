#'Make a GRangesList object from an AS event
#'
#'This function makes a GRangesList object from an AS event.
#'
#'@param event AS event ID in ASpedia format
#'@return a GRangesList object
#'@keywords internal
#'@importFrom GenomicRanges makeGRangesListFromDataFrame
makeGRangesListFromEvent <- function(event, gtf){

    event <- unlist(strsplit(event, split = ":"))
    genesymbol <- event[1]
    geneid <- gtf$gene_id[grep(genesymbol, gtf$gene_name)[1]]
    type <- event[2]
    seqname <- event[3]
    strand <- ifelse(as.numeric(event[4]) < as.numeric(event[length(event)]),
                     "+", "-")
    coords <- as.numeric(event[4:length(event)])

    if(strand == "-"){
        coords <- rev(coords)
    }

    if((type == "A5SS" & strand == "+") | (type == "A3SS" & strand == "-")){
        df <- data.frame(seqname = seqname, start = coords[c(1, 1, 4)],
                         end = coords[c(2, 3, 5)], strand = strand,
                         exon = c("exon_short", "exon_long", "exon_flanking"),
                         GeneID = geneid, geneSymbol = genesymbol,
                         stringsAsFactors = FALSE)
        grl <- makeGRangesListFromDataFrame(df, split.field = "exon",
                                            keep.extra.columns = TRUE)
        grl <- grl[c("exon_long", "exon_short", "exon_flanking")]
    }

    if((type == "A3SS" & strand == "+") | (type == "A5SS" & strand == "-")){
        df <- data.frame(seqname = seqname, start = coords[c(1, 3, 4)],
                         end = coords[c(2, 5, 5)], strand = strand,
                         exon = c("exon_flanking", "exon_long", "exon_short"),
                         GeneID = geneid, geneSymbol = genesymbol,
                         stringsAsFactors = FALSE)
        grl <- makeGRangesListFromDataFrame(df, split.field = "exon",
                                            keep.extra.columns = TRUE)
        grl <- grl[c("exon_long", "exon_short", "exon_flanking")]
    }

    if(type == "SE"){
        df <- data.frame(seqname = seqname, start = coords[c(3, 1, 5)],
                         end = coords[c(4, 2, 6)], strand = strand,
                         exon = c("exon_target", "exon_upstream",
                                  "exon_downstream"),
                         GeneID = geneid, geneSymbol = genesymbol,
                         stringsAsFactors = FALSE)
        grl <- makeGRangesListFromDataFrame(df, split.field = "exon",
                                            keep.extra.columns = TRUE)
        grl <- grl[c("exon_target", "exon_upstream", "exon_downstream")]
    }

    if(type == "MXE"){
        df <- data.frame(seqname = seqname, start = coords[c(1, 3, 5, 7)],
                         end = coords[c(2, 4, 6, 8)], strand = strand,
                         exon = c("exon_upstream", "exon_1", "exon_2",
                                  "exon_downstream"),
                         GeneID = geneid, geneSymbol = genesymbol,
                         stringsAsFactors = FALSE)
        grl <- makeGRangesListFromDataFrame(df, split.field = "exon",
                                            keep.extra.columns = TRUE)
        grl <- grl[c("exon_1", "exon_2", "exon_upstream", "exon_downstream")]
    }

    if(type == "RI"){
        df <- data.frame(seqname = seqname, start = coords[c(1, 1, 3)],
                         end = coords[c(4, 2, 4)], strand = strand,
                         exon = c("exon_ir", "exon_upstream",
                                  "exon_downstream"),
                         GeneID = geneid, geneSymbol = genesymbol,
                         stringsAsFactors = FALSE)
        grl <- makeGRangesListFromDataFrame(df, split.field = "exon",
                                            keep.extra.columns = TRUE)
        grl <- grl[c("exon_ir", "exon_upstream", "exon_downstream")]
    }

    return(list(grl = grl, type = type))

}
