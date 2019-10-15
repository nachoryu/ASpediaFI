#'AS event annotation
#'
#'This function converts information of AS events into ASpedia format.
#'
#'@param ASlist a list of AS events
#'@param gtf a GRanges object representing GTF
#'@return A list of annotated AS events
#'@keywords internal
#'@importFrom limma strsplit2
#'@noRd
annotate <- function(ASlist, gtf){

    a5ss <- ASlist$A5SS
    a3ss <- ASlist$A3SS
    se <- ASlist$SE
    mxe <- ASlist$MXE
    ri <- ASlist$RI

    id_to_name <- data.frame(id = gtf$gene_id, name = gtf$gene_name,
                             stringsAsFactors = FALSE)
    id_to_name <- unique(id_to_name)
    rownames(id_to_name) <- id_to_name$id

    #Annotate A5SS events
    if(length(a5ss)){
        short.start <- strsplit2(a5ss[,"ShortEX"], split = "-")[,1]
        short.end <- strsplit2(a5ss[,"ShortEX"], split = "-")[,2]
        long.start <- strsplit2(a5ss[,"LongEX"], split = "-")[,1]
        long.end <- strsplit2(a5ss[,"LongEX"], split = "-")[,2]
        nei.start <- strsplit2(a5ss[,"NeighborEX"], split = "-")[,1]
        nei.end <- strsplit2(a5ss[,"NeighborEX"], split = "-")[,2]
        gene.symbol <- id_to_name[a5ss[,"EnsID"], "name"]

        aspedia.id <- ifelse(a5ss[,"Strand"] == "+",
                             paste(gene.symbol, "A5SS", a5ss[,"Nchr"],
                                   long.start, short.end,
                                   long.end, nei.start, nei.end, sep = ":"),
                             paste(gene.symbol, "A5SS",  a5ss[,"Nchr"],
                                   long.end, short.start,
                                   long.start, nei.end, nei.start, sep = ":"))

        a5ss <- rbind(cbind(a5ss, EventID = aspedia.id))
        rownames(a5ss) <- 1:nrow(a5ss)
    }

    #Annotate A3SS events
    if(length(a3ss)){
        short.start <- strsplit2(a3ss[,"ShortEX"], split = "-")[,1]
        short.end <- strsplit2(a3ss[,"ShortEX"], split = "-")[,2]
        long.start <- strsplit2(a3ss[,"LongEX"], split = "-")[,1]
        long.end <- strsplit2(a3ss[,"LongEX"], split = "-")[,2]
        nei.start <- strsplit2(a3ss[,"NeighborEX"], split = "-")[,1]
        nei.end <- strsplit2(a3ss[,"NeighborEX"], split = "-")[,2]
        gene.symbol <- id_to_name[a3ss[,"EnsID"], "name"]

        aspedia.id <- ifelse(a3ss[,"Strand"] == "+",
                             paste(gene.symbol, "A3SS", a3ss[,"Nchr"],
                                   nei.start, nei.end,
                                   long.start, short.start, long.end,
                                   sep = ":"),
                             paste(gene.symbol, "A3SS",
                                   a3ss[,"Nchr"], nei.end, nei.start, long.end,
                                   short.end, long.start, sep = ":"))

        a3ss <- rbind(cbind(a3ss, EventID = aspedia.id))
        rownames(a3ss) <- 1:nrow(a3ss)
    }

    #Annotate SE events
    if(length(se)){
        first.start <- strsplit2(se[,"1stEX"], split = "-")[,1]
        first.end <- strsplit2(se[,"1stEX"], split = "-")[,2]
        up.start <- strsplit2(se[,"UpEX"], split = "-")[,1]
        up.end <- strsplit2(se[,"UpEX"], split = "-")[,2]
        down.start <- strsplit2(se[,"DownEX"], split = "-")[,1]
        down.end <- strsplit2(se[,"DownEX"], split = "-")[,2]
        gene.symbol <- id_to_name[se[,"EnsID"], "name"]

        aspedia.id <- ifelse(se[,"Strand"] == "+",
                             paste(gene.symbol, "SE", se[,"Nchr"],
                                   down.start, down.end,
                                   first.start, first.end, up.start,
                                   up.end, sep = ":"),
                             paste(gene.symbol, "SE", se[,"Nchr"],
                                   up.end, up.start, first.end,
                                   first.start, down.end, down.start,
                                   sep = ":"))

        se <- rbind(cbind(se, EventID = aspedia.id))
        rownames(se) <- 1:nrow(se)
    }

    #Annotate MXE events
    if(length(mxe)){
        first.start <- strsplit2(mxe[,"1stEX"], split = "-")[,1]
        first.end <- strsplit2(mxe[,"1stEX"], split = "-")[,2]
        second.start <- strsplit2(mxe[,"2ndEX"], split = "-")[,1]
        second.end <- strsplit2(mxe[,"2ndEX"], split = "-")[,2]
        up.start <- strsplit2(mxe[,"UpEX"], split = "-")[,1]
        up.end <- strsplit2(mxe[,"UpEX"], split = "-")[,2]
        down.start <- strsplit2(mxe[,"DownEX"], split = "-")[,1]
        down.end <- strsplit2(mxe[,"DownEX"], split = "-")[,2]
        gene.symbol <- id_to_name[mxe[,"EnsID"], "name"]

        aspedia.id <- ifelse(mxe[,"Strand"] == "+",
                             paste(gene.symbol, "MXE", mxe[,"Nchr"],
                                   down.start, down.end,
                                   first.start, first.end, second.start,
                                   second.end, up.start, up.end, sep = ":"),
                             paste(gene.symbol, "MXE", mxe[,"Nchr"],
                                   up.end, up.start, second.end,
                                   second.start, first.end, first.start,
                                   down.end, down.start, sep = ":"))

        mxe <- rbind(cbind(mxe, EventID = aspedia.id))
        rownames(mxe) <- 1:nrow(mxe)
    }

    #Annotate RI events
    if(length(ri)){
        up.start <- strsplit2(ri[,"UpEX"], split = "-")[,1]
        up.end <- strsplit2(ri[,"UpEX"], split = "-")[,2]
        down.start <- strsplit2(ri[,"DownEX"], split = "-")[,1]
        down.end <- strsplit2(ri[,"DownEX"], split = "-")[,2]
        re.start <- strsplit2(ri[,"RetainEX"], split = "-")[,1]
        re.end <- strsplit2(ri[,"RetainEX"], split = "-")[,2]
        gene.symbol <- id_to_name[ri[,"EnsID"], "name"]

        aspedia.id <- ifelse(ri[,"Strand"] == "+",
                             paste(gene.symbol, "RI", ri[,"Nchr"],
                                   down.start, down.end,
                                   up.start, up.end, sep = ":"),
                             paste(gene.symbol, "RI", ri[,"Nchr"],
                                   up.end, up.start, down.end,
                                   down.start, sep = ":"))

        ri <- rbind(cbind(ri, EventID = aspedia.id))
        rownames(ri) <- 1:nrow(ri)
    }

    return(list(A5SS = a5ss, A3SS = a3ss, SE = se, MXE = mxe, RI = ri))
}
