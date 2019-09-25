#'Detection of retained introns
#'
#'This function detects retained intron (RI) events.
#'
#'@param exon.info A matrix containing exon locus
#'@param intron.info A matrix containing intron locus
#'@param alt.intron.info A matrix containing flanking introns of
#'alternative exons
#'@param tx.gene A table of transcripts including ID, gene names, etc.
#'@return A matrix containing retained intron events
#'@details This function is borrowed from the \code{IVAS} package.
#'@references Han, S. et al. (2017). Genome wide discovery of genetic variants
#'affecting alternative splicing patterns in human using bioinformatics method.
#'\emph{Genes & Genomics}, 39.
#'@keywords internal
detectRI <- function(exon.info, intron.info, alt.intron.info, tx.gene = NULL){
    RI.exons <- NULL
    alt.intron.info <- alt.intron.info[1,]
    exon.mat <- exon.info
    intron.mat <- intron.info
    RI.range <- alt.intron.info
    RIjunc1 <- sort(RI.range)[2] + 1
    RIjunc2 <- sort(RI.range)[1] - 1
    RI.exons <- rbind(unique(exon.info[exon.info[,"end"] > RIjunc1 &
                                           exon.info[,"start"] < RIjunc2,]))
    if (length(RI.exons) == 0) return (NULL)
    colnames(RI.exons) <- c("start", "end")
    get.RI.intron <- function(p.RI){
        ri.st.end <- as.integer(strsplit(p.RI,"-")[[1]])
        intron.info[intron.info[,"start"] - 1 > ri.st.end[1] &
                        intron.info[,"end"] < ri.st.end[2],]
    }
    RI.intron <- do.call(rbind,lapply(paste(RI.exons[,"start"],
                                            RI.exons[,"end"], sep = "-"),
                                      get.RI.intron))
    colnames(RI.intron) <- c("start","end")
    new.RIjunc1 <- unique(RI.intron[,2]) + 1
    new.RIjunc2 <- unique(RI.intron[,1]) - 1
    colnames(RI.exons) <- c("start","end")
    targetEX <- paste(RI.exons[,"start"], RI.exons[,"end"], sep = "-")
    downEX <- rbind(exon.mat[is.element(exon.mat[,"end"], new.RIjunc2),])
    upEX <- rbind(exon.mat[is.element(exon.mat[,"start"], new.RIjunc1),])
    downEX <- rbind(downEX[is.element(rownames(downEX),
                                      intersect(rownames(downEX),
                                                rownames(upEX))),])
    upEX <- rbind(upEX[is.element(rownames(upEX),
                                  intersect(rownames(downEX),
                                            rownames(upEX))),])
    if (nrow(downEX) == 0 | nrow(upEX) == 0)    return (NULL)
    p.downEX <- paste(downEX[,"start"], downEX[,"end"], sep = "-")
    p.upEX <- paste(upEX[,"start"], upEX[,"end"], sep = "-")
    RI.result <- NULL
    for (i in 1:length(p.upEX)){
        for (j in 1:length(p.downEX)){
            tested.tx <- tx.gene[tx.gene[tx.gene[,"TXID"] ==
                                             rownames(downEX)[j], "TXSTART"] <=
                                     downEX[j,"start"] &
                                     tx.gene[tx.gene[,"TXID"] ==
                                                 rownames(downEX)[j],"TXEND"] >=
                                     upEX[i,"end"],"TXID"]
            tested.down.num <- rownames(downEX) == rownames(downEX)[j]
            tested.down.ex <- rbind(downEX[tested.down.num,])
            rownames(tested.down.ex) <- rownames(downEX)[tested.down.num]
            tested.down.ex <- rbind(tested.down.ex[is.element(rownames(tested.down.ex),
                                                              tested.tx),])
            if (nrow(tested.down.ex) == 0)    next
            tested.down.ex <- rbind(tested.down.ex[tested.down.ex[,"end"] <
                                                       upEX[i,"start"],])
            if (nrow(tested.down.ex) == 0) next
            tested.down.ex <- rbind(tested.down.ex[order(tested.down.ex[,1],
                                                         decreasing=TRUE)[1],])
            if (upEX[i,"start"] <= tested.down.ex[,"end"]) next
            up.down.EX <- paste(tested.down.ex[,"start"], upEX[i,"end"],
                                sep = "-")
            doex.des <- unique(rbind(exon.info[exon.info[,"end"] ==
                                                   tested.down.ex[,"end"],]))
            upex.des <- unique(rbind(exon.info[exon.info[,"start"] ==
                                                   upEX[i,"start"],]))
            doex.des <- paste(sort(paste(doex.des[,"start"], doex.des[,"end"],
                                         sep = "-")), collapse = ",")
            upex.des <- paste(sort(paste(upex.des[,"start"], upex.des[,"end"],
                                         sep = "-")), collapse = ",")
            RIex.des <- paste(sort(unique(c(up.down.EX, targetEX))),
                              collapse = ",")
            RI.result <- rbind(RI.result, cbind(rbind(up.down.EX),
                                                paste(tested.down.ex,
                                                      collapse = "-"),
                                                p.upEX[i], RIex.des, doex.des,
                                                upex.des, "RI", "possible"))
        }
    }
    if (length(RI.result) == 0) return (NULL)
    colnames(RI.result) <- c("RetainEX", "DownEX", "UpEX", "Retain_des",
                             "Do_des", "Up_des", "Types", "status")
    spliced.int <- paste(as.double(do.call(rbind,
                                           strsplit(RI.result[,"DownEX"],
                                                    "-"))[,2]) + 1,
                         as.double(do.call(rbind,
                                           strsplit(RI.result[,"UpEX"],
                                                    "-"))[,1]) - 1, sep = "-")
    Re.int.test <- is.element(spliced.int, paste(intron.mat[,"start"],
                                                 intron.mat[,"end"], sep = "-"))
    Re.ex <- paste(exon.mat[,"start"], exon.mat[,"end"], sep = "-")
    Re.ex.test <- RI.result[,"RetainEX"] %in% Re.ex
    RI.result[Re.int.test & Re.ex.test, "status"] <- "exist"
    RI.result <- rbind(RI.result[,c("RetainEX", "DownEX", "UpEX", "Types",
                                    "status")])
    rownames(RI.result) <- 1:nrow(RI.result)
    return(unique(RI.result))
}
