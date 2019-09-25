#'Detection of AS events from TxDb
#'
#'This function detects AS events from a TxDb object.
#'
#'@param TxDb A TxDb object
#'@param num.cores The number of cores for parallel processing
#'@return A list of AS events
#'@details This function is modified from \code{Splicingfinder} in the
#'\code{IVAS} package.
#'@references Han, S. et al. (2017). Genome wide discovery of genetic variants
#'affecting alternative splicing patterns in human using bioinformatics method.
#'\emph{Genes & Genomics}, 39.
#'@keywords internal
#'@importFrom BiocParallel bplapply SnowParam
#'@importFrom IVAS findAlternative chrseparate
#'@importFrom GenomicFeatures exonsBy intronsByTranscript
#'@importFrom GenomicRanges start end seqnames
#'@importFrom biomaRt select
detectEventFromTxDb <- function(TxDb, num.cores = 1){
    CalAlt <- function(altInvalue){
        if (length(altInvalue) == 0)
            return (NULL)
        exonInfo <- altInvalue["exonRange"]
        intronInfo <- altInvalue["intronRange"]
        tx.gene <- altInvalue[["tableBygene"]]
        altintron <- altInvalue[["alterIntron"]]
        rownames(tx.gene) <- tx.gene[,"TXID"]
        predictedSQTL <- NULL
        exon.start <- unlist(start(exonInfo[[1]]))
        exon.end <- unlist(end(exonInfo[[1]]))
        exon.mat <- cbind(exon.start, exon.end)
        exon.locus <- paste(exon.start, exon.end, sep = "-")
        intron.start <- unlist(start(intronInfo[[1]]))
        intron.end <- unlist(end(intronInfo[[1]]))
        intron.mat <- cbind(intron.start, intron.end)
        intron.locus <- paste(intron.start, intron.end, sep = "-")
        alt.mat.int <- cbind(as.integer(start(altintron)),
                             as.integer(end(altintron)))
        names(exon.locus) <- names(exon.end)
        names(intron.locus) <- names(intron.end)
        strandinfo <- as.character(altInvalue[[1]]@strand@values)
        Nchr <- as.character(altInvalue[[1]]@seqnames@values)
        colnames(exon.mat) <- c("start", "end")
        colnames(intron.mat) <- c("start", "end")
        colnames(alt.mat.int) <- c("start", "end")
        uni.tx <- unique(rownames(exon.mat))
        final.exon.mat <- NULL
        for (i in 1:length(uni.tx)){
            each.tx.info <- rbind(exon.mat[rownames(exon.mat) == uni.tx[i],])
            each.tx.mat <- rbind(each.tx.info[order(each.tx.info[,"start"]),])
            if (nrow(each.tx.mat) == 1) rownames(each.tx.mat) <- uni.tx[i]
            final.exon.mat <- rbind(final.exon.mat, each.tx.mat)
        }
        final.intron.mat <- NULL
        for (i in 1:length(uni.tx)){
            each.tx.info <- rbind(intron.mat[rownames(intron.mat) == uni.tx[i],])
            each.tx.mat <- rbind(each.tx.info[order(each.tx.info[,"start"]),])
            if (nrow(each.tx.mat) == 1) rownames(each.tx.mat) <- uni.tx[i]
            final.intron.mat <- rbind(final.intron.mat, each.tx.mat)
        }
        final.result.list <- NULL
        total.list.names <- NULL
        result.n <- 1
        for(i in 1:nrow(alt.mat.int)){
            ES.result <- detectSEandMXE(final.exon.mat, final.intron.mat,
                                rbind(alt.mat.int[i,]))
            ES.result <- rbind(ES.result[ES.result[,"status"] == "exist",])
            if (length(ES.result) != 0){
                final.result.list[[result.n]] <- ES.result
                total.list.names <- c(total.list.names,
                                      paste("ES.",
                                            length(grep("ES",
                                                        total.list.names)) + 1,
                                            sep = ""))
                result.n <- result.n + 1
            }
            ASS.result <- detectASS(final.exon.mat, final.intron.mat,
                                    rbind(alt.mat.int[i,]))
            ASS.result <- rbind(ASS.result[ASS.result[,"status"] == "exist",])
            if(length(ASS.result) != 0){
                final.result.list[[result.n]] <- ASS.result
                total.list.names <- c(total.list.names,
                                      paste("ASS.",
                                            length(grep("ASS",
                                                        total.list.names)) + 1,
                                            sep = ""))
                result.n <- result.n + 1
            }
            RI.result <- detectRI(final.exon.mat, final.intron.mat,
                                  rbind(alt.mat.int[i,]), tx.gene=tx.gene)
            RI.result <- rbind(RI.result[RI.result[,"status"] == "exist",])
            if(length(RI.result) != 0){
                final.result.list[[result.n]] <- RI.result
                total.list.names <- c(total.list.names,
                                      paste("RI.",
                                            length(grep("RI",
                                                        total.list.names)) + 1,
                                            sep = ""))
                result.n <- result.n + 1
            }
        }
        names(final.result.list) <- total.list.names
        cluster.result <- list(NULL, NULL, NULL)
        names(cluster.result) <- c("ES", "ASS", "RI")
        if (!is.null(final.result.list)){
            for (each.list.num in 1:length(final.result.list)){
                result.type <- unlist(strsplit(names(final.result.list[each.list.num]),
                                               "[.]"))[1]
                cluster.result[[result.type]] <- unique(rbind(cluster.result[[result.type]],
                                                              final.result.list[[each.list.num]]))
            }
        }
        if (length(cluster.result[["ES"]]) != 0){
            row.names(cluster.result[["ES"]]) <- 1:nrow(rbind(cluster.result[["ES"]]))
            cluster.result[["ES"]] <- cbind(unique(tx.gene[,"GENEID"]),
                                                  unique(tx.gene[,"TXCHROM"]),
                                                  unique(tx.gene[,"TXSTRAND"]),
                                                  cluster.result[["ES"]])
            colnames(cluster.result[["ES"]])[1:3] <- c("EnsID", "Nchr",
                                                             "Strand")
        }
        if (length(cluster.result[["ASS"]]) != 0){
            row.names(cluster.result[["ASS"]]) <- 1:nrow(rbind(cluster.result[["ASS"]]))
            cluster.result[["ASS"]] <- cbind(unique(tx.gene[,"GENEID"]),
                                             unique(tx.gene[,"TXCHROM"]),
                                             unique(tx.gene[,"TXSTRAND"]),
                                             cluster.result[["ASS"]])
            colnames(cluster.result[["ASS"]])[1:3] <- c("EnsID", "Nchr",
                                                        "Strand")
        }
        if (length(cluster.result[["RI"]]) != 0){
            row.names(cluster.result[["RI"]]) <- 1:nrow(rbind(cluster.result[["RI"]]))
            cluster.result[["RI"]] <- cbind(unique(tx.gene[,"GENEID"]),
                                            unique(tx.gene[,"TXCHROM"]),
                                            unique(tx.gene[,"TXSTRAND"]),
                                            cluster.result[["RI"]])
            colnames(cluster.result[["RI"]])[1:3] <- c("EnsID", "Nchr",
                                                       "Strand")
        }
        return(cluster.result)
    }
    mulspl <- function(j){
        strandinfo <- unique(each.chr.txTable[each.chr.txTable[,"GENEID"] ==
                                                  up2gene[j], "TXSTRAND"])
        Altvalue <- findAlternative(up2gene[j], each.chr.txTable,
                                    each.chr.exon.range, each.chr.intron.range,
                                    Total.chr[i])
        Alt.result <- CalAlt(Altvalue)
        Alt.result
    }
    trans.exon.range <- exonsBy(TxDb, by="tx")
    trans.intron.range <- intronsByTranscript(TxDb)
    txTable <- select(TxDb, keys=names(trans.exon.range),
                      columns=c("TXCHROM", "TXNAME", "GENEID", "TXSTART",
                                "TXEND", "TXSTRAND"), keytype="TXID")
    txTable <- gsub(" ", "", as.matrix(txTable))
    Total.chr <- as.character(sort(unique(txTable[,"TXCHROM"])))
    ES.finl.result <- NULL
    ASS.finl.result <- NULL
    IR.finl.result <- NULL
    Alt.result <- NULL
    j = NULL
    MP <- SnowParam(workers = num.cores, type="SOCK")
    for (i in 1:length(Total.chr)){
        print (paste("-------------------Processing : ", Total.chr[i],
                     " -------------------", sep = ""))
        each.chr.db <- chrseparate(TxDb, Total.chr[i])
        each.chr.names <- names(unlist(trans.exon.range))[as.character(seqnames(unlist(trans.exon.range))) == Total.chr[i]]
        each.chr.exon.range <- trans.exon.range[is.element(names(trans.exon.range),
                                                           each.chr.names),]
        each.chr.intron.range <- trans.intron.range[is.element(names(trans.intron.range),
                                                               each.chr.names),]
        transnum <- names(each.chr.exon.range)
        each.chr.txTable <- txTable[txTable[,"TXCHROM"] == Total.chr[i],]
        tx.num <- table(each.chr.txTable[,"GENEID"])
        not.one.gene <- which(tx.num > 1)
        up2gene <- names(not.one.gene)
        pa.result <- bplapply(seq_len(length(up2gene)), mulspl, BPPARAM = MP)
        ES.finl.result <- rbind(ES.finl.result, do.call(rbind,
                                                        lapply(pa.result,
                                                               function(x) x$"ES")))
        ASS.finl.result <- rbind(ASS.finl.result,do.call(rbind,
                                                         lapply(pa.result,
                                                                function(x) x$"ASS")))
        IR.finl.result <- rbind(IR.finl.result, do.call(rbind,
                                                        lapply(pa.result,
                                                               function(x) x$"RI")))
    }
    SE.final <- rbind(ES.finl.result[ES.finl.result[,"Types"] == "SE",])
    if(length(SE.final)){
        SE.final <- cbind(Index = paste0("SE", 1:nrow(SE.final)), SE.final)
    }
    MXE.final <- rbind(ES.finl.result[ES.finl.result[,"Types"] == "MXE",])
    if(length(MXE.final)){
        MXE.final <- cbind(Index = paste0("MXE", 1:nrow(MXE.final)), MXE.final)
    }
    A5SS.final <- rbind(ASS.finl.result[ASS.finl.result[,"Types"] == "A5SS",])
    if(length(A5SS.final)){
        A5SS.final <- cbind(Index = paste0("A5SS", 1:nrow(A5SS.final)),
                            A5SS.final)
    }
    A3SS.final <- rbind(ASS.finl.result[ASS.finl.result[,"Types"] == "A3SS",])
    if(length(A3SS.final)){
        A3SS.final <- cbind(Index = paste0("A3SS", 1:nrow(A3SS.final)),
                            A3SS.final)
    }
    RI.final <- IR.finl.result
    if(length(RI.final)){
        RI.final <- cbind(Index = paste0("RI", 1:nrow(RI.final)), RI.final)
    }
    return(list(A3SS = A3SS.final, A5SS = A5SS.final, SE = SE.final,
                MXE = MXE.final, RI = RI.final))
}

