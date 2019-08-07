#'Detection of AS events from TxDb
#'
#'This function detects AS events from a TxDb object.
#'
#'@param ASlist a list of AS events
#'@param TxDb A TxDb object
#'@return A list of AS events
#'@details This function is modified from \code{Splicingfinder} in the
#'\code{IVAS} package.
#'@references Han, S. et al. (2017). Genome wide discovery of genetic variants
#'affecting alternative splicing patterns in human using bioinformatics method.
#'\emph{Genes & Genomics}, 39.
#'@keywords internal
#'@importFrom IVAS findAlternative chrseparate
#'@importFrom GenomicFeatures exonsBy intronsByTranscript
#'@importFrom biomaRt select
#'@importFrom limma strsplit2
filterEvent <- function(ASlist, TxDb){

    ###Extract information relevant to AS genes
    se <- rbind(ASlist[["SE"]])
    mxe <- rbind(ASlist[["MXE"]])
    a5ss <- rbind(ASlist[["A5SS"]])
    a3ss <- rbind(ASlist[["A3SS"]])
    ri <- rbind(ASlist[["RI"]])

    tested.geneid <- unique(c(se[,"EnsID"], mxe[,"EnsID"], a5ss[,"EnsID"],
                              a3ss[,"EnsID"], ri[,"EnsID"]))
    total.chr <- unique(c(se[,"Nchr"], mxe[,"Nchr"], a5ss[,"Nchr"],
                          a3ss[,"Nchr"], ri[,"Nchr"]))
    TxDb <- chrseparate(TxDb, total.chr)
    total.exon.range <- exonsBy(TxDb, by = "tx")
    total.intron.range <- intronsByTranscript(TxDb)
    txTable <- select(TxDb, keys = names(total.exon.range),
                      columns = c("TXCHROM", "TXNAME", "GENEID",
                                  "TXSTART", "TXEND", "TXSTRAND"),
                      keytype = "TXID")
    sub.txTable <- rbind(txTable[is.element(txTable$GENEID, tested.geneid),])
    sub.exon.range <- total.exon.range[is.element(names(total.exon.range),
                                                  as.matrix(sub.txTable[, "TXID"])),]
    sub.intron.range <- total.intron.range[is.element(names(total.intron.range),
                                                      as.matrix(sub.txTable[, "TXID"])),]

    #Valdiation for exon skipping
    se <- rbind(se[se[,"1stEX"] != "NA" & se[,"2ndEX"] == "NA",])
    se.genes <- unique(se[,"EnsID"])
    se.sane <- lapply(1:length(se.genes), function(i){
        each.se <- rbind(se[se[,"EnsID"] == se.genes[i],])
        chr.se <- unique(each.se[,"Nchr"])
        each.se.sane <- NULL
        for(chr in chr.se){
            each.se.sub <- rbind(each.se[each.se[,"Nchr"] == chr,])
            each.tx.info <- sub.txTable[sub.txTable$GENEID ==
                                            unique(each.se.sub[,"EnsID"]),]
            each.tx.info <- each.tx.info[each.tx.info$TXCHROM == chr,]
            Exon.info <- findAlternative(unique(each.se.sub[,"EnsID"]),
                                         sub.txTable, sub.exon.range,
                                         sub.intron.range, chr)$exonRange
            exon.start <- unlist(start(Exon.info))
            exon.end <- unlist(end(Exon.info))
            exon.mat <- cbind(exon.start, exon.end)
            each.exon.locus <- paste(exon.start, exon.end, sep = "-")
            names(each.exon.locus) <- names(exon.end)
            sane <- sapply(1:length(each.se.sub[,"1stEX"]), function(each.nums){
                index.num <- each.se.sub[each.nums, "Index"]
                each.targets <- each.se.sub[each.nums, "1stEX"]
                each.se.result <- rbind(each.se.sub[each.nums, ])
                s.tar.ex <- each.se.result[,"1stEX"]
                s.do.ex <- each.se.result[,"DownEX"]
                s.up.ex <- each.se.result[,"UpEX"]
                test.do.tx <- names(each.exon.locus[is.element(each.exon.locus,
                                                               s.do.ex) | is.element(each.exon.locus, s.do.ex)])
                test.up.tx <- names(each.exon.locus[is.element(each.exon.locus,
                                                               s.up.ex) | is.element(each.exon.locus, s.up.ex)])
                test.tx.ex <- each.exon.locus[is.element(names(each.exon.locus),
                                                         intersect(test.do.tx,
                                                                   test.up.tx))]
                included.tx <- unique(names(test.tx.ex)[is.element(test.tx.ex,
                                                                   s.tar.ex)])
                skipped.tx <- unique(names(test.tx.ex)[!is.element(names(test.tx.ex),
                                                                   included.tx)])
                included.tx <- each.tx.info[is.element(as.matrix(each.tx.info["TXID"]),
                                                       included.tx), "TXNAME"]
                skipped.tx <- each.tx.info[is.element(as.matrix(each.tx.info["TXID"]),
                                                      skipped.tx), "TXNAME"]
                length(included.tx) != 0 & length(skipped.tx) != 0
            })
            each.se.sane <- rbind(each.se.sane, rbind(each.se.sub[sane,]))
        }
        each.se.sane
    }
    )
    se.sane <- do.call(rbind, se.sane)
    se.sane <- rbind(se.sane[!duplicated(paste(se.sane[,"Nchr"],
                                               se.sane[,"1stEX"],
                                         se.sane[,"DownEX"], se.sane[,"UpEX"],
                                         sep = ":")),])
    se.sane <- rbind(se.sane[,c("EnsID", "Nchr", "Strand", "1stEX", "DownEX",
                                "UpEX")])

    #Validation for mutually exclusive exon
    mxe <- rbind(mxe[mxe[,"1stEX"] != "NA" & mxe[,"2ndEX"] != "NA",])
    mxe.genes <- unique(mxe[,"EnsID"])
    mxe.sane <- lapply(1:length(mxe.genes), function(i){
        each.mxe <- rbind(mxe[mxe[,"EnsID"] == mxe.genes[i],])
        chr.mxe <- unique(each.mxe[,"Nchr"])
        each.mxe.sane <- NULL
        for(chr in chr.mxe){
            each.mxe.sub <- rbind(each.mxe[each.mxe[,"Nchr"] == chr,])
            each.tx.info <- sub.txTable[sub.txTable$GENEID ==
                                            unique(each.mxe.sub[,"EnsID"]),]
            each.tx.info <- each.tx.info[each.tx.info$TXCHROM == chr,]
            Exon.info <- findAlternative(unique(each.mxe.sub[,"EnsID"]),
                                         sub.txTable, sub.exon.range,
                                         sub.intron.range, chr)$exonRange
            exon.start <- unlist(start(Exon.info))
            exon.end <- unlist(end(Exon.info))
            exon.mat <- cbind(exon.start, exon.end)
            each.exon.locus <- paste(exon.start, exon.end, sep = "-")
            names(each.exon.locus) <- names(exon.end)
            sane <- sapply(1:length(each.mxe.sub[,"1stEX"]), function(each.nums){
                index.num <- each.mxe.sub[each.nums, "Index"]
                each.targets <- each.mxe.sub[each.nums, "1stEX"]
                each.mxe.result <- rbind(each.mxe.sub[each.nums,])
                s.first.ex <- each.mxe.result[,"1stEX"]
                s.second.ex <- each.mxe.result[,"2ndEX"]
                s.do.ex <- each.mxe.result[,"DownEX"]
                s.up.ex <- each.mxe.result[,"UpEX"]
                test.do.tx <- names(each.exon.locus[is.element(each.exon.locus,
                                                               s.do.ex)])
                test.up.tx <- names(each.exon.locus[is.element(each.exon.locus,
                                                               s.up.ex)])
                test.tx.ex <- each.exon.locus[is.element(names(each.exon.locus),
                                                         intersect(test.do.tx,
                                                                   test.up.tx))]
                first.tx <- unique(names(test.tx.ex)[is.element(test.tx.ex,
                                                                s.first.ex)])
                second.tx <- unique(names(test.tx.ex)[is.element(test.tx.ex,
                                                                 s.second.ex)])
                included.tx <- first.tx[!is.element(first.tx, second.tx)]
                skipped.tx <- second.tx[!is.element(second.tx, first.tx)]
                included.tx <- each.tx.info[is.element(as.matrix(each.tx.info["TXID"]),
                                                       included.tx), "TXNAME"]
                skipped.tx <- each.tx.info[is.element(as.matrix(each.tx.info["TXID"]),
                                                      skipped.tx), "TXNAME"]
                length(included.tx) != 0 & length(skipped.tx) != 0
            })
            each.mxe.sane <- rbind(each.mxe.sane, rbind(each.mxe.sub[sane,]))

        }
        each.mxe.sane
    }
    )
    mxe.sane <- do.call(rbind, mxe.sane)
    mxe.sane <- rbind(mxe.sane[!duplicated(paste(mxe.sane[,"Nchr"],
                                                 mxe.sane[,"1stEX"],
                                                 mxe.sane[,"2ndEX"],
                                                 mxe.sane[,"UpEX"],
                                                 mxe.sane[,"DownEX"],
                                                 sep = ":")),])
    mxe.sane <- rbind(mxe.sane[,c("EnsID", "Nchr", "Strand", "1stEX", "2ndEX",
                            "DownEX", "UpEX")])

    #Validation for ASS
    ass <- rbind(a5ss, a3ss)
    nomatch.start <- strsplit2(ass[,"ShortEX"], split = "-")[,1] !=
        strsplit2(ass[,"LongEX"], split = "-")[,1]
    nomatch.end <- strsplit2(ass[,"ShortEX"], split = "-")[,2] !=
        strsplit2(ass[,"LongEX"], split = "-")[,2]
    violated <- nomatch.start & nomatch.end
    ass <- rbind(ass[!violated,])
    ass.genes <- unique(ass[,"EnsID"])
    ass.sane <- lapply(1:length(ass.genes), function(i){
        each.ass <- rbind(ass[ass[,"EnsID"] == ass.genes[i],])
        chr.ass <- unique(each.ass[,"Nchr"])
        each.ass.sane <- NULL
        for(chr in chr.ass){
            each.ass.sub <- rbind(each.ass[each.ass[,"Nchr"] == chr,])
            each.tx.info <- sub.txTable[sub.txTable$GENEID ==
                                            unique(each.ass.sub[,"EnsID"]),]
            each.tx.info <- each.tx.info[each.tx.info$TXCHROM == chr,]
            Exon.info <- findAlternative(unique(each.ass.sub[,"EnsID"]),
                                         sub.txTable, sub.exon.range,
                                         sub.intron.range, chr)$exonRange
            exon.start <- unlist(start(Exon.info))
            exon.end <- unlist(end(Exon.info))
            exon.mat <- cbind(exon.start, exon.end)
            each.exon.locus <- paste(exon.start, exon.end, sep = "-")
            names(each.exon.locus) <- names(exon.end)
            each.ass.sane.sub <- lapply(1:length(each.ass.sub[,"ShortEX"]), function(each.nums) {
                short.ok <- FALSE
                long.ok <- FALSE
                index.num <- each.ass.sub[each.nums, "Index"]
                each.targets <- each.ass.sub[each.nums, "ShortEX"]
                each.ass.result <- rbind(each.ass.sub[each.nums, ])
                s.short.ex <- each.ass.result[,"ShortEX"]
                s.long.ex <- each.ass.result[,"LongEX"]
                short.nei.ex <- each.ass.result[,"ShortNeighborEX"]
                long.nei.ex <- each.ass.result[,"LongNeighborEX"]
                test.snei.tx <- c(intersect(names(each.exon.locus[(is.element(each.exon.locus, s.short.ex))]),
                                            names(each.exon.locus[(is.element(each.exon.locus, short.nei.ex))])),
                                  intersect(names(each.exon.locus[(is.element(each.exon.locus, s.long.ex))]),
                                            names(each.exon.locus[(is.element(each.exon.locus, short.nei.ex))])))
                test.snei.ex <- each.exon.locus[is.element(names(each.exon.locus), test.snei.tx)]
                snei.shortex <- unique(intersect(names(test.snei.ex)[is.element(test.snei.ex, s.short.ex)],
                                                 names(test.snei.ex)[is.element(test.snei.ex, short.nei.ex)]))
                snei.longex <- unique(intersect(names(test.snei.ex)[is.element(test.snei.ex, s.long.ex)],
                                                names(test.snei.ex)[is.element(test.snei.ex, short.nei.ex)]))
                snei.included.tx <- each.tx.info[is.element(as.matrix(each.tx.info["TXID"]),
                                                            snei.longex), "TXNAME"]
                snei.skipped.tx <- each.tx.info[is.element(as.matrix(each.tx.info["TXID"]),
                                                           snei.shortex), "TXNAME"]

                temp <- NULL
                if(short.nei.ex == long.nei.ex){

                    if(length(snei.included.tx) & length(snei.skipped.tx) > 0){
                        each.ass.result <- cbind(each.ass.result, NeighborEX =
                                                     each.ass.result[,"ShortNeighborEX"])
                        each.ass.result <- rbind(each.ass.result[,c("EnsID",
                                                                    "Nchr",
                                                              "Strand",
                                                              "ShortEX",
                                                              "LongEX",
                                                              "NeighborEX",
                                                              "Types")])
                        temp <- each.ass.result
                    }
                } else {
                    if(length(snei.skipped.tx) * length(snei.included.tx) > 0){
                        skipped.ok = FALSE
                        included.ok = FALSE
                        skipped <- each.tx.info[each.tx.info$TXNAME %in%
                                                    snei.skipped.tx, "TXID"]
                        for(x in skipped){
                            exons <- as.matrix(data.frame(Exon.info[[as.character(x)]],
                                                          stringsAsFactors = F))
                            exons <- cbind(exons, locus = paste(exons[,"start"],
                                                                exons[,"end"],
                                                                sep = "-"))
                            short.num <- as.numeric(rownames(exons[exons[,"locus"] == s.short.ex,]))
                            nei.num <- as.numeric(rownames(exons[exons[,"locus"] == short.nei.ex,]))
                            if(length(short.num) & length(nei.num)){
                                if(abs(short.num - nei.num) == 1){
                                    skipped.ok = TRUE
                                }
                            }
                        }
                        included <- each.tx.info[each.tx.info$TXNAME %in% snei.included.tx, "TXID"]
                        for(x in included){
                            exons <- as.matrix(data.frame(Exon.info[[as.character(x)]],
                                                          stringsAsFactors = F))
                            exons <- cbind(exons, locus = paste(exons[,"start"],
                                                                exons[,"end"],
                                                                sep = "-"))
                            long.num <- as.numeric(rownames(exons[exons[,"locus"] == s.long.ex,]))
                            nei.num <- as.numeric(rownames(exons[exons[,"locus"] == short.nei.ex,]))
                            if(length(long.num) & length(nei.num)){
                                if(abs(long.num - nei.num) == 1){
                                    included.ok = TRUE
                                }
                            }
                        }
                        if(skipped.ok & included.ok){
                            short.ok = TRUE
                        }
                    }

                    test.lnei.tx <- c(intersect(names(each.exon.locus[(is.element(each.exon.locus, s.short.ex))]),
                                                names(each.exon.locus[(is.element(each.exon.locus, long.nei.ex))])),
                                      intersect(names(each.exon.locus[(is.element(each.exon.locus, s.long.ex))]),
                                                names(each.exon.locus[(is.element(each.exon.locus, long.nei.ex))])))
                    test.lnei.ex <- each.exon.locus[is.element(names(each.exon.locus), test.lnei.tx)]
                    lnei.shortex <- unique(intersect(names(test.lnei.ex)[is.element(test.lnei.ex, s.short.ex)],
                                                     names(test.lnei.ex)[is.element(test.lnei.ex, long.nei.ex)]))
                    lnei.longex <- unique(intersect(names(test.lnei.ex)[is.element(test.lnei.ex, s.long.ex)],
                                                    names(test.lnei.ex)[is.element(test.lnei.ex, long.nei.ex)]))
                    lnei.included.tx <- each.tx.info[is.element(as.matrix(each.tx.info["TXID"]),
                                                                lnei.longex), "TXNAME"]
                    lnei.skipped.tx <- each.tx.info[is.element(as.matrix(each.tx.info["TXID"]),
                                                               lnei.shortex), "TXNAME"]

                    if(length(lnei.skipped.tx) * length(lnei.included.tx) > 0){
                        skipped.ok = FALSE
                        included.ok = FALSE
                        skipped <- each.tx.info[each.tx.info$TXNAME %in% lnei.skipped.tx, "TXID"]
                        for(x in skipped){
                            exons <- as.matrix(data.frame(Exon.info[[as.character(x)]],
                                                          stringsAsFactors = F))
                            exons <- cbind(exons, locus = paste(exons[,"start"],
                                                                exons[,"end"],
                                                                sep = "-"))
                            short.num <- as.numeric(rownames(exons[exons[,"locus"] == s.short.ex,]))
                            nei.num <- as.numeric(rownames(exons[exons[,"locus"] == long.nei.ex,]))
                            if(length(nei.num) & length(short.num)){
                                if(abs(short.num - nei.num) == 1){
                                    skipped.ok = TRUE
                                }
                            }
                        }
                        included <- each.tx.info[each.tx.info$TXNAME %in% lnei.included.tx, "TXID"]
                        for(x in included){
                            exons <- as.matrix(data.frame(Exon.info[[as.character(x)]],
                                                          stringsAsFactors = F))
                            exons <- cbind(exons, locus = paste(exons[,"start"],
                                                                exons[,"end"],
                                                                sep = "-"))
                            long.num <- as.numeric(rownames(exons[exons[,"locus"] == s.long.ex,]))
                            nei.num <- as.numeric(rownames(exons[exons[,"locus"] == long.nei.ex,]))
                            if(length(long.num) & length(nei.num)){
                                if(abs(long.num - nei.num) == 1){
                                    included.ok = TRUE
                                }
                            }
                        }
                        if(skipped.ok & included.ok){
                            long.ok = TRUE
                        }
                    }

                    if(short.ok){
                        each.ass.short <- each.ass.result
                        each.ass.short <- cbind(each.ass.short, NeighborEX = each.ass.result[,"ShortNeighborEX"])
                        each.ass.short <- rbind(each.ass.short[,c("EnsID", "Nchr", "Strand", "ShortEX", "LongEX",
                                                            "NeighborEX", "Types")])
                        temp <- rbind(temp, rbind(each.ass.short))
                    }

                    if(long.ok){
                        each.ass.long <- each.ass.result
                        each.ass.long <- cbind(each.ass.result, NeighborEX = each.ass.result[,"LongNeighborEX"])
                        each.ass.long <- rbind(each.ass.long[,c("EnsID", "Nchr", "Strand", "ShortEX", "LongEX",
                                                                "NeighborEX", "Types")])
                        temp <- rbind(temp, rbind(each.ass.long))
                    }
                }
                temp
            })
            each.ass.sane.sub <- do.call(rbind, each.ass.sane.sub)
        }
        each.ass.sane <- rbind(each.ass.sane, each.ass.sane.sub)
        each.ass.sane
    })
    ass.sane <- do.call(rbind, ass.sane)
    ass.sane <- rbind(ass.sane[!duplicated(paste(ass.sane[,"Nchr"],
                                                 ass.sane[,"ShortEX"],
                                                 ass.sane[,"LongEX"],
                                                 ass.sane[,"NeighborEX"],
                                                 sep = ":")),])
    asstypes <- ifelse(ass.sane[,"Strand"] == "+", ass.sane[,"Types"],
                             ifelse(ass.sane[,"Types"] == "A5SS",
                                    "A3SS", "A5SS"))
    a5ss.sane <- rbind(ass.sane[asstypes == "A5SS", c("EnsID", "Nchr", "Strand",
                                                "ShortEX" ,"LongEX",
                                                "NeighborEX")])
    a3ss.sane <- rbind(ass.sane[asstypes == "A3SS", c("EnsID", "Nchr", "Strand",
                                                "ShortEX" ,"LongEX",
                                                "NeighborEX")])

    #Validation for IR
    ri.genes <- unique(ri[,"EnsID"])
    ri.sane <- lapply(1:length(ri.genes), function(i){
        each.ri <- rbind(ri[ri[,"EnsID"] == ri.genes[i],])
        chr.ri <- unique(each.ri[,"Nchr"])
        each.ri.sane <- NULL
        for(chr in chr.ri){
            each.ri.sub <- rbind(each.ri[each.ri[,"Nchr"] == chr,])
            each.tx.info <- sub.txTable[sub.txTable$GENEID == unique(each.ri.sub[,"EnsID"]),]
            each.tx.info <- each.tx.info[each.tx.info$TXCHROM == chr,]
            Exon.info <- findAlternative(unique(each.ri.sub[,"EnsID"]),
                                         sub.txTable, sub.exon.range,
                                         sub.intron.range, chr)$exonRange
            exon.start <- unlist(start(Exon.info))
            exon.end <- unlist(end(Exon.info))
            exon.mat <- cbind(exon.start, exon.end)
            each.exon.locus <- paste(exon.start, exon.end, sep = "-")
            names(each.exon.locus) <- names(exon.end)
            sane <- sapply(1:length(each.ri.sub[,"RetainEX"]),function(each.nums){
                index.num <- each.ri.sub[each.nums, "Index"]
                each.targets <- each.ri.sub[each.nums, "RetainEX"]
                each.ri.result <- rbind(each.ri.sub[each.nums, ])
                s.re.ex <- each.ri.result[,"RetainEX"]
                s.do.ex <- each.ri.result[,"DownEX"]
                s.up.ex <- each.ri.result[,"UpEX"]
                test.do.tx <- names(each.exon.locus[is.element(each.exon.locus,
                                                               s.do.ex)])
                test.up.tx <- names(each.exon.locus[is.element(each.exon.locus,
                                                               s.up.ex)])
                test.re.tx <- names(each.exon.locus[is.element(each.exon.locus,
                                                               s.re.ex)])
                test.tx.ex <- each.exon.locus[is.element(names(each.exon.locus),
                                                         c(intersect(test.do.tx,
                                                                     test.up.tx),

                                                           test.re.tx))]
                do.tx <- unique(names(test.tx.ex)[is.element(test.tx.ex,
                                                             s.do.ex)])
                up.tx <- unique(names(test.tx.ex)[is.element(test.tx.ex,
                                                             s.up.ex)])
                re.tx <- unique(names(test.tx.ex)[is.element(test.tx.ex,
                                                             s.re.ex)])
                included.tx <- unique(intersect(do.tx, up.tx))
                skipped.tx <- re.tx
                included.tx <- each.tx.info[is.element(as.matrix(each.tx.info["TXID"]),
                                                       included.tx), "TXNAME"]
                skipped.tx <- each.tx.info[is.element(as.matrix(each.tx.info["TXID"]),
                                                      skipped.tx), "TXNAME"]
                length(included.tx) != 0 & length(skipped.tx) != 0
            })
            each.ri.sane <- rbind(each.ri.sane, rbind(each.ri.sub[sane,]))
        }
        each.ri.sane
    }
    )
    ri.sane <- do.call(rbind, ri.sane)
    ri.sane <- ri.sane[!duplicated(paste(ri.sane[,"Nchr"], ri.sane[,"RetainEX"],
                                         ri.sane[,"DownEX"],
                                         ri.sane[,"UpEX"], sep = ":")),]
    ri.sane <- rbind(ri.sane[,c("EnsID", "Nchr", "Strand", "RetainEX",
                                "DownEX", "UpEX")])

    return(list(A5SS = a5ss.sane, A3SS = a3ss.sane, SE = se.sane,
                MXE = mxe.sane, RI = ri.sane))

}
