#'Detection of skipped and mutually exclusive exons
#'
#'This function detects skipped (SE) and mutually exclusive exons (MXE) events.
#'
#'@param exon.info A matrix containing exon locus
#'@param intron.info A matrix containing intron locus
#'@param alt.intron.info A matrix containing flanking introns of
#'alternative exons
#'@return A matrix containing SE and MXE events
#'@details This function is borrowed from the \code{IVAS} package.
#'@references Han, A. et al. (2017). Genome wide discovery of genetic variants
#'affecting alternative splicing patterns in human using bioinformatics method.
#'\emph{Genes & Genomics}, 39.
#'@keywords internal
#'@importFrom stats complete.cases
detectSEandMXE <- function(exon.info, intron.info, alt.intron.info){
    tx.exon.ranges <- cbind(tapply(exon.info[,"start"], rownames(exon.info),
                                   min),
                            tapply(exon.info[,"start"], rownames(exon.info),
                                   max))
    colnames(tx.exon.ranges) <- c("start", "end")
    alt.intron.info <- alt.intron.info[1,]
    fi.exon.num <- grep(alt.intron.info["start"] - 1, exon.info[,"end"])
    se.exon.num <- grep(alt.intron.info["end"] + 1, exon.info[,"start"])
    SE.EX <- NULL
    for (i in 1:length(fi.exon.num)){
        fi.exon <- exon.info[fi.exon.num[i],]
        if (length(which(fi.exon["start"] > intron.info[,"start"] - 1 &
                         fi.exon["end"] < intron.info[,"end"] - 1)) != 0){
            SE.EX <- rbind(SE.EX,fi.exon)
        }
    }
    for (i in 1:length(se.exon.num)){
        se.exon <- exon.info[se.exon.num[i],]
        if (length(which(se.exon["start"] > intron.info[,"start"] - 1 &
                         se.exon["end"] < intron.info[,"end"] - 1)) != 0){
            SE.EX <- rbind(SE.EX,se.exon)
        }
    }
    if(length(SE.EX) == 0) return (NULL)
    SE.EX <- rbind(SE.EX[SE.EX[,"start"] != SE.EX[,"end"],])
    if(length(SE.EX) == 0) return (NULL)
    SE.EX <- unique(SE.EX)
    get.test.tx <- function(tx.cal){
        each.tx.cal <- as.double(strsplit(tx.cal,"-")[[1]])
        each.test.tx <- rownames(tx.exon.ranges)[tx.exon.ranges[,"start"] <
                                                     each.tx.cal[1] &
                                                     tx.exon.ranges[,"end"] >
                                                     each.tx.cal[2]]
    }
    test.tx <- unique(unlist(lapply(paste(SE.EX[,"start"], SE.EX[,"end"],
                                          sep = "-"), get.test.tx)))
    exon.info <- rbind(exon.info[is.element(rownames(exon.info), test.tx),])
    if (nrow(exon.info) == 1) rownames(exon.info) <- test.tx
    intron.info <- rbind(intron.info[is.element(rownames(intron.info),
                                                test.tx),])
    if (nrow(intron.info) == 1) rownames(intron.info) <- test.tx
    rownames(SE.EX) <- 1:nrow(SE.EX)
    get.ES.int.num <- function(tx.cal){
        each.tx.se <- as.double(strsplit(tx.cal,"-")[[1]])
        int.num <- which(intron.info[,"start"] + 1 < each.tx.se[1] &
                             intron.info[,"end"]-1 > each.tx.se[2])
        int.num
    }
    ES.int.num <- unique(unlist(lapply(paste(SE.EX[,"start"], SE.EX[,"end"],
                                             sep = "-"), get.ES.int.num)))
    if (length(ES.int.num) == 0) return(NULL)
    ES.int.mat <- cbind(as.double(rownames(intron.info)[ES.int.num]),
                        rbind(intron.info[ES.int.num,]))
    rm.int.num <- NULL
    for (i in 1:nrow(ES.int.mat)){
        num.over.exons <- rownames(exon.info)[exon.info[,"start"] >
                                                  ES.int.mat[i, "start"] - 1 &
                                                  exon.info[, "end"] <
                                                  ES.int.mat[i, "end"] + 1]
        t.num.over.exons <- table(num.over.exons)
        if (length(which(t.num.over.exons > 2)) != 0){
            rm.int.num <- c(rm.int.num, i)
        }
    }
    if(length(rm.int.num) != 0) {ES.int.mat <- rbind(ES.int.mat[-rm.int.num,])}
    if(length(ES.int.mat) == 0) return (NULL)
    colnames(ES.int.mat) <- c("txid", "start", "end")
    min.int <- min(ES.int.mat[,"start"])
    max.int <- max(ES.int.mat[,"end"])
    final.SE.num <- which(exon.info[,"start"] > min.int - 1 &
                              exon.info[,"end"] < max.int + 1)
    final.SE.EX <- cbind(as.double(rownames(exon.info)[final.SE.num]),
                         rbind(exon.info[final.SE.num,]))
    colnames(final.SE.EX) <- c("txid", "start", "end")
    get.ES.int.num2 <- function(fse){
        each.fse <- as.double(strsplit(fse, "-")[[1]])
        names(each.fse) <- c("start", "end")
        which(intron.info[,"start"] + 1 < each.fse["start"] &
                  intron.info[,"end"] - 1 > each.fse["end"])
    }
    ES.int.num <- lapply(paste(final.SE.EX[,"start"], final.SE.EX[,"end"],
                               sep = "-"), get.ES.int.num2)
    ES.int.num <- unique(unlist(ES.int.num))
    ES.int.mat <- cbind(as.double(rownames(intron.info)[ES.int.num]),
                        rbind(intron.info[ES.int.num,]))
    colnames(ES.int.mat) <- c("txid", "start", "end")
    colnames(final.SE.EX) <- c("txid", "start", "end")
    final.SE.EX <- rbind(final.SE.EX[is.element(final.SE.EX[,"txid"],
                                                as.double(test.tx)),])
    t.final.SE.EX <- table(final.SE.EX[,"txid"])
    up2.txid <- as.double(names(t.final.SE.EX)[t.final.SE.EX > 2])
    up1.txid <- as.double(names(t.final.SE.EX)[t.final.SE.EX == 2])
    if (length(ES.int.mat) == 0 | length(final.SE.EX) == 0) return (NULL)
    get.alt.result <- function(total.ES.int){
        ES.int <- as.double(strsplit(total.ES.int, "-")[[1]])
        names(ES.int) <- colnames(ES.int.mat)
        alt.ex <- exon.info[rownames(exon.info)==ES.int["txid"] &
                                (exon.info[,"end"] == ES.int["start"] - 1 |
                                     exon.info[,"start"] == ES.int["end"] + 1),]
        alt.ex.result <- cbind(paste(alt.ex[,"start"], alt.ex[,"end"],
                                     sep = "-"))
    }
    alt.result <- lapply(paste(ES.int.mat[,1], ES.int.mat[,2], ES.int.mat[,3],
                               sep = "-"), get.alt.result)
    pre.alt.result <- NULL
    if (length(alt.result) != 0){
        for (i in 1:length(alt.result)){
            pre.alt.result <- rbind(pre.alt.result, cbind(alt.result[[i]][1],
                                                          alt.result[[i]][2]))
        }
    }
    pre.alt.result <- pre.alt.result[!is.element(pre.alt.result[,1],
                                                 paste(final.SE.EX[,"start"],
                                                       final.SE.EX[,"end"],
                                                       sep = "-")) &
                                         !is.element(pre.alt.result[,2],
                                                     paste(final.SE.EX[,"start"],
                                                           final.SE.EX[,"end"],
                                                           sep = "-")),]
    alt.result <- rbind(unique(pre.alt.result))
    alt.result <- alt.result[complete.cases(alt.result),,drop = FALSE]
    skipped.flank.ex <- alt.result
    up2.es.mat <- NULL
    up2.es.mat.2 <- NULL
    if (length(up1.txid) != 0){
        for (i in 1:length(up1.txid)){
            up1.final.SE.EX <- final.SE.EX[final.SE.EX[,"txid"] == up1.txid[i],]
            min.start <- as.double(min(up1.final.SE.EX[,"start"]))
            max.end <- as.double(max(up1.final.SE.EX[,"end"]))
            le.int <- intron.info[rownames(intron.info) == up1.txid[i] &
                                      intron.info[,"end"] ==
                                      min.start - 1, "start"] - 1
            ri.int <- intron.info[rownames(intron.info) == up1.txid[i] &
                                      intron.info[,"start"] == max.end + 1,
                                  "end"] + 1
            le.ex <- unique(rbind(exon.info[rownames(exon.info) == up1.txid[i] &
                                                exon.info[,"end"] == le.int,]))
            ri.ex <- unique(rbind(exon.info[rownames(exon.info) == up1.txid[i] &
                                                exon.info[,"start"] == ri.int,]))
            le.ex <- paste(le.ex[,"start"], le.ex[,"end"], sep = "-")
            ri.ex <- paste(ri.ex[,"start"], ri.ex[,"end"], sep = "-")
            le.ex.des <- paste(sort(unique(paste(exon.info[exon.info[,"end"] ==
                                                               le.int,"start"],
                                                 exon.info[exon.info[,"end"] ==
                                                               le.int,"end"],
                                                 sep = "-"))), collapse = ",")
            ri.ex.des <- paste(sort(unique(paste(exon.info[exon.info[,"start"] ==
                                                               ri.int,"start"],
                                                 exon.info[exon.info[,"start"] ==
                                                               ri.int,"end"],
                                                 sep = "-"))),collapse = ",")
            if (length(le.ex) != 0 & length(ri.ex) != 0){
                targets.exs <- rbind(c(paste(up1.final.SE.EX[1, -1],
                                             collapse = "-"),
                                       paste(up1.final.SE.EX[2, -1],
                                             collapse = "-")))
                up2.es.mat <- rbind(up2.es.mat, cbind(targets.exs,
                                                      paste(le.ex,
                                                            collapse = ","),
                                                      paste(ri.ex,
                                                            collapse = ","),
                                                      targets.exs, le.ex.des,
                                                      ri.ex.des, "SE"))
                if (length(alt.result) != 0){
                    up2.es.mat.2 <- cbind(paste(up1.final.SE.EX[1, -1],
                                                collapse = "-"),
                                          paste(up1.final.SE.EX[2, -1],
                                                collapse = "-"),
                                          rbind(alt.result), "SE")
                }
            }
        }
    }
    get.inclu.flack.ex <- function(total.in.ex){
        in.ex <- as.double(strsplit(total.in.ex, "-")[[1]])
        names(in.ex) <- c("txid", "start", "end")
        ES.mat.2 <- NULL
        if (length(alt.result) != 0){
            ES.mat.2 <- cbind(paste(in.ex[-1], collapse = "-"), "NA", alt.result,
                              "SE")
        }
        ES.mat.2 <- rbind(ES.mat.2, up2.es.mat.2)
        if (length(ES.mat.2) != 0){
            get.ES.mat.2 <- function(each.num){
                each.ES.mat.2 <- ES.mat.2[each.num,]
                pre.do.des <- unique(rbind(exon.info[exon.info[,"end"] ==
                                                         unlist(strsplit(each.ES.mat.2[3], "-"))[2],]))
                pre.up.des <- unique(rbind(exon.info[exon.info[,"start"] ==
                                                         unlist(strsplit(each.ES.mat.2[4], "-"))[1],]))
                each.ES.mat.2.do.des <- paste(sort(paste(pre.do.des[,"start"],
                                                         pre.do.des[,"end"],
                                                         sep = "-")),
                                              collapse = ",")
                each.ES.mat.2.up.des <- paste(sort(paste(pre.up.des[,"start"],
                                                         pre.up.des[,"end"],
                                                         sep = "-")),
                                              collapse = ",")
                c(each.ES.mat.2[-5], each.ES.mat.2[1:2], each.ES.mat.2.do.des,
                  each.ES.mat.2.up.des, each.ES.mat.2[5])
            }
            ES.mat.2 <- do.call(rbind, lapply(1:nrow(ES.mat.2), get.ES.mat.2))
        }
        le.intron <- rbind(intron.info[intron.info[, "end"] + 1 ==
                                           in.ex["start"],])
        le.exon.num <- is.element(exon.info[,"end"], le.intron[,"start"] - 1)
        le.exon <- cbind(rownames(exon.info)[le.exon.num],
                         rbind(exon.info[le.exon.num,]))
        colnames(le.exon) <- c("txnum", "le.start", "le.end")
        ri.intron <- rbind(intron.info[intron.info[,"start"] - 1 ==
                                           in.ex["end"],])
        ri.exon.num <- is.element(exon.info[,"start"], ri.intron[,"end"] + 1)
        ri.exon <- cbind(rownames(exon.info)[ri.exon.num],
                         rbind(exon.info[ri.exon.num,]))
        colnames(ri.exon) <- c("txnum", "ri.start", "ri.end")
        inExTx <- in.ex["txid"]
        le.exon <- rbind(le.exon[is.element(le.exon[,"txnum"], inExTx),])
        ri.exon <- rbind(ri.exon[is.element(ri.exon[,"txnum"], inExTx),])
        le.ri.exon <- unique(S4Vectors::merge(le.exon, ri.exon, by.x="txnum",
                                   by.y="txnum")[,2:5])
        if (nrow(le.ri.exon) != 0){
            get.ES.result <- function(ex.nums){
                each.le.ri.exon <- rbind(le.ri.exon[ex.nums,])
                le.des <- unique(rbind(exon.info[exon.info[,"end"] ==
                                                     each.le.ri.exon[,"le.end"],]))
                ri.des <- unique(rbind(exon.info[exon.info[,"start"] ==
                                                     each.le.ri.exon[,"ri.start"],]))
                p.le.des <- paste(sort(paste(le.des[,"start"], le.des[,"end"],
                                             sep = "-")), collapse = ",")
                p.ri.des <- paste(sort(paste(ri.des[,"start"], ri.des[,"end"],
                                             sep = "-")), collapse = ",")
                cbind(paste(in.ex[-1], collapse = "-"), paste(each.le.ri.exon[,"le.start"],
                                                              each.le.ri.exon[,"le.end"],
                                                              sep = "-"),
                      paste(each.le.ri.exon[,"ri.start"], each.le.ri.exon[,"ri.end"],
                            sep = "-"),
                      paste(in.ex[-1], collapse = "-"), "NA",
                      p.le.des, p.ri.des, "SE")
            }
            ES.result <- do.call(rbind, lapply(1:nrow(le.ri.exon),
                                               get.ES.result))
            ES.result <- cbind(ES.result[,1], "NA", rbind(ES.result[,-1]))
            ES.result <- rbind(ES.result, up2.es.mat)
            final.ES.result <- NULL
            if (length(ES.mat.2) != 0){
                for (i in 1:nrow(ES.result)){
                    for (j in 1:nrow(ES.mat.2)){
                        final.ES.result <- rbind(final.ES.result,
                                                 c(ES.result[i, -c(3:9)],
                                                   ES.mat.2[j, 3],
                                                   ES.result[i, 4],
                                                   ES.result[i, -c(3:9)],
                                                   ES.mat.2[j, 7],
                                                   ES.result[i, 8],
                                                   ES.result[i, 9]))
                        final.ES.result <- rbind(final.ES.result,
                                                 c(ES.result[i, -c(4:9)],
                                                   ES.mat.2[j, 4],
                                                   ES.result[i, 5:7],
                                                   ES.mat.2[j, 8],
                                                   ES.result[i, 9]))
                    }
                }
            }
            ES.result <- unique(rbind(ES.result, final.ES.result, ES.mat.2))
            ES.result
        }

    }
    inclu.flack.ex <- lapply(paste(final.SE.EX[,"txid"], final.SE.EX[,"start"],
                                   final.SE.EX[,"end"], sep = "-"),
                             get.inclu.flack.ex)
    ES.mat <- do.call(rbind,inclu.flack.ex)
    if (length(ES.mat) == 0) return(NULL)
    ES.mat <- unique(ES.mat)
    sec.mxe.mat <- rbind(ES.mat[ES.mat[,2] != "NA",])
    u.1st.ex <- unique(ES.mat[ES.mat[,2] == "NA",1])
    u.do.ex <- unique(ES.mat[ES.mat[,2] == "NA",3])
    u.up.ex <- unique(ES.mat[ES.mat[,2] == "NA",4])
    final.pre.ES.mat <- NULL
    for (i in 1:length(u.1st.ex)){
        for (j in 1:length(u.do.ex)){
            if (u.1st.ex[i] == u.do.ex[j])    next
            do.des <- unique(ES.mat[ES.mat[,3] == u.do.ex[j], 7])
            pre.ES.mat <- do.call(rbind, lapply(u.up.ex, function(each.up.ex){
                if (u.1st.ex[i] != each.up.ex){
                    up.des <- unique(ES.mat[ES.mat[,4] == each.up.ex,8])
                    cbind(u.1st.ex[i], "NA", u.do.ex[j], each.up.ex,
                          u.1st.ex[i], "NA", do.des,up.des, "SE")
                }
            }))
            final.pre.ES.mat <- rbind(final.pre.ES.mat, pre.ES.mat)
        }
    }
    ES.mat <- unique(rbind(final.pre.ES.mat, sec.mxe.mat))
    colnames(ES.mat) <- c("TarEX", "2ndEX", "DownEX", "UpEX", "Tar_des",
                          "2nd_des", "Do_des", "Up_des", "Types")
    final.MXE.cal <- NULL
    pre.MXE.mat <- NULL
    MXE.ES.mat <- NULL
    final.down.ex <- unique(do.call(rbind, strsplit(ES.mat[,"DownEX"], "-")))
    final.up.ex <- unique(do.call(rbind, strsplit(ES.mat[,"UpEX"], "-")))
    colnames(final.down.ex) <- c("start", "end")
    colnames(final.up.ex) <- c("start", "end")
    u.final.SE.EX <- cbind("1", rbind(unique(final.SE.EX[, c("start", "end")])))
    colnames(u.final.SE.EX) <- c("txid", "start", "end")
    get.final.MXE.cal <- function(total.es){
        MXE.mat <- NULL
        each.es <- as.double(strsplit(total.es, "-")[[1]])
        names(each.es) <- c("txid", "start", "end")
        pre.SE.intron <- rbind(ES.mat[is.element(ES.mat[,"TarEX"],
                                                 paste(each.es["start"],
                                                       each.es["end"],
                                                       sep = "-")),])
        if (length(pre.SE.intron) != 0){
            SE.intron <- cbind(do.call(rbind, strsplit(pre.SE.intron[,"DownEX"], "-"))[,2],
                               do.call(rbind, strsplit(pre.SE.intron[,"TarEX"], "-"))[,1])
            SE.intron <- rbind(SE.intron, cbind(do.call(rbind, strsplit(pre.SE.intron[,"TarEX"], "-"))[,2],
                                                do.call(rbind, strsplit(pre.SE.intron[,"UpEX"], "-"))[,1]))
            u.SE.intron <- rbind(unique(SE.intron))
            get.MXE.ex.test <- function(mxe.intron){
                mxe.intron <- as.integer(unlist(strsplit(mxe.intron,"-")))
                do.mxe <- final.down.ex[as.integer(final.down.ex[,"start"]) >
                                            mxe.intron[1] - 1 &
                                            as.integer(final.down.ex[,"end"]) <
                                            mxe.intron[2] + 1,]
                up.mxe <- final.up.ex[as.integer(final.up.ex[,"start"]) >
                                          mxe.intron[1] - 1 &
                                          as.integer(final.up.ex[,"end"]) <
                                          mxe.intron[2] + 1,]
                rbind(do.mxe,up.mxe)
            }
            MXE.ex.test <- lapply(paste(u.SE.intron[,1], u.SE.intron[,2],
                                        sep = "-"), get.MXE.ex.test)
            MXE.ex.test <- do.call(rbind, MXE.ex.test)
            if (length(MXE.ex.test) != 0){
                MXE.ex.test <- unique(rbind(MXE.ex.test))
                get.MXE.mat <- function(total.other.ex){
                    other.ex <- unlist(strsplit(total.other.ex,"-"))
                    if (as.double(each.es["end"]) < as.double(other.ex[1])){
                        MXE.cal <- cbind(paste(each.es[-1], collapse = "-"),
                                         paste(other.ex, collapse = "-"))
                    }
                    else if (as.double(each.es["start"]) > as.double(other.ex[2])){
                        MXE.cal <- cbind(paste(other.ex, collapse = "-"),
                                         paste(each.es[-1], collapse = "-"))
                    }
                }
                MXE.mat <- do.call(rbind, lapply(paste(MXE.ex.test[,1],
                                                       MXE.ex.test[,2],
                                                       sep = "-"),
                                                 get.MXE.mat))
            }
        }
        MXE.mat
    }
    final.MXE.cal <- do.call(rbind, lapply(paste(u.final.SE.EX[,1],
                                                 u.final.SE.EX[,2],
                                                 u.final.SE.EX[,3],
                                                 sep = "-"),
                                           get.final.MXE.cal
    ))
    if (length(final.MXE.cal) != 0){
        MXE.cal <- final.MXE.cal
        MXE.cal <- unique(MXE.cal)
        if(length(MXE.cal) != 0){
            get.MXE.ES.mat <- function(total.mem){
                mem <- strsplit(total.mem,"-")[[1]]
                ES.do <- do.call(rbind,strsplit(ES.mat[,"DownEX"], "-"))
                ES.up <- do.call(rbind,strsplit(ES.mat[,"UpEX"], "-"))
                MXE.ES <- unique(rbind(ES.mat[as.double(ES.do[,1]) <
                                                  min(as.double(mem)) &
                                                  as.double(ES.up[,2]) >
                                                  max(as.double(mem)),
                                              c("DownEX", "UpEX", "Do_des",
                                                "Up_des", "Types")]))
                if (length(MXE.ES) != 0){
                    tar.exs <- matrix(rep(c(paste(mem[1:2], collapse = "-"),
                                            paste(mem[3:4], collapse = "-")),
                                          nrow(MXE.ES)), nrow=nrow(MXE.ES),
                                      byrow = TRUE)
                    cbind(tar.exs, rbind(MXE.ES[,c("DownEX", "UpEX")]),
                          tar.exs, rbind(MXE.ES[,c("Do_des", "Up_des",
                                                   "Types")]))
                }
            }
            MXE.ES.mat <- do.call(rbind, lapply(paste(MXE.cal[,1], MXE.cal[,2],
                                                      sep = "-"),
                                                get.MXE.ES.mat))
            if (length(MXE.ES.mat) != 0){
                colnames(MXE.ES.mat) <- c("1stEX", "2ndEX", "DownEX",
                                          "UpEX", "1st_des", "2nd_des",
                                          "Do_des", "Up_des", "Types")
                MXE.ES.mat[,"Types"] <- "MXE"
            }
        }
        ES.mat <- rbind(ES.mat, unique(MXE.ES.mat))
    }
    ES.mat <- unique(ES.mat)
    down.ex.e <- as.double(do.call(rbind,strsplit(ES.mat[, "DownEX"], "-"))[,2])
    fi.ex.s.e <- do.call(rbind,strsplit(ES.mat[, "TarEX"],"-"))
    up.ex.e <- as.double(do.call(rbind,strsplit(ES.mat[, "UpEX"],"-"))[,1])
    ES.mat <- rbind(ES.mat[down.ex.e < as.double(fi.ex.s.e[,1]) &
                               as.double(fi.ex.s.e[,2]) < up.ex.e,])
    ES.mat <- cbind(ES.mat,"possible")
    colnames(ES.mat) <- c("1stEX", "2ndEX", "DownEX", "UpEX", "1st_des",
                          "2nd_des", "Do_des", "Up_des", "Types", "status")
    po.test.mat <- paste(do.call(rbind, strsplit(ES.mat[,"DownEX"], "-"))[,2],
                         do.call(rbind, strsplit(ES.mat[,"1stEX"], "-"))[,1],
                         sep = "-")
    po.test.mat <- cbind(po.test.mat,
                         paste(do.call(rbind, strsplit(ES.mat[,"1stEX"],
                                                       "-"))[,2],
                               do.call(rbind, strsplit(ES.mat[,"UpEX"],
                                                      "-"))[,1], sep = "-"))
    po.test.mat <- cbind(po.test.mat,
                         paste(do.call(rbind, strsplit(ES.mat[,"DownEX"],
                                                       "-"))[,2],
                               do.call(rbind,strsplit(ES.mat[,"UpEX"],
                                                      "-"))[,1], sep = "-"))
    po.test.mat <- cbind(po.test.mat, rbind(ES.mat[,3:4]))
    intron.info.mat <- paste(intron.info[,"start"] - 1, intron.info[,"end"] + 1,
                             sep = "-")
    names(intron.info.mat) <- rownames(intron.info)
    exon.info.mat <- paste(exon.info[,"start"], exon.info[,"end"], sep = "-")
    names(exon.info.mat) <- rownames(exon.info)
    po.over.result <- lapply(1:nrow(po.test.mat), function(total.ptm.num){
        ptm <- po.test.mat[total.ptm.num,]
        ptm.test.int <- which(is.element(names(intron.info.mat)[which(intron.info.mat == ptm[1])],
                                         names(intron.info.mat)[which(intron.info.mat == ptm[2])]) == "TRUE")
        ptm.test.ex <- intersect(names(exon.info.mat)[which(exon.info.mat == ptm[4])],
                                 names(exon.info.mat)[which(exon.info.mat == ptm[5])])
        ptm.test.ski.int <- names(intron.info.mat)[which(intron.info.mat == ptm[3])]
        if(length(intersect(ptm.test.ex,ptm.test.ski.int)) != 0 & length(ptm.test.int) != 0){
            total.ptm.num
        }
    })
    po.over.result <- is.element(1:nrow(ES.mat), unlist(po.over.result))
    po.over.result <- po.over.result & (ES.mat[,"2ndEX"] == "NA")
    ES.mat[po.over.result, "status"] <- "exist"
    se.exons.num <- ES.mat[, "2ndEX"] != "NA" & ES.mat[,"Types"] == "SE"
    se.MXE.num <- ES.mat[,"2ndEX"] != "NA" & ES.mat[,"Types"] == "MXE"
    if(length(which(se.exons.num=="TRUE")) != 0){
        es.se.test <- NULL
        es.se.ex <- rbind(ES.mat[se.exons.num,])
        for (i in 1:nrow(es.se.ex)){
            es.se.1 <- which(is.element(names(exon.info.mat)[which(exon.info.mat == es.se.ex[i, 1])],
                                        names(exon.info.mat)[which(exon.info.mat == es.se.ex[i,2])]) == "TRUE")
            es.se.2 <- which(intron.info.mat == paste(unlist(strsplit(es.se.ex[i,3], "-"))[2],
                                                      unlist(strsplit(es.se.ex[i,4], "-"))[1],
                                                      sep = "-"))
            if (length(es.se.1) != 0 & length(es.se.2) != 0){
                es.se.test <- c(es.se.test, 1==1)
            }
        }
        ES.mat[which(se.exons.num == "TRUE")[es.se.test], "status"] <- "exist"
    }
    if(length(which(se.MXE.num == "TRUE")) != 0){
        mxe.se <- rbind(ES.mat[se.MXE.num,])
        mxe.se.1 <- paste(do.call(rbind, strsplit(mxe.se[,"1stEX"], "-"))[,2],
                          do.call(rbind, strsplit(mxe.se[,"UpEX"], "-"))[,1],
                          sep = "-")
        mxe.se.2 <- paste(do.call(rbind, strsplit(mxe.se[,"DownEX"], "-"))[,2],
                          do.call(rbind, strsplit(mxe.se[,"2ndEX"], "-"))[,1],
                          sep = "-")
        mxe.se.3 <- paste(do.call(rbind, strsplit(mxe.se[,"DownEX"], "-"))[,2],
                          do.call(rbind, strsplit(mxe.se[,"1stEX"], "-"))[,1],
                          sep = "-")
        mxe.se.4 <- paste(do.call(rbind,strsplit(mxe.se[, "2ndEX"], "-"))[,2],
                          do.call(rbind,strsplit(mxe.se[, "UpEX"], "-"))[,1],
                          sep = "-")
        mxe.se.test <- is.element(mxe.se.1, intron.info.mat) &
            is.element(mxe.se.2, intron.info.mat) &
            is.element(mxe.se.3, intron.info.mat) &
            is.element(mxe.se.4, intron.info.mat)
        ES.mat[which(se.MXE.num == "TRUE")[mxe.se.test], "status"] <- "exist"
    }
    rownames(ES.mat) <- 1:nrow(ES.mat)
    colnames(ES.mat) <- c("1stEX", "2ndEX", "DownEX", "UpEX", "1st_des",
                          "2nd_des", "Do_des", "Up_des", "Types", "status")
    ES.mat <- rbind(ES.mat[,c("1stEX", "2ndEX", "DownEX", "UpEX", "Types",
                              "status")])
    return (unique(ES.mat))
}
