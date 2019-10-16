#'AS event quantification
#'
#'This function calculates PSI of AS events from input bam files.
#'
#'@param AS.list a list of AS events
#'@param sample.info a data frame containing names, bam file paths, and
#'groups (optional) of samples
#'@param read.type a type of RNA-seq reads ("single" or "paired")
#'@param read.length read length
#'@param insert.size insert size
#'@param min.reads a minimum number of reads mapped to a given exon
#'@param num.cores the number of cores for parallel processing
#'@return A matrix containing PSI of AS events
#'@details This function is modified from the \code{RatioFromReads} function
#'in the \code{IMAS} package. \code{sample.info} must include two columns
#'containing names and bam file paths of samples, respectively. The first,
#'second, and third columns must correspond to names, bam file paths, and groups
#'of samples.
#'@references Han, S. et al. (2017). IMAS: Integrative analysis of
#'Multi-omics data for Alternative Splicing. R package version 1.8.0.
#'@keywords internal
#'@import SummarizedExperiment
#'@importFrom BiocParallel bplapply SnowParam
#'@noRd
quantify <- function(AS.list, sample.info, read.type = "paired",
                     read.length, insert.size, min.reads = 3,
                     num.cores = 1){
    splitSplice <- function(EX1, EX2) {
        s.EX1 <- strsplit(unlist(strsplit(EX1, ",")), "-")
        s.EX2 <- strsplit(unlist(strsplit(EX2, ",")), "-")
        s.EX1 <- matrix(as.integer(do.call(rbind, s.EX1)), ncol = 2)
        s.EX2 <- matrix(as.integer(do.call(rbind, s.EX2)), ncol = 2)
        colnames(s.EX1) <- c("start", "end")
        colnames(s.EX2) <- c("start", "end")
        final.re <- unlist(lapply(s.EX1[, "end"], function(each.fi.s.EX) {
            pre.re <- lapply(s.EX2[, "start"], function(each.se.s.EX) {
                paste(sort(c(each.fi.s.EX, each.se.s.EX)), collapse = "-")
            })
            do.call(rbind, pre.re)
        }))
        return(unique(final.re))
    }
    coorEX <- function(spl.re, g.Info, read.length, inse, min.reads = 5,
                       AStype) {
        Normalized.values <- function(exons.l, inse) {
            normal.max.in <- 0
            normal.max.skip <- 0
            normal.in <- 0
            normal.skip <- 0
            jun.in.normal <- 0
            jun.skip.normal <- 0
            if (AStype == "SE" | AStype == "MXE") {
                max.p.reads <- exons.l[seq_len(3)] + (read.l - 1)
                fi.in.value <- abs(min((1 + inse + read.l) -
                                           (exons.l[1] + 1), 0))
                fi.nex.value <- abs(min(exons.l[2] - (read.l + inse), 0))
                se.in.value <- abs(min((1 + inse + read.l) -
                                           (exons.l[2] + 1), 0))
                se.nex.value <- abs(min(exons.l[3] - (read.l + inse), 0))
                fi.in.minus <- fi.in.value + fi.nex.value
                se.in.minus <- se.in.value + se.nex.value
                se.fi <- abs(min(abs(read.l + inse - exons.l[2]) -
                                     (read.l - 1), 0))
                sk.in.value <- abs(min((1 + inse + read.l) - (exons.l[1] + 1),
                                       0))
                sk.nex.value <- abs(min(exons.l[3] - (read.l + inse), 0))
                sk.minus.value <- sk.in.value + sk.nex.value
                normal.fi.in <- max((max.p.reads[1] - fi.in.minus), 0)
                normal.se.in <- max((max.p.reads[2] - se.in.minus), 0)
                normal.max.in <- max(normal.fi.in + normal.se.in - se.fi, 0)
                normal.max.skip <- max(max.p.reads[1] - sk.minus.value, 0)
                exon.r.a <- abs(min(1 + inse - exons.l[1], 0))
                exon.r.b <- abs(min(exons.l[2] - (read.l + inse), 0))
                exon.r.c <- abs(min((1 + inse) - (exons.l[2] + exons.l[3]), 0))
                exon.r.d <- abs(min(read.l + inse - exons.l[2], 0))
                exon.r.e <- abs(min((inse + 1) - exons.l[3], 0))
                jun.pair.a <- abs(min(read.l - 1, exon.r.a))
                jun.pair.b <- abs(min(read.l - 1, exons.l[1],
                                      exon.r.b))
                jun.pair.c <- abs(min(read.l - 1, exon.r.c))
                jun.pair.d <- abs(min(read.l - 1, exon.r.d))
                jun.pair.e <- abs(min(read.l - 1, exon.r.e))
                jun.sk.pa <- sum(jun.pair.a, jun.pair.e)
                jun.in.pa <- sum(jun.sk.pa, jun.pair.b, jun.pair.c,
                                 jun.pair.d)
                normal.in <- max(normal.max.in - jun.in.pa, 0)
                normal.skip <- max(normal.max.skip - jun.sk.pa, 0)
                jun.in.normal <- max(2 * (2 * (read.l - 2 * an.size + 1)) -
                                         (2 * abs(min(exons.l[2] - (read.l -
                                                                        2 + 1),
                                                      0)) +
                                              abs(min(abs(read.l + inse -
                                                              exons.l[2]) -
                                                          (read.l - 1), 0))), 0)
                jun.skip.normal <- max(2 * (read.l - 2 * an.size + 1), 0)
            }
            else if (AStype == "ASS") {
                normal.in <- 2 * exons.l[2]
                jun.in.normal <- max(2 * (read.l - 2 * an.size + 1), 0)
                jun.skip.normal <- max(2 * (read.l - 2 * an.size + 1), 0)
            }
            else if (AStype == "RI") {
                jun.in.normal <- max(2 * (read.l - 2 * an.size + 1), 0)
                normal.skip <- 2 * (exons.l[2] + (read.l - 1))
            }
            return.mat <- c(normal.max.in, normal.max.skip, normal.in,
                            normal.skip, jun.in.normal, jun.skip.normal)
            cn <- c("pair.in", "pair.sk", "pairwojun.in", "pairwojun.sk",
                    "jun.in", "jun.sk")
            names(return.mat) <- cn
            return(return.mat)
        }
        g.1.info <- g.Info[[1]]
        g.2.info <- g.Info[[2]]
        paired.r <- spl.re$pairedInfo
        exon.r <- spl.re$exonInfo
        junction.r <- spl.re$junctionInfo
        total.exon.l <- 0
        if (sum(c(paired.r, exon.r, junction.r)) == 0) {
            return("NA")
        }
        total.exons <- do.call(rbind, strsplit(names(exon.r),
                                               "-"))
        total.exon.l <- apply(total.exons, 1, function(x) diff(as.double(x) +
                                                                   1))
        g.1.p.r <- paired.r[is.element(names(paired.r), g.1.info$paired)]
        g.2.p.r <- paired.r[is.element(names(paired.r), g.2.info$paired)]
        g.1.e.r <- exon.r[is.element(names(exon.r), g.1.info$exon[2])]
        g.2.e.r <- 0
        j.nm <- names(junction.r)
        g.1.j.r <- junction.r[is.element(j.nm, g.1.info$junction)]
        g.2.j.r <- junction.r[is.element(j.nm, g.2.info$junction)]
        if (!length(g.1.p.r) & !length(g.2.p.r) & !length(g.1.e.r)) {
            if (any(length(g.1.j.r)) & any(length(g.2.j.r))) {
                g.1.p.r <- 0
                g.2.p.r <- 0
                g.1.e.r <- 0
            }
            else return("NA")
        }
        if (!any(length(g.1.e.r))) {
            g.1.e.r <- 0
        }
        if (!any(length(g.1.j.r)))
            g.1.j.r <- 0
        if (!any(length(g.2.j.r)))
            g.2.j.r <- 0
        read.l <- as.integer(read.length)
        an.size = 1
        normal.values <- Normalized.values(total.exon.l, inse)
        if (AStype == "SE" | AStype == "MXE") {
            if (AStype == "MXE") {
                fi.Nvalues <- Normalized.values(total.exon.l[c(1, 2, 4)], inse)
                se.N.values <- Normalized.values(total.exon.l[c(1, 3, 4)], inse)
                fi.Nvalues[c("pair.sk", "pairwojun.sk",
                             "jun.sk")] <- se.N.values[c("pair.in",
                                                         "pairwojun.in",
                                                         "jun.in")]
                normal.values <- fi.Nvalues
            }
            total.reads <- sum(g.1.p.r, g.1.j.r, g.2.p.r, g.2.j.r)
            g.1.nor.num <- sum(normal.values[c("pairwojun.in",
                                               "jun.in")])
            g.2.nor.num <- sum(normal.values[c("pairwojun.sk",
                                               "jun.sk")])
            group.1.read.count <- sum(g.1.p.r, g.1.j.r)
            group.2.read.count <- sum(g.2.p.r, g.2.j.r)
            read.num.test.1 <- sum(g.1.p.r) == 0 & sum(g.1.j.r) == 0
            read.num.test.2 <- sum(g.2.p.r) == 0 & sum(g.2.j.r) == 0
            if ((read.num.test.1 & read.num.test.2) | total.reads < min.reads){
                return("NA")
            }
            if (sum(g.1.j.r) == 0 | sum(g.2.j.r) == 0) {
                g.1.nor.num <- normal.values["pair.in"]
                g.2.nor.num <- normal.values["pair.sk"]
                group.1.read.count <- sum(g.1.p.r)
                group.2.read.count <- sum(g.2.p.r)
            }
            else if (sum(g.1.p.r) == 0 | sum(g.2.p.r) == 0) {
                g.1.nor.num <- normal.values["jun.in"]
                g.2.nor.num <- normal.values["jun.sk"]
                group.1.read.count <- sum(g.1.j.r)
                group.2.read.count <- sum(g.2.j.r)
            }
        }
        else if (AStype == "RI") {
            g.1.nor.num <- normal.values["jun.in"]
            g.2.nor.num <- normal.values["pairwojun.sk"]
            group.1.read.count <- sum(g.1.j.r)
            group.2.read.count <- sum(g.1.e.r)
        }
        else if (AStype == "ASS") {
            g.1.nor.num <- normal.values["jun.in"]
            g.2.nor.num <- normal.values["jun.sk"]
            group.1.read.count <- sum(g.1.j.r)
            group.2.read.count <- sum(g.2.j.r)
            if (g.1.e.r > 5) {
                g.1.nor.num <- sum(normal.values[c("pairwojun.in",
                                                   "jun.in")])
                group.1.read.count <- sum(g.1.e.r, g.1.j.r)
            }
        }
        ratio.g.1 <- group.1.read.count/g.1.nor.num
        ratio.g.2 <- group.2.read.count/g.2.nor.num
        group.total <- sum(ratio.g.1, ratio.g.2)
        group.1.2.ratio <- ifelse(AStype == "RI",
                                  ratio.g.2/group.total,
                                  1 - (ratio.g.2/group.total))
        if (group.1.read.count == 0 & group.2.read.count == 0)
            return("NA")
        if (group.1.read.count < min.reads & group.2.read.count < min.reads){
            group.1.2.ratio <- "NA"
        }
        return(group.1.2.ratio)
    }
    Each.Cal.ratio <- function(bamfiles = NULL, splicingInfo = NULL, splitEnv,
                               cEnv, ins, min.reads, read.length, read.type,
                               parm){
        ExReads <- function(t.ex, t.sp, g.e1, g.e2, g1, g2,
                            s1, s2, ch, er, met, alt, ins, min.reads) {
            coor.re <- NULL
            pre.bam.re <- lapply(seq_along(bamfiles), function(ebam) {
                T.r <- SplicingReads(bamfiles[ebam], t.ex, t.sp,
                                     er, ch, met, ins)
                if(length(T.r$exonInfo) == 0){
                    return("NA")
                }
                if(length(T.r$pairedInfo) == 0){
                    T.r$pairedInfo <- NULL
                }
                group.1.list <- list(g1, g.e1, s1)
                group.2.list <- list(g2, g.e2, s2)
                names(group.1.list) <- c("paired", "exon", "junction")
                names(group.2.list) <- c("paired", "exon", "junction")
                t.g.li <- list(group.1.list, group.2.list)
                names(t.g.li) <- c("Inclu", "Skip")
                coor.re <- cEnv$coorEX(T.r, t.g.li, read.length, ins, min.reads,
                                       alt)
                return(coor.re)
            })
            pre.bam.re <- do.call(cbind, pre.bam.re)
            pre.bam.re[is.na(pre.bam.re)] <- "NA"
            return(pre.bam.re)
        }
        SE.te <- function(SE.num) {
            SE.re <- rbind(SE.result[SE.num,])
            e.dw.st <- do.call(rbind, strsplit(SE.re[,"DownEX"], "-"))
            e.up.en <- do.call(rbind, strsplit(SE.re[,"UpEX"], "-"))
            each.ran <- rbind(unique(cbind(e.dw.st[,1], e.up.en[,2])))
            fi.pr <- paste(SE.re[, c("DownEX", "1stEX")], collapse = "~")
            se.pr <- paste(SE.re[, c("1stEX", "UpEX")], collapse = "~")
            g.1.p <- rbind(fi.pr, se.pr)
            g.2.p <- paste(SE.re[, c("DownEX", "UpEX")], collapse = "~")
            fi.s1 <- splitEnv$splitSplice(SE.re[,"DownEX"], SE.re[,"1stEX"])
            fi.s2 <- splitEnv$splitSplice(SE.re[,"1stEX"], SE.re[,"UpEX"])
            g.1.s <- c(fi.s1, fi.s2)
            g.2.s <- splitEnv$splitSplice(SE.re[,"DownEX"], SE.re[,"UpEX"])
            t.ex <- SE.re[,c("DownEX", "1stEX", "UpEX")]
            t.sp <- c(g.1.s, g.2.s)
            e.chr <- SE.re[,"Nchr"]
            pr.re <- ExReads(t.ex, t.sp, t.ex, t.ex, g.1.p, g.2.p, g.1.s,
                             g.2.s, e.chr, each.ran, read.type, "SE", ins,
                             min.reads)
            pr.re
        }
        MXE.te <- function(MXE.num) {
            MXE.re <- rbind(MXE.result[MXE.num,])
            do.st <- do.call(rbind, strsplit(MXE.re[,"DownEX"], "-"))
            up.en <- do.call(rbind, strsplit(MXE.re[,"UpEX"], "-"))
            do.st <- do.st[,1]
            up.en <- up.en[,2]
            each.ran <- rbind(unique(cbind(do.st, up.en)))
            fi.p.1 <- paste(MXE.re[, c("DownEX", "1stEX")], collapse = "~")
            fi.p.2 <- paste(MXE.re[, c("DownEX", "2ndEX")], collapse = "~")
            se.p.1 <- paste(MXE.re[, c("1stEX", "UpEX")], collapse = "~")
            se.p.2 <- paste(MXE.re[, c("2ndEX", "UpEX")], collapse = "~")
            g.1.p <- rbind(fi.p.1, se.p.1)
            g.2.p <- rbind(fi.p.2, se.p.2)
            fi.s1 <- splitEnv$splitSplice(MXE.re[,"DownEX"], MXE.re[,"1stEX"])
            fi.s2 <- splitEnv$splitSplice(MXE.re[,"1stEX"], MXE.re[,"UpEX"])
            se.s1 <- splitEnv$splitSplice(MXE.re[,"DownEX"], MXE.re[,"2ndEX"])
            se.s2 <- splitEnv$splitSplice(MXE.re[,"2ndEX"], MXE.re[,"UpEX"])
            g.1.s <- c(fi.s1, fi.s2)
            g.2.s <- c(se.s1, se.s2)
            t.sp <- c(g.1.s, g.2.s)
            g.1.e <- c(MXE.re[,"DownEX"], MXE.re[,"1stEX"], MXE.re[,"UpEX"])
            g.2.e <- c(MXE.re[,"DownEX"], MXE.re[,"2ndEX"], MXE.re[,"UpEX"])
            t.ex <- MXE.re[, c("DownEX", "1stEX", "2ndEX", "UpEX")]
            e.chr <- MXE.re[,"Nchr"]
            pr.re <- ExReads(t.ex, t.sp, g.1.e, g.2.e, g.1.p, g.2.p, g.1.s,
                             g.2.s, e.chr, each.ran, read.type,
                             "MXE", ins, min.reads)
            pr.re
        }
        RI.te <- function(RI.num) {
            ea.re <- rbind(RI.result[RI.num, ])
            s.down <- strsplit(ea.re[,"DownEX"], "-")
            s.up <- strsplit(ea.re[,"UpEX"], "-")
            ex.sp <- paste(unlist(s.down)[2], unlist(s.up)[1], sep = "-")
            do.st <- do.call(rbind, s.down)[, 1]
            up.en <- do.call(rbind, s.up)[, 2]
            each.ran <- rbind(unique(cbind(do.st, up.en)))
            g.1.p <- paste(ea.re[, c("DownEX", "UpEX")], collapse = "~")
            g.1.p <- rbind(g.1.p)
            fi.pr.2 <- paste(c(ea.re[,"DownEX"], ex.sp), collapse = "~")
            se.pr.2 <- paste(c(ex.sp, ea.re[,"UpEX"]), collapse = "~")
            g.2.p <- rbind(fi.pair.2, se.pair.2)
            g.1.s <- c(splitEnv$splitSplice(ea.re[,"DownEX"], ea.re[,"UpEX"]))
            g.2.s <- "NA"
            t.ex <- c(ea.re[,"DownEX"], ex.sp, ea.re[,"UpEX"])
            e.chr <- ea.re[,"Nchr"]
            pr.re <- ExReads(t.ex, g.1.s, t.ex, t.ex, g.1.p,
                             g.2.p, g.1.s, g.2.s, e.chr, each.ran, "exon",
                             "RI", ins, min.reads)
            pr.re
        }
        ASS.te <- function(ASS.num) {
            ea.re <- rbind(ASS.result[ASS.num, ])
            ea.ty <- ifelse(ea.re[,"Types"] == "A5SS", 1, -1)
            ea.strand <- ifelse(ea.re[,"Strand"] == "+", 1, -1)
            ASS.nei <- unlist(strsplit(ea.re[,"NeighborEX"], "-"))
            sh.ex <- unlist(strsplit(ea.re[,"ShortEX"], "-"))
            lo.ex <- unlist(strsplit(ea.re[,"LongEX"], "-"))
            sh.ex.1 <- sh.ex[1]
            sh.ex.2 <- sh.ex[2]
            lo.ex.1 <- lo.ex[1]
            lo.ex.2 <- lo.ex[2]
            if (ea.ty * ea.strand == 1) {
                min.r <- max(as.integer(sh.ex)) - 20
                max.r <- min(as.integer(ASS.nei)) + 20
                ASS.sp <- paste(sh.ex.2, lo.ex.2, sep = "-")
                names(ASS.sp) <- "Alt.ex"
                g.2.p <- c(ea.re[,"ShortEX"], ea.re[,"NeighborEX"])
                g.2.p <- rbind(paste(g.2.p, collapse = "~"))
                g.1.s <- c(splitEnv$splitSplice(ea.re[,"LongEX"],
                                                ea.re[,"NeighborEX"]))
                g.2.s <- c(splitEnv$splitSplice(ea.re[,"ShortEX"],
                                                ea.re[,"NeighborEX"]))
                g.1.e <- c(ASS.sp, ea.re[,"NeighborEX"])
                g.2.e <- c(ASS.sp, ea.re[,"NeighborEX"])
                t.ex <- c(ASS.sp, ea.re[,"NeighborEX"])
            }
            else if (ea.ty * ea.strand == -1) {
                min.r <- max(as.integer(ASS.nei)) - 20
                max.r <- min(as.integer(sh.ex)) + 20
                ASS.sp <- paste(lo.ex.1, sh.ex.1, sep = "-")
                names(ASS.sp) <- "Alt.ex"
                g.2.p <- c(ea.re[,"NeighborEX"], ea.re[,"ShortEX"])
                g.2.p <- rbind(paste(g.2.p, collapse = "~"))
                g.1.s <- c(splitEnv$splitSplice(ea.re[,"NeighborEX"],
                                                ea.re[,"LongEX"]))
                g.2.s <- c(splitEnv$splitSplice(ea.re[,"NeighborEX"],
                                                ea.re[,"ShortEX"]))
                g.1.e <- c(ea.re[,"NeighborEX"], ASS.sp)
                g.2.e <- c(ea.re[,"NeighborEX"], ASS.sp)
                t.ex <- c(ea.re[,"NeighborEX"], ASS.sp)
            }
            each.ran <- cbind(min.r, max.r)
            colnames(each.ran) <- c("start", "end")
            g.1.p <- NULL
            t.sp <- c(g.1.s, g.2.s)
            e.chr <- unique(ea.re[,"Nchr"])
            pr.re <- ExReads(t.ex, t.sp, g.1.e, g.2.e, g.1.p,
                             g.2.p, g.1.s, g.2.s, e.chr, each.ran, "exon",
                             "ASS", ins, min.reads)
            pr.re
        }
        if (!any(length(splicingInfo)))
            return(NULL)
        SE.ratio <- NULL
        MXE.ratio <- NULL
        RI.ratio <- NULL
        ASS.ratio <- NULL
        if (!is.null(splicingInfo[["SE"]])) {
            print("Calculating PSI of SE events")
            SE.result <- rbind(splicingInfo$SE)
            SE.ratio <- bplapply(seq_len(nrow(SE.result)), SE.te,
                                 BPPARAM = parm)
            SE.ratio <- do.call(rbind, SE.ratio)
            rownames(SE.ratio) <- SE.result[,"EventID"]
        }
        if(!is.null(splicingInfo[["MXE"]])){
            print("Calculating PSI of MXE events")
            MXE.result <- rbind(splicingInfo$MXE)
            MXE.ratio <- suppressWarnings(bplapply(seq_len(nrow(MXE.result)),
                                                   MXE.te, BPPARAM = parm))
            MXE.ratio <- do.call(rbind, MXE.ratio)
            rownames(MXE.ratio) <- MXE.result[,"EventID"]
        }
        if (!is.null(splicingInfo[["RI"]])) {
            print("Calculating PSI of RI events")
            RI.result <- rbind(splicingInfo$RI)
            RI.ratio <- bplapply(seq_len(nrow(RI.result)), RI.te,
                                 BPPARAM = parm)
            RI.ratio <- do.call(rbind, RI.ratio)
            rownames(RI.ratio) <- RI.result[,"EventID"]
        }
        if (!is.null(splicingInfo[["A5SS"]]) |
            !is.null(splicingInfo[["A3SS"]])) {
            print("Calculating PSI of ASS events")
            ASS.result <- rbind(splicingInfo[["A5SS"]], splicingInfo[["A3SS"]])
            ASS.result <- cbind(ASS.result,
                                Types = c(rep("A5SS",
                                              nrow(splicingInfo[["A5SS"]])),
                                          rep("A3SS",
                                              nrow(splicingInfo[["A3SS"]]))))
            ASS.ratio <- bplapply(seq_len(nrow(ASS.result)), ASS.te,
                                  BPPARAM = parm)
            ASS.ratio <- do.call(rbind, ASS.ratio)
            rownames(ASS.ratio) <- ASS.result[,"EventID"]
        }
        final.re <- rbind(ASS.ratio, SE.ratio, MXE.ratio, RI.ratio)
        rowid <- suppressWarnings(rownames(final.re))
        final.re <- suppressWarnings(apply(final.re, 2, as.numeric))
        rownames(final.re) <- rowid
        return(final.re)
    }
    splitEnv <- environment(splitSplice)
    cEnv <- environment(coorEX)
    ins <- insert.size
    ea.re <- NULL
    MXE.num <- NULL
    SE.num <- NULL
    RI.num <- NULL
    ASS.num <- NULL
    g.1.p <- NULL
    g.2.p <- NULL
    fi.pr <- NULL
    se.pr <- NULL
    pr.re <- NULL
    fi.pair.2 <- NULL
    se.pair.2 <- NULL
    final.re <- NULL
    called.packages <- c("GenomicRanges", "GenomicFeatures")
    if(ncol(sample.info) >= 3){
        colnames(sample.info)[seq_len(3)] <- c("name", "path", "condition")
    } else{
        colnames(sample.info)[seq_len(2)] <- c("name", "path")
    }
    sample.info$name <- as.character(sample.info$name)
    sample.info$path <- as.character(sample.info$path)
    sample.files <- rbind(sample.info$path)
    sample.names <- sample.info$name
    parm <- SnowParam(workers = num.cores, type = "SOCK")
    final.ra <- Each.Cal.ratio(sample.files, AS.list, splitEnv, cEnv, ins,
                               min.reads, read.length, read.type, parm)
    colnames(final.ra) <- sample.names
    coldat <- data.frame(name = sample.info$name,
                         path = sample.info$path,
                         stringsAsFactors = FALSE)
    if("condition" %in% colnames(sample.info)){
        coldat$group <- sample.info$condition
    }
    final.res <- SummarizedExperiment(assays = list(psi = final.ra),
                                      colData = coldat)
    return(final.res)
}
