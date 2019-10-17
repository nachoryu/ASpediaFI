#' Detection of alternative splice sites
#'
#' This function detects alternative 5'/3' splice sites (ALSS) events.
#'
#' @param exon.info A matrix containing exon locus
#' @param intron.info A matrix containing intron locus
#' @param alt.intron.info A matrix containing flanking introns of
#' alternative exons
#' @return A matrix containing ALSS events
#' @details This function is borrowed from the \code{IVAS} package.
#' @references Han, S. et al. (2017). Genome wide discovery of genetic variants
#' affecting alternative splicing patterns in human using bioinformatics method.
#' \emph{Genes & Genomics}, 39.
#' @keywords internal
#' @noRd
detectALSS <- function(exon.info, intron.info, alt.intron.info) {

    alt.intron.info <- alt.intron.info[1, ]
    over.ex <- rbind(exon.info[exon.info[, "end"] > alt.intron.info["end"] +
                    1 & exon.info[, "start"] < alt.intron.info["start"] - 1, ])
    left.ex.num <- which(exon.info[, "end"] == alt.intron.info["start"] - 1)
    right.ex.num <- which(exon.info[, "start"] == alt.intron.info["end"] +
                            1)
    left.ex <- unique(cbind(as.double(rownames(exon.info)[left.ex.num]),
                            rbind(exon.info[left.ex.num, ])))
    right.ex <- unique(cbind(as.double(rownames(exon.info)[right.ex.num]),
                            rbind(exon.info[right.ex.num, ])))
    colnames(left.ex) <- c("id", "start", "end")
    colnames(right.ex) <- c("id", "start", "end")
    left.ex.result <- NULL
    for (i in seq_len(nrow(left.ex))) {
        lx <- intron.info[intron.info[, "start"] < left.ex[i, "start"] &
                            intron.info[, "end"] > left.ex[i, "end"], ]
        if (length(lx) != 0) {
            left.ex.result <- rbind(left.ex.result, cbind(left.ex[i, "id"],
                                                            rbind(lx)))
        }
    }
    right.ex.result <- NULL
    for (i in seq_len(nrow(right.ex))) {
        lx <- intron.info[intron.info[, "start"] < right.ex[i, "start"] &
                            intron.info[, "end"] > right.ex[i, "end"], ]
        if (length(lx) != 0) {
            right.ex.result <- rbind(right.ex.result, cbind(right.ex[i, "id"],
                                                            rbind(lx)))
        }
    }
    fi.te.num <- exon.info[, "start"] < alt.intron.info["start"] - 1 &
                    exon.info[, "end"] > alt.intron.info["start"] - 1
    fi.te.ex <- rbind(exon.info[fi.te.num, ])
    rownames(fi.te.ex) <- rownames(exon.info)[fi.te.num]
    fi.te.ex.num <- !is.element(paste(fi.te.ex[, "start"], fi.te.ex[, "end"]),
                                paste(over.ex[, "start"], over.ex[, "end"]))
    fi.nm <- rownames(fi.te.ex)[fi.te.ex.num]
    fi.te.ex <- rbind(fi.te.ex[fi.te.ex.num, ])
    rownames(fi.te.ex) <- fi.nm
    se.te.num <- intron.info[, "start"] < alt.intron.info["start"] &
                    intron.info[, "end"] > alt.intron.info["start"]
    se.te.int <- rbind(cbind(as.integer(rownames(intron.info)[se.te.num]),
                        rbind(intron.info[se.te.num, ])))
    se.te.test <- !is.element(paste(se.te.int[, 1], se.te.int[, "start"],
                                se.te.int[, "end"]),
                                paste(left.ex.result[, 1],
                                    left.ex.result[, "start"],
                                    left.ex.result[, "end"]))
    se.te.int <- rbind(se.te.int[se.te.test, ])
    th.te.num <- exon.info[, "start"] < alt.intron.info["end"] + 1 &
                    exon.info[, "end"] > alt.intron.info["end"] + 1
    th.te.ex <- rbind(exon.info[th.te.num, ])
    rownames(th.te.ex) <- rownames(exon.info)[th.te.num]
    th.te.ex.num <- !is.element(paste(th.te.ex[, "start"], th.te.ex[, "end"]),
                                paste(over.ex[, "start"], over.ex[, "end"]))
    th.rm <- rownames(th.te.ex)[th.te.ex.num]
    th.te.ex <- rbind(th.te.ex[th.te.ex.num, ])
    rownames(th.te.ex) <- th.rm
    fo.te.num <- intron.info[, "start"] < alt.intron.info["end"] &
                    intron.info[, "end"] > alt.intron.info["end"]
    fo.te.int <- rbind(cbind(as.integer(rownames(intron.info)[fo.te.num]),
                        rbind(intron.info[fo.te.num, ])))
    fo.te.test <- !is.element(paste(fo.te.int[, 1], fo.te.int[, "start"],
                                    fo.te.int[, "end"]),
                                paste(right.ex.result[, 1],
                                        right.ex.result[, "start"],
                                        right.ex.result[, "end"]))
    fo.te.int <- rbind(fo.te.int[fo.te.test, ])
    ALSS.final.result <- NULL
    if (length(fi.te.ex) != 0 | length(se.te.int) != 0) {
        tar.tx.num <- which(exon.info[, "end"] == alt.intron.info["start"] - 1)
        up.ex.start <- intron.info[is.element(intron.info[, "start"] - 1,
                                                exon.info[tar.tx.num, "end"]),
                                    "end"] + 1
        up.tx.exon.num <- which((is.element(rownames(exon.info),
                                            names(tar.tx.num)) &
                                    is.element(exon.info[, "start"],
                                                up.ex.start)) == "TRUE")
        tar.tx <- rbind(exon.info[tar.tx.num, ])
        up.tx <- rbind(exon.info[up.tx.exon.num, ])
        if (nrow(up.tx) == 1 | nrow(tar.tx) == 1) {
            rownames(up.tx) <- rownames(exon.info)[up.tx.exon.num]
            rownames(tar.tx) <- rownames(exon.info)[tar.tx.num]
        }
        each.up.tx <- tapply(up.tx[, "start"], rownames(up.tx), min)
        up.rn <- rownames(up.tx)[is.element(paste(rownames(up.tx),
                                                    up.tx[, "start"]),
                                            paste(names(each.up.tx),
                                                    each.up.tx))]
        up.tx <- rbind(up.tx[is.element(paste(rownames(up.tx),
                                                up.tx[, "start"]),
                                        paste(names(each.up.tx), each.up.tx)),
                                ])
        rownames(up.tx) <- up.rn
        over.names <- intersect(rownames(tar.tx), rownames(up.tx))
        tar.tx <- rbind(tar.tx[over.names, ])
        up.tx <- rbind(up.tx[over.names, ])
        tar.ex <- paste(tar.tx[, "start"], tar.tx[, "end"], sep = "-")
        tar.up.ex <- paste(up.tx[, "start"], up.tx[, "end"], sep = "-")
        ALSS.final.result <- NULL
        if (length(fi.te.ex) != 0) {
            long.tx <- rbind(fi.te.ex)
            long.int.num <- which(is.element(intron.info[, "start"] - 1,
                                                long.tx[, "end"]) == "TRUE")
            long.up.tx <- rbind(exon.info[is.element(rownames(exon.info),
                                    rownames(intron.info)[long.int.num]) &
                                        is.element(exon.info[, "start"],
                                    intron.info[long.int.num, "end"] + 1), ])
            if (nrow(long.up.tx) == 1) {
                rnforup <- rownames(exon.info)[is.element(
                                                rownames(exon.info),
                                                rownames(
                                                    intron.info)[long.int.num]
                                                ) &
                                                is.element(
                                                    exon.info[, "start"],
                                                    intron.info[long.int.num,
                                                                    "end"] + 1)]
                rownames(long.up.tx) <- rnforup
            }
            each.up.tx <- tapply(long.up.tx[, "start"], rownames(long.up.tx),
                                    min)
            long.up.tx <- rbind(long.up.tx[is.element(paste(
                                                        rownames(long.up.tx),
                                                        long.up.tx[, "start"]),
                                                        paste(names(each.up.tx),
                                                                each.up.tx)), ])
            rownames(long.up.tx) <- rownames(long.up.tx)[is.element(
                                        paste(rownames(long.up.tx),
                                            long.up.tx[, "start"]),
                                        paste(names(each.up.tx),
                                            each.up.tx))]
            over.names <- intersect(rownames(long.tx), rownames(long.up.tx))
            long.tx <- rbind(long.tx[over.names, ])
            long.up.tx <- rbind(long.up.tx[over.names, ])
            long.ex <- paste(long.tx[, "start"], long.tx[, "end"], sep = "-")
            long.up.ex <- paste(long.up.tx[, "start"], long.up.tx[, "end"],
                                sep = "-")
            if (length(long.up.ex) != 0) {
                for (i in seq_len(length(tar.ex))) {
                    for (j in seq_len(length(long.ex))) {
                        s.tar.ex <- unlist(strsplit(tar.ex[i], "-"))
                        s.tar.up.ex <- unlist(strsplit(tar.up.ex[i], "-"))
                        s.long.ex <- unlist(strsplit(long.ex[j], "-"))
                        if (s.tar.ex[2] >= s.long.ex[1] & s.tar.up.ex[1] <=
                            s.long.ex[2]) {
                                next
                        }
                        s.total.tar.ex <- do.call(rbind, strsplit(tar.ex, "-"))
                        s.total.tar.up.ex <- do.call(rbind, strsplit(tar.up.ex,
                                                                        "-"))
                        Re.rm <- as.double(s.total.tar.up.ex[, 1]) <=
                            s.long.ex[2] &
                                as.double(s.total.tar.ex[, 2]) >= s.long.ex[1]
                        tarex.nums <- grep(s.tar.ex[2], tar.ex)
                        tarex.des <- paste(sort(unique(
                                    tar.ex[is.element(seq_len(length(tar.ex)),
                                        tarex.nums) & !Re.rm])), collapse = ",")
                        tarex.up.des <- paste(sort(unique(
                                    tar.up.ex[is.element(
                                        seq_len(length(tar.ex)), tarex.nums) &
                                            !Re.rm])), collapse = ",")
                        if (length(which((tarex.nums & !Re.rm) == TRUE)) == 0) {
                            next
                        }
                        longex.nums <- grep(s.long.ex[2], long.ex)
                        longex.des <- paste(sort(unique(long.ex[longex.nums])),
                                            collapse = ",")
                        longex.up.des <- paste(
                                            sort(unique(
                                                    long.up.ex[longex.nums])),
                                            collapse = ",")
                        tar.mat <- cbind(tar.ex[i], long.ex[j], tar.up.ex[i],
                                            long.up.ex[j], tarex.des,
                                            longex.des, tarex.up.des,
                                            longex.up.des, "A5SS")
                        ALSS.final.result <- rbind(ALSS.final.result, tar.mat)
                    }
                }
            }
        }
        if (length(se.te.int) != 0) {
            se.te.tx <- rownames(intron.info)[is.element(
                                                paste(intron.info[, "start"],
                                                        intron.info[, "end"]),
                                                paste(se.te.int[, "start"],
                                                        se.te.int[, "end"]))]
            long.tx <- rbind(exon.info[is.element(rownames(exon.info),
                                                            se.te.tx) &
                                            is.element(exon.info[, "end"],
                                                se.te.int[, "start"] - 1), ])
            long.up.tx <- rbind(exon.info[is.element(rownames(exon.info),
                                                        se.te.tx) &
                                            is.element(exon.info[, "start"],
                                                    se.te.int[, "end"] + 1), ])
            if (nrow(long.tx) == 1 | nrow(long.up.tx) == 1) {
                rownames(long.tx) <- rownames(exon.info)[is.element(
                                                            rownames(exon.info),
                                                                se.te.tx) &
                                                            is.element(
                                                                exon.info[,
                                                                    "end"],
                                                                se.te.int[,
                                                                    "start"]
                                                                        - 1)]
                rownames(long.up.tx) <- rownames(exon.info)[is.element(
                                                                rownames(
                                                                    exon.info),
                                                                se.te.tx) &
                                                            is.element(
                                                                exon.info[,
                                                                    "start"],
                                                                se.te.int[,
                                                                    "end"] + 1)]
            }
            each.up.tx <- tapply(long.up.tx[, "start"], rownames(long.up.tx),
                                    min)
            long.up.tx <- rbind(long.up.tx[is.element(
                                            paste(rownames(long.up.tx),
                                                    long.up.tx[, "start"]),
                                            paste(names(each.up.tx),
                                                    each.up.tx)), ])
            rownames(long.up.tx) <- rownames(long.up.tx)[is.element(
                                        paste(rownames(long.up.tx),
                                            long.up.tx[, "start"]),
                                        paste(names(each.up.tx),
                                            each.up.tx))]
            over.names <- intersect(rownames(long.tx), rownames(long.up.tx))
            long.tx <- rbind(long.tx[over.names, ])
            long.up.tx <- rbind(long.up.tx[over.names, ])
            long.ex <- paste(long.tx[, "start"], long.tx[, "end"], sep = "-")
            long.up.ex <- paste(long.up.tx[, "start"], long.up.tx[, "end"],
                                    sep = "-")
            if (length(long.up.ex) != 0) {
                for (i in seq_len(length(tar.ex))) {
                    for (j in seq_len(length(long.ex))) {
                        if (as.double(strsplit(long.ex[j], "-")[[1]][2]) >
                            as.double(strsplit(tar.ex[i], "-")[[1]][1])) {
                                s.tar.ex <- as.double(unlist(strsplit(tar.ex[i],
                                                                        "-")))
                                s.long.ex <- as.double(unlist(strsplit(
                                                                long.ex[j],
                                                                    "-")))
                                s.long.up.ex <- unlist(strsplit(long.up.ex[j],
                                                                    "-"))
                            if (s.long.up.ex[1] <= s.tar.ex[2] & s.long.ex[2] >=
                                s.tar.ex[1]) {
                                next
                            }
                            s.total.long.ex <- do.call(rbind, strsplit(long.ex,
                                                                        "-"))
                            s.total.long.up.ex <- do.call(rbind,
                                                            strsplit(long.up.ex,
                                                                        "-"))
                            Re.rm <- as.double(s.total.long.ex[, 2]) >=
                                                    s.tar.ex[1] &
                                        as.double(s.total.long.up.ex[, 1]) <=
                                                    s.tar.ex[2]
                            tarex.nums <- grep(s.tar.ex[2], tar.ex)
                            tarex.des <- paste(sort(unique(tar.ex[tarex.nums])),
                                                collapse = ",")
                            tarex.up.des <- paste(
                                                sort(unique(
                                                    tar.up.ex[tarex.nums])),
                                                collapse = ",")
                            longex.nums <- grep(s.long.ex[2], long.ex)
                            longex.des <- paste(
                                sort(unique(long.ex[is.element(
                                    seq_len(length(long.ex)), longex.nums) &
                                        !Re.rm])), collapse = ",")
                            longex.up.des <- paste(
                                sort(unique(long.up.ex[is.element(
                                    seq_len(length(long.ex)), longex.nums) &
                                        !Re.rm])), collapse = ",")
                            if (length(which((is.element(
                                    seq_len(length(long.ex)), longex.nums) &
                                        !Re.rm) == TRUE)) == 0) {
                                next
                            }
                            tar.mat <- cbind(long.ex[j], tar.ex[i],
                                                long.up.ex[j], tar.up.ex[i],
                                                longex.des, tarex.des,
                                                longex.up.des,
                                                tarex.up.des, "A5SS")
                            ALSS.final.result <- rbind(ALSS.final.result,
                                                        tar.mat)
                        }
                    }
                }
            }
        }
    }
    if (length(th.te.ex) != 0 | length(fo.te.int) != 0) {
        tar.tx.num <- which(exon.info[, "start"] == alt.intron.info["end"] + 1)
        down.ex.start <- intron.info[is.element(intron.info[, "end"] + 1,
                                                    exon.info[tar.tx.num,
                                                                "start"]),
                                    "start"] - 1
        down.tx.exon.num <- which((is.element(rownames(exon.info),
                                    names(tar.tx.num)) &
                                    is.element(exon.info[, "end"],
                                        down.ex.start)) == "TRUE")
        tar.tx <- rbind(exon.info[tar.tx.num, ])
        down.tx <- rbind(exon.info[down.tx.exon.num, ])
        if (nrow(down.tx) == 1 | nrow(tar.tx) == 1) {
            rownames(down.tx) <- rownames(exon.info)[down.tx.exon.num]
            rownames(tar.tx) <- rownames(exon.info)[tar.tx.num]
        }
        each.down.tx <- tapply(down.tx[, "start"], rownames(down.tx), max)
        down.rn <- rownames(down.tx)[is.element(paste(rownames(down.tx),
                                                        down.tx[, "start"]),
                                                paste(names(each.down.tx),
                                                        each.down.tx))]
        down.tx <- rbind(down.tx[is.element(paste(rownames(down.tx),
                                                    down.tx[, "start"]),
                                            paste(names(each.down.tx),
                                                    each.down.tx)), ])
        rownames(down.tx) <- down.rn
        over.names <- intersect(rownames(tar.tx), rownames(down.tx))
        tar.tx <- rbind(tar.tx[over.names, ])
        down.tx <- rbind(down.tx[over.names, ])
        tar.ex <- paste(tar.tx[, "start"], tar.tx[, "end"], sep = "-")
        tar.down.ex <- paste(down.tx[, "start"], down.tx[, "end"], sep = "-")
        if (length(th.te.ex) != 0) {
            long.tx <- rbind(th.te.ex)
            long.int.num <- which(is.element(intron.info[, "end"] + 1,
                                                long.tx[, "start"]) == "TRUE")
            long.down.num <- is.element(rownames(exon.info),
                                        rownames(intron.info)[long.int.num]) &
                                is.element(exon.info[, "end"],
                                    intron.info[long.int.num, "start"] - 1)
            long.down.tx <- rbind(exon.info[long.down.num, ])
            rownames(long.down.tx) <- rownames(exon.info)[long.down.num]
            each.down.tx <- tapply(long.down.tx[, "start"],
                                    rownames(long.down.tx), max)
            long.down.rn <- rownames(long.down.tx)[is.element(
                                paste(rownames(long.down.tx),
                                    long.down.tx[, "start"]),
                                paste(names(each.down.tx),
                                    each.down.tx))]
            long.down.tx <- rbind(long.down.tx[is.element(
                                paste(rownames(long.down.tx),
                                    long.down.tx[, "start"]),
                                paste(names(each.down.tx),
                                    each.down.tx)), ])
            rownames(long.down.tx) <- long.down.rn
            over.names <- intersect(rownames(long.tx), rownames(long.down.tx))
            long.tx <- rbind(long.tx[over.names, ])
            long.down.tx <- rbind(long.down.tx[over.names, ])
            long.ex <- paste(long.tx[, "start"], long.tx[, "end"], sep = "-")
            long.down.ex <- paste(long.down.tx[, "start"],
                                    long.down.tx[, "end"],
                                    sep = "-")
            if (length(long.down.ex) != 0) {
                for (i in seq_len(length(tar.ex))) {
                    for (j in seq_len(length(long.ex))) {
                        s.tar.ex <- as.double(unlist(strsplit(tar.ex[i], "-")))
                        s.long.ex <- as.double(unlist(strsplit(long.ex[j],
                                                                "-")))
                        s.tar.down.ex <-
                                as.double(unlist(strsplit(tar.down.ex[i],
                                                            "-")))
                        if (s.tar.down.ex[2] >= s.long.ex[1] & s.tar.ex[1] <=
                            s.long.ex[2]) {
                            next
                        }
                        s.total.tar.down.ex <- do.call(rbind,
                                                        strsplit(tar.down.ex,
                                                                        "-"))
                        s.total.tar.ex <- do.call(rbind, strsplit(tar.ex, "-"))
                        Re.rm <- as.double(s.total.tar.down.ex[, 2]) >=
                                            s.long.ex[1] &
                                    as.double(s.total.tar.ex[, 1]) <=
                                                s.long.ex[2]
                        tarex.nums <- grep(s.tar.ex[1], tar.ex)
                        tarex.des <- paste(
                                        sort(unique(
                                            tar.ex[is.element(
                                                seq_len(length(tar.ex)),
                                                            tarex.nums) &
                                                        !Re.rm])),
                                        collapse = ",")
                        tarex.down.des <- paste(
                                        sort(unique(
                                            tar.down.ex[is.element(
                                                seq_len(length(tar.ex)),
                                                    tarex.nums) & !Re.rm])),
                                        collapse = ",")
                        if (length(which((tarex.nums & !Re.rm) == TRUE)) == 0) {
                            next
                        }
                        longex.nums <- grep(s.long.ex[1], long.ex)
                        longex.des <- paste(sort(unique(long.ex[longex.nums])),
                                            collapse = ",")
                        longex.down.des <- paste(sort(unique(
                                                    long.down.ex[longex.nums])),
                                                    collapse = ",")
                        tar.mat <- cbind(tar.ex[i], long.ex[j], tar.down.ex[i],
                                            long.down.ex[j], tarex.des,
                                            longex.des, tarex.down.des,
                                            longex.down.des, "A3SS")
                        ALSS.final.result <- rbind(ALSS.final.result, tar.mat)
                    }
                }
            }
        }
        if (length(fo.te.int) != 0) {
            fo.te.tx <- rownames(intron.info)[is.element(
                                                paste(intron.info[, "start"],
                                                        intron.info[, "end"]),
                                                paste(fo.te.int[, "start"],
                                                        fo.te.int[, "end"]))]
            long.tx <- rbind(exon.info[is.element(rownames(exon.info),
                                                    fo.te.tx) &
                                        is.element(exon.info[, "start"],
                                                    fo.te.int[, "end"] + 1), ])
            long.down.tx <- rbind(exon.info[is.element(rownames(exon.info),
                                                        fo.te.tx) &
                                            is.element(exon.info[, "end"],
                                                        fo.te.int[, "start"] -
                                                                        1), ])
            if (nrow(long.down.tx) == 1 | nrow(long.tx) == 1) {
                rownames(long.down.tx) <- rownames(
                                            exon.info)[is.element(
                                                rownames(exon.info), fo.te.tx) &
                                                        is.element(exon.info[,
                                                                        "end"],
                                                            fo.te.int[,
                                                                "start"] - 1)]
                rownames(long.tx) <- rownames(
                                        exon.info)[is.element(
                                            rownames(exon.info), fo.te.tx) &
                                                is.element(exon.info[, "start"],
                                                        fo.te.int[, "end"] + 1)]
            }
            each.down.tx <- tapply(long.down.tx[, "start"],
                                    rownames(long.down.tx), max)
            long.down.tx <- rbind(long.down.tx[is.element(
                                                paste(rownames(long.down.tx),
                                                    long.down.tx[, "start"]),
                                                paste(names(each.down.tx),
                                                    each.down.tx)), ])
            rownames(long.down.tx) <- rownames(long.down.tx)[is.element(
                                                paste(rownames(long.down.tx),
                                                    long.down.tx[, "start"]),
                                                paste(names(each.down.tx),
                                                    each.down.tx))]
            over.names <- intersect(rownames(long.tx), rownames(long.down.tx))
            long.tx <- rbind(long.tx[over.names, ])
            long.down.tx <- rbind(long.down.tx[over.names, ])
            long.ex <- paste(long.tx[, "start"], long.tx[, "end"], sep = "-")
            long.down.ex <- paste(long.down.tx[, "start"],
                                    long.down.tx[, "end"], sep = "-")
            if (length(long.down.ex) != 0) {
                for (i in seq_len(length(tar.ex))) {
                    for (j in seq_len(length(long.ex))) {
                        if (as.double(strsplit(long.ex[j], "-")[[1]][1]) <
                            as.double(strsplit(tar.ex[i], "-")[[1]][2])) {
                            s.tar.ex <- as.double(unlist(strsplit(tar.ex[i],
                                                                        "-")))
                            tarex.nums <- grep(s.tar.ex[1], tar.ex)
                            tarex.des <- paste(sort(unique(tar.ex[tarex.nums])),
                                                collapse = ",")
                            tarex.down.des <- paste(sort(unique(
                                                    tar.down.ex[tarex.nums])),
                                                collapse = ",")
                            s.long.ex <- as.double(unlist(strsplit(long.ex[j],
                                                                        "-")))
                            s.total.long.ex <- do.call(rbind, strsplit(long.ex,
                                                                        "-"))
                            s.down.long.ex <- as.double(unlist(strsplit(
                                                            long.down.ex[j],
                                                                        "-")))
                            s.total.down.long.ex <- do.call(rbind, strsplit(
                                                            long.down.ex,
                                                                        "-"))
                            if (s.down.long.ex[2] >= s.tar.ex[1] &
                                s.long.ex[1] <=
                                s.tar.ex[2]) {
                                next
                            }
                            Re.rm <- as.double(s.total.down.long.ex[, 2]) >=
                                s.tar.ex[1] &
                                as.double(s.total.long.ex[, 1]) <= s.tar.ex[2]
                            longex.nums <- grep(s.long.ex[1], long.ex)
                            longex.des <- paste(sort(unique(
                                long.ex[is.element(seq_len(length(long.ex)),
                                                        longex.nums) &
                                                                    !Re.rm])),
                                                    collapse = ",")
                            longex.down.des <- paste(sort(unique(
                                                    long.down.ex[is.element(
                                                        seq_len(length(
                                                            long.ex)),
                                                    longex.nums) & !Re.rm])),
                                                        collapse = ",")
                            if (length(which((is.element(seq_len(
                                                            length(long.ex)),
                                                            longex.nums) &
                                                !Re.rm) == TRUE)) == 0) {
                                next
                            }
                            tar.mat <- cbind(long.ex[j], tar.ex[i],
                                                long.down.ex[j],
                                                tar.down.ex[i], longex.des,
                                                tarex.des, longex.down.des,
                                                tarex.down.des, "A3SS")
                            ALSS.final.result <- rbind(ALSS.final.result,
                                                        tar.mat)
                        }
                    }
                }
            }
        }
    }
    if (length(ALSS.final.result) == 0) {
        return(NULL)
    }
    ALSS.final.result <- unique(ALSS.final.result)
    ALSS.final.result <- cbind(ALSS.final.result, "possible")
    colnames(ALSS.final.result) <- c("ShortEX", "LongEX", "ShortNeighborEX",
                                        "LongNeighborEX", "Short_des",
                                        "Long_des", "ShortNeighbor_des",
                                        "LongNeighbor_des",
                                        "Types", "status")
    A3.num <- which(ALSS.final.result[, "Types"] == "A3SS")
    A5.num <- which(ALSS.final.result[, "Types"] == "A5SS")
    A3.final.result <- rbind(ALSS.final.result[A3.num, ])
    A5.final.result <- rbind(ALSS.final.result[A5.num, ])
    int.test <- paste(intron.info[, "start"] - 1, intron.info[, "end"] + 1,
                        sep = "-")
    if (length(A3.num) != 0) {
        A3.long.test.int <- paste(do.call(rbind,
                                    strsplit(A3.final.result[,
                                                "ShortNeighborEX"], "-"))[, 2],
                                    do.call(rbind,
                                        strsplit(A3.final.result[, "ShortEX"],
                                                        "-"))[, 1], sep = "-")
        A3.short.test.int <- paste(do.call(rbind,
                                    strsplit(A3.final.result[,
                                                "LongNeighborEX"], "-"))[, 2],
                                    do.call(rbind,
                                        strsplit(A3.final.result[, "LongEX"],
                                                    "-"))[, 1], sep = "-")
        A3.final.result[is.element(A3.long.test.int, int.test) |
                            is.element(A3.short.test.int, int.test),
                        "status"] <- "exist"
    }
    if (length(A5.num) != 0) {
        A5.long.test.int <- paste(do.call(rbind,
                                    strsplit(A5.final.result[,
                                        "ShortEX"], "-"))[, 2],
                                    do.call(rbind,
                                        strsplit(A5.final.result[,
                                            "ShortNeighborEX"], "-"))[, 1],
                                                sep = "-")
        A5.short.test.int <- paste(do.call(rbind,
                                    strsplit(A5.final.result[, "LongEX"],
                                                "-"))[, 2],
                                    do.call(rbind,
                                        strsplit(A5.final.result[,
                                            "LongNeighborEX"], "-"))[, 1],
                                                sep = "-")
        A5.final.result[is.element(A5.long.test.int, int.test) |
                            is.element(A5.short.test.int, int.test),
                            "status"] <- "exist"
    }
    ALSS.final.result <- rbind(A5.final.result, A3.final.result)
    ALSS.final.result <- rbind(ALSS.final.result[!is.na(
                                ALSS.final.result[, "ShortNeighborEX"]) &
                                    !is.na(
                                        ALSS.final.result[, "LongNeighborEX"]),
                                                    ])
    ALSS.final.result <- rbind(ALSS.final.result[, c("ShortEX", "LongEX",
                                                        "ShortNeighborEX",
                                                        "LongNeighborEX",
                                                        "Types", "status")])
    return(unique(ALSS.final.result))
}
