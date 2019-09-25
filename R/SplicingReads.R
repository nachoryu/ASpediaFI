#'Counting junction and paired-end reads
#'
#'This function counts junction and paired-end reads.
#'
#'@param bam.file a path to bam file
#'@param test.exon a data frame containing AS target exon and its neighbors
#'@param spli.jun A data frame containing junction information
#'@param e.ran A range of parsingr eads from a bam file
#'@param chr A chromosome number
#'@param read.type a type of RNA-Seq reads
#'@param insert.size An insert.sizert size
#'@return A list of counts of junction reads and paired-end reads
#'@details This function is borrowed from the \code{IMAS} package.
#'@references Han, S. et al. (2017). IMAS: Integrative analysis of Multi-omics
#'data for Alternative Splicing. R package version 1.8.0.
#'@import GenomicAlignments
#'@importFrom GenomicRanges GRanges
#'@importFrom IRanges IRanges
#'@importFrom Rsamtools ScanBamParam
#'@importFrom S4Vectors elementMetadata
#'@keywords internal
SplicingReads <- function(bam.file = NULL, test.exon = NULL, spli.jun = NULL,
                          e.ran = NULL, chr = NULL, read.type = "paired",
                          insert.size = 40){
    ReadCover <- function(t.reads) {
        cigar.s <- cigar(t.reads)
        pos.s <- start(t.reads)
        ex.reads <- extractAlignmentRangesOnReference(cigar.s, pos.s)
        names(ex.reads) <- names(t.reads)
        return(ex.reads)
    }
    CoverEnv <- environment(ReadCover)
    N.test <- function(Nreads){
        p.se.nms <- paste(names(Nreads), elementMetadata(Nreads)[,"seq"])
        t.p.se.nms <- table(p.se.nms)
        nums.t <- which(t.p.se.nms == 1)
        nm.t <- names(t.p.se.nms)[nums.t]
        nm.t <- which(is.element(t.p.se.nms, nums.t) == TRUE)
        Nreads <- Nreads[nm.t, ]
        Nco <- CoverEnv$ReadCover(Nreads)
        spliInfo <- lapply(Nco, function(ea.re) {
            st <- sort(start(ea.re))[-1]
            en <- sort(end(ea.re))[-length(ea.re)]
            paste(en, st, sep = "-")
        })
        spliInfo <- do.call(rbind, spliInfo)
        cu.info <- is.element(spliInfo, spli.jun)
        spliInfo <- spliInfo[cu.info, ]
        t.spliInfo <- table(spliInfo)
        spliInfo <- as.integer(t.spliInfo)
        names(spliInfo) <- names(t.spliInfo)
        f.li <- list(names(spliInfo), spliInfo)
        names(f.li) <- c("re.nms", "junctionInfo")
        return(f.li)
    }
    ex.test <- function(Pco.Gran, over.mat, pre.re) {
        Tnms <- names(Pco.Gran)[over.mat[, "queryHits"]]
        Tnms <- cbind(Tnms, nums = over.mat[, "queryHits"])
        E.reads <- rbind(Tnms[!is.element(Tnms[, "Tnms"], unlist(pre.re)),])
        ov.E.r <- is.element(over.mat[, "queryHits"], E.reads[,"nums"])
        ov.E.r <- table(over.mat[ov.E.r, "subjectHits"])
        names(ov.E.r) <- p.ex[as.integer(names(ov.E.r))]
        return(ov.E.r)
    }
    exEnv <- environment(ex.test)
    pa.test <- function(Preads, read.type) {
        count.reads <- function(bet.num) {
            fi.nums <- over.mat[,"subjectHits"] == bet.num[1]
            se.nums <- over.mat[,"subjectHits"] == bet.num[2]
            f.mat <- rbind(over.mat[fi.nums | se.nums,])
            f.nms <- names(Pco.Gran)[f.mat[,"queryHits"]]
            t.f.nms <- table(f.nms)
            t.f.nms <- t.f.nms[t.f.nms == 2]
            if (!any(seq_along(t.f.nms)))
                return(NULL)
            return(names(t.f.nms))
        }
        count.op <- function(te.li) {
            names(te.li) <- c("In", "Sk")
            in1.n <- as.integer(te.li$In[1,])
            in2.n <- as.integer(te.li$In[2,])
            sk.n <- as.integer(te.li$Sk[1,])
            in1.nm <- paste(p.ex[in1.n], collapse = "~")
            in2.nm <- paste(p.ex[in2.n], collapse = "~")
            sk.nm <- paste(p.ex[sk.n], collapse = "~")
            in1 <- count.reads(in1.n)
            in2 <- count.reads(in2.n)
            sk <- count.reads(sk.n)
            f.sk <- sk[!is.element(sk, c(in1, in2))]
            f.li <- list(in1, in2, sk)
            names(f.li) <- c(in1.nm, in2.nm, sk.nm)
            return(f.li)
        }
        Pco <- CoverEnv$ReadCover(Preads)
        Pco.Gran <- GRanges("*", unlist(Pco))
        over.mat <- findOverlaps(Pco.Gran, Ex.range)
        over.mat <- as.matrix(over.mat)
        if (read.type == "exon") {
            ex.re <- exEnv$ex.test(Pco.Gran, over.mat, NULL)
            return(ex.re)
        }
        in.mat <- rbind(c(1, 2), c(2, 3))
        sk.mat <- rbind(c(1, 3))
        pre.re <- count.op(list(in.mat, sk.mat))
        if (length(Ex.range) == 4) {
            in.mat <- rbind(c(2, 3), c(3, 4))
            sk.mat <- rbind(c(2, 4))
            pre.re2 <- count.op(list(in.mat, sk.mat))
            Nov <- !is.element(names(pre.re2), names(pre.re))
            pre.re2 <- pre.re2[Nov]
            pre.re3 <- count.op(list(NULL, rbind(c(1, 4))))
            pre.re3 <- pre.re3[length(pre.re3)]
            pre.re <- apply(cbind(pre.re, pre.re2, pre.re3),
                            1, function(x) unname(unlist(x)))
        }
        ex.re <- exEnv$ex.test(Pco.Gran, over.mat, pre.re)
        pre.re <- sapply(pre.re, function(x) length(x))
        f.re <- list(pre.re, ex.re)
        names(f.re) <- c("pairedInfo", "exonInfo")
        return(f.re)
    }
    final.exon.result <- NULL
    if (!any(ncol(test.exon)) & any(length(test.exon))) {
        test.exon <- do.call(rbind, strsplit(test.exon, "-"))
        colnames(test.exon) <- c("start", "end")
        st <- as.double(test.exon[,"start"])
        en <- as.double(test.exon[,"end"])
        Ex.range <- GRanges("*", IRanges(st, en))
        p.ex <- paste(test.exon[,"start"], test.exon[,"end"], sep = "-")
    }
    if (!any(ncol(spli.jun)) & any(length(spli.jun))) {
        total.spli.sites <- unique(do.call(rbind, strsplit(spli.jun, "-")))
        colnames(total.spli.sites) <- c("start", "end")
    }
    final.re <- list(NULL, NULL, NULL)
    names(final.re) <- c("pairedInfo", "exonInfo", "junctionInfo")
    I.ran <- IRanges(min(as.integer(e.ran)), max(as.integer(e.ran)))
    which.ragnes <- GRanges(seqnames = chr, ranges = I.ran)
    what.param <- c("qname", "rname", "pos", "mapq", "cigar",
                    "seq", "qual")
    param <- ScanBamParam(which = which.ragnes, what = what.param, tag = "MD")
    T.reads <- readGAlignments(bam.file, use.names = TRUE, param = param)
    if (length(T.reads) < 10)
        return(final.re)
    Nreads <- T.reads[njunc(T.reads) == 1, ]
    N.re <- N.test(Nreads)
    not.N <- !is.element(names(T.reads), N.re$re.nms)
    Not.N.reads <- T.reads[not.N, ]
    P.re <- pa.test(Not.N.reads, read.type)
    N.re <- N.re$junctionInfo
    if(read.type == "paired"){
        final.re <- list(P.re$pairedInfo, P.re$exonInfo, N.re)
    }
    else{
        final.re <- list(NULL, P.re, N.re)
    }
    names(final.re) <- c("pairedInfo", "exonInfo", "junctionInfo")
    return(final.re)
}
