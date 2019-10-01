#'Discriminative Random Walk with Restart (DRaWR)
#'
#'This function runs the DRaWR on query gene set.
#'
#'@param query.genes a data frame containing query genes and their weights
#'@param universe a data frame containing entire genes
#'@param network a data frame containing edges and their weights and types
#'@param restart a restart probability
#'@param num.folds number of folds for cross-validation
#'@param num.feats number of feature nodes to be kept in the final network
#'@return A list of AS events
#'@details This function is modified from the \code{DRaWR} function in the
#'\code{DRaWR} package.
#'@importFrom ROCR performance prediction plot
#'@importFrom Matrix t colSums
#'@importFrom DRaWR threeCol2listMat threeCol2MaxMat RWR
#'@importFrom stats aggregate
#'@importFrom graphics legend
#'@importFrom methods new slot
#'@references Blatti, C. et al. (2016). Characterizing gene sets using
#'discrminative random walks with restart on heterogeneous biological networks
#'\emph{Bioinformatics}, 32.
#'@keywords internal
drawr <- function (query.genes, universe, network, restart,
                   num.folds, num.feats) {

    maxiters = 50
    thresh = 1e-04
    edges = network
    colnames(edges) = c("src", "target", "weight", "type")
    edges$weight <- as.numeric(as.character(edges$weight))
    wmin = min(edges$weight)
    wmins = edges[which(edges$weight == wmin), ]
    if (wmin <= 0) {
        show("Invalid Edge Weights")
        return(-1)
    }
    wmax = max(edges$weight)
    wmaxs = edges[which(edges$weight == wmax), ]
    property_types = c("AS", "Pathway")
    all_etypes = as.character(unique(edges$type))
    prop_etypes = intersect(all_etypes, property_types)
    ntypes = length(prop_etypes)
    features = unique(edges[which(edges$type %in% prop_etypes),
                            c("src", "type")])
    featnodes = as.character(features[, "src"])
    nfeats = length(featnodes)
    nkeep = min(num.feats, nfeats)

    #Make network undirected
    restable = NULL
    tmpsrc = c(as.character(edges$src), as.character(edges$target))
    tmptarget = c(as.character(edges$target), as.character(edges$src))
    tmpweight = c(edges$weight, edges$weight)
    tmptype = c(as.character(edges$type), as.character(edges$type))
    edges = data.frame(tmpsrc, tmptarget, tmpweight, tmptype,
                       stringsAsFactors = FALSE)
    colnames(edges) = c("src", "target", "weight", "type")

    #Normalize edge type
    typetable = aggregate(weight ~ type, data = edges, FUN = sum)
    rownames(typetable) = as.character(typetable[, 1])
    edges$weight = edges$weight/typetable[as.character(edges$type),
                                          "weight"]

    forcelarge = 0
    node_estimate = length(unique(c(as.character(edges[, 1]),
                                    as.character(edges[, 2]))))
    boolSparceMat = (node_estimate^2 < .Machine$integer.max) *
        (1 - forcelarge)
    transmat = NULL
    nodenames = NULL
    colsum = NULL
    if (boolSparceMat) {
        n1Matrix = threeCol2MaxMat(as.character(edges[,"src"]),
                                   as.character(edges[,"target"]),
                                   as.numeric(edges[,"weight"]))
        nodenames = rownames(n1Matrix)
        colsum = colSums(n1Matrix)
        transmat = t(t(n1Matrix)/colsum)
        rm(n1Matrix)
    }
    else {
        ll = threeCol2listMat(as.character(edges[,"src"]),
                              as.character(edges[,"target"]),
                              as.numeric(edges[,"weight"]))
        transmat = ll
        nodenames = ll$avals
        colsum = ll$colsums
        rm(ll)
    }
    nnodes = length(nodenames)
    rownames(universe) = as.character(universe[,1])
    if (dim(universe)[2] < 2) {
        universe = cbind(as.character(universe[,1]),
                         rep(1, length(universe[,1])))
    }
    uniIDs = sort(intersect(nodenames, as.character(unique(universe[,1]))))
    if (length(uniIDs) < 1) {
        return(-1)
    }

    restable = NULL

    evaltabstart = cbind(nodenames, -1, 0)
    colnames(evaltabstart) = c("node", "type","universe")
    rownames(evaltabstart) = nodenames
    evaltabstart[uniIDs, "universe"] = 1
    evaltabstart[as.character(features[, 1]), "type"] = as.character(features[,
                                                                              2])
    blankvec = structure(rep(0, nnodes), names = nodenames)
    startvec = blankvec + 1/nnodes
    biter = 0

    #Baseline
    query = blankvec
    midxs = match(uniIDs, nodenames)
    query[midxs] = as.numeric(universe[uniIDs, 2])
    rwr_res = RWR(boolSparceMat, transmat, restart, query, startvec,
                  maxiters, thresh)
    biter = rwr_res$iter
    evaltabstart = cbind(evaltabstart, as.numeric(rwr_res$vec))
    colnames(evaltabstart)[4] = "baseline"

    rownames(query.genes) = as.character(query.genes[, 1])
    queryIDs = sort(intersect(uniIDs, as.character(unique(query.genes[, 1]))))
    nquery = length(queryIDs)

    #Cross-validation
    folds = sample(cut(seq(1, nquery), breaks = num.folds, labels = FALSE))

    perf_b <- new("performance", x.name = "False positive rate",
                  y.name = "True positive rate", alpha.name = "Cutoff")
    perf_s1 <- new("performance", x.name = "False positive rate",
                   y.name = "True positive rate", alpha.name = "Cutoff")
    perf_d <- new("performance", x.name = "False positive rate",
                  y.name = "True positive rate", alpha.name = "Cutoff")
    perf_s2 <- new("performance", x.name = "False positive rate",
                   y.name = "True positive rate", alpha.name = "Cutoff")

    for (iter in 1:num.folds) {
        train_idxs = which(folds != iter)
        train_nidxs = queryIDs[train_idxs]
        ntrain = length(train_nidxs)
        test_idxs = which(folds == iter)
        testuni = as.character(setdiff(uniIDs, train_nidxs))
        if (num.folds == 1) {
            test_idxs = train_idxs
            testuni = as.character(uniIDs)
        }
        test_nidxs = queryIDs[test_idxs]
        ntest = length(test_nidxs)
        ntestuni = length(testuni)
        if (ntrain * ntest * 1 * ntestuni == 0) {
            row1 = c(iter, "baseline", NA)
            row2 = c(iter, "baseline", NA)
            row3 = c(iter, "baseline", NA)
            row4 = c(iter, "baseline", NA)
            restable = rbind(restable, row1, row2, row3, row4)
        }
        query = blankvec
        midxs = match(train_nidxs, nodenames)
        query[midxs] = as.numeric(query.genes[train_nidxs, 2])
        rwr_res = RWR(boolSparceMat, transmat, restart, query, startvec,
                      maxiters, thresh)
        qiter = rwr_res$iter
        evaltab = cbind(evaltabstart, 0, 0, as.numeric(rwr_res$vec[nodenames]))
        colnames(evaltab) = c("node", "type", "universe", "baseline", "train",
                              "test", "stage1")
        diff = as.numeric(evaltab[,"stage1"]) - as.numeric(evaltab[,"baseline"])
        evaltab = cbind(evaltab, diff)
        evaltab[train_nidxs, "train"] = 1
        evaltab[test_nidxs, "test"] = 1
        model = prediction(as.numeric(evaltab[testuni, "baseline"]),
                           evaltab[testuni, "test"])
        auc = performance(model, "auc")
        perf_baseline = performance(model, "tpr", "fpr")
        aucval = round(as.numeric(slot(auc, "y.values")), 3)
        row = c(iter, "baseline", aucval)
        restable = rbind(restable, row)
        model = prediction(as.numeric(evaltab[testuni, "stage1"]),
                                 evaltab[testuni, "test"])
        auc = performance(model, "auc")
        perf_stage1 = performance(model, "tpr", "fpr")
        aucval = round(as.numeric(slot(auc, "y.values")), 3)
        row = c(iter, "stage1", aucval)
        restable = rbind(restable, row)
        model = prediction(as.numeric(evaltab[testuni, "diff"]),
                                 evaltab[testuni, "test"])
        auc = performance(model, "auc")
        perf_diff = performance(model, "tpr", "fpr")
        aucval = round(as.numeric(slot(auc, "y.values")), 3)
        row = c(iter, "diff", aucval)
        restable = rbind(restable, row)
        colnames(restable) = c("iter", "stage", "aucval")
        if (nfeats == 0) {
            row4 = c(iter, "stage2", NA)
            restable = rbind(restable, row4)
        }
        keep = rep(-1, nnodes)
        evaltab = cbind(evaltab, keep)
        ss = sort(as.numeric(evaltab[featnodes, "diff"]), decreasing = TRUE,
                  index.return = TRUE)
        sortedfeats = featnodes[ss$ix]
        keepfeats = sortedfeats[1:nkeep]
        evaltab[featnodes, "keep"] = 0
        evaltab[keepfeats, "keep"] = 1
        newedges = edges[which(edges$type %in% setdiff(all_etypes,
                                                       prop_etypes)), ]
        keepidxs = which(edges[, "src"] %in% keepfeats |
                             edges[, "target"] %in% keepfeats)
        newedges = rbind(newedges, edges[keepidxs, ])
        tmpnames = unique(c(as.character(newedges[, 1]),
                            as.character(newedges[, 2])))
        train_nidxs2 = intersect(train_nidxs, tmpnames)
        test_nidxs2 = intersect(test_nidxs, tmpnames)
        testuni2 = intersect(testuni, tmpnames)
        if (length(train_nidxs2) * length(test_nidxs2) *
            1 * length(testuni2) == 0) {
            row4 = c(iter, "stage2", NA)
            restable = rbind(restable, row4)
        }
        typetable2 = NULL
        typetable2 = aggregate(weight ~ type, data = newedges,
                               FUN = sum)
        rownames(typetable2) = as.character(typetable2[,1])
        newedges$weight = newedges$weight/typetable2[as.character(newedges$type),
                                                     "weight"]
        transmat2 = NULL
        nodenames2 = NULL
        colsum2 = NULL
        boolSparceMat2 = (length(tmpnames)^2 < .Machine$integer.max) *
            (1 - forcelarge)
        if (boolSparceMat2) {
            n2Matrix = threeCol2MaxMat(as.character(newedges[,"src"]),
                                       as.character(newedges[,"target"]),
                                       as.numeric(newedges[,"weight"]))
            nodenames2 = rownames(n2Matrix)
            colsum2 = colSums(n2Matrix)
            transmat2 = t(t(n2Matrix)/colsum2)
            rm(n2Matrix)
        }
        else {
            ll2 = threeCol2listMat(as.character(newedges[,"src"]),
                                   as.character(newedges[,"target"]),
                                   as.numeric(newedges[,"weight"]))
            transmat2 = ll2
            nodenames2 = ll2$avals
            colsum2 = ll2$colsums
            rm(ll2)
        }
        nnodes2 = length(nodenames2)
        rm(newedges)
        blankvec2 = structure(rep(0, nnodes2), names = nodenames2)
        query2 = blankvec2
        midxs = match(train_nidxs2, nodenames2)
        query2[midxs] = as.numeric(query.genes[train_nidxs2, 2])
        startvec2 = blankvec2 + 1/nnodes2
        rwr_res = RWR(boolSparceMat2, transmat2, restart,
                      query2, startvec2, maxiters, thresh)
        q2iter = rwr_res$iter
        stage2 = rep(0, nnodes)
        evaltab = cbind(evaltab, stage2)
        evaltab[nodenames2, "stage2"] = rwr_res$vec[nodenames2]
        model = prediction(as.numeric(evaltab[testuni2, "stage2"]),
                                 evaltab[testuni2, "test"])
        auc = performance(model, "auc")
        perf_stage2 = performance(model, "tpr", "fpr")
        aucval = round(as.numeric(slot(auc, "y.values")), 3)
        row = c(iter, "stage2", aucval)
        restable = rbind(restable, row)

        perf_d@x.values[[iter]] <- unlist(perf_diff@x.values)
        perf_d@y.values[[iter]] <- unlist(perf_diff@y.values)
        perf_d@alpha.values[[iter]] <- unlist(perf_diff@alpha.values)

        perf_s2@x.values[[iter]] <- unlist(perf_stage2@x.values)
        perf_s2@y.values[[iter]] <- unlist(perf_stage2@y.values)
        perf_s2@alpha.values[[iter]] <- unlist(perf_stage2@alpha.values)
    }

    rownames(restable) <- 1:nrow(restable)

    #Stage 1 network
    query = blankvec
    midxs = match(queryIDs, nodenames)
    query[midxs] = as.numeric(query.genes[queryIDs, 2])
    rwr_res = RWR(boolSparceMat, transmat, restart, query, startvec, maxiters, thresh)
    qiter = rwr_res$iter
    evaltab = cbind(evaltabstart, 0, 0, as.numeric(rwr_res$vec[nodenames]))
    colnames(evaltab) = c("node", "type", "universe", "baseline", "train",
                          "test", "stage1")
    diff = as.numeric(evaltab[, "stage1"]) - as.numeric(evaltab[, "baseline"])
    evaltab = cbind(evaltab, diff)
    evaltab[queryIDs, "train"] = 1
    evaltab[queryIDs, "test"] = 1
    keep = rep(-1, nnodes)
    evaltab = cbind(evaltab, keep)
    ss = sort(as.numeric(evaltab[featnodes, "diff"]), decreasing = TRUE,
              index.return = TRUE)
    sortedfeats = featnodes[ss$ix]
    keepfeats = sortedfeats[1:nkeep]
    evaltab[featnodes, "keep"] = 0
    evaltab[keepfeats, "keep"] = 1
    newedges = edges[which(edges$type %in% setdiff(all_etypes,
                                                   prop_etypes)), ]
    keepidxs = which(edges[, "src"] %in% keepfeats |
                         edges[, "target"] %in% keepfeats)
    newedges = rbind(newedges, edges[keepidxs, ])
    tmpnames = unique(c(as.character(newedges[, 1]),
                        as.character(newedges[, 2])))
    queryIDs2 = intersect(queryIDs, tmpnames)
    typetable2 = NULL
    typetable2 = aggregate(weight ~ type, data = newedges,
                           FUN = sum)
    rownames(typetable2) = as.character(typetable2[,1])
    newedges$weight = newedges$weight/typetable2[as.character(newedges$type),
                                                 "weight"]
    transmat2 = NULL
    nodenames2 = NULL
    colsum2 = NULL
    boolSparceMat2 = (length(tmpnames)^2 < .Machine$integer.max) *
        (1 - forcelarge)
    if (boolSparceMat2) {
        n2Matrix = threeCol2MaxMat(as.character(newedges[,"src"]),
                                   as.character(newedges[,"target"]),
                                   as.numeric(newedges[,"weight"]))
        nodenames2 = rownames(n2Matrix)
        colsum2 = colSums(n2Matrix)
        transmat2 = t(t(n2Matrix)/colsum2)
        rm(n2Matrix)
    }
    else {
        ll2 = threeCol2listMat(as.character(newedges[,"src"]),
                               as.character(newedges[,"target"]),
                               as.numeric(newedges[,"weight"]))
        transmat2 = ll2
        nodenames2 = ll2$avals
        colsum2 = ll2$colsums
        rm(ll2)
    }
    nnodes2 = length(nodenames2)
    rm(newedges)
    blankvec2 = structure(rep(0, nnodes2), names = nodenames2)
    query2 = blankvec2
    midxs = match(queryIDs2, nodenames2)
    query2[midxs] = as.numeric(query.genes[queryIDs2, 2])
    startvec2 = blankvec2 + 1/nnodes2
    rwr_res = RWR(boolSparceMat2, transmat2, restart,
                  query2, startvec2, maxiters, thresh)
    q2iter = rwr_res$iter
    stage2 = rep(0, nnodes)
    evaltab = cbind(evaltab, stage2)
    evaltab[nodenames2, "stage2"] = rwr_res$vec[nodenames2]

    #Stagewise performance and ROC curves
    perftable <- data.frame(row.names = 1:nrow(restable),
                            fold = as.numeric(restable[,"iter"]),
                            stage = restable[,"stage"],
                            auc = as.numeric(restable[,"aucval"]),
                            stringsAsFactors = FALSE)
    perftable <- aggregate(auc ~ stage, data = perftable, mean)

    auc1 = perftable[perftable$stage == "diff", "auc"]
    auc2 = perftable[perftable$stage == "stage2", "auc"]
    plot(perf_d, avg = "threshold", col = "skyblue", lwd = 2)
    plot(perf_s2, avg = "threshold", col = "orange", lwd = 2, add = TRUE)
    legend("bottomright", legend = c(paste0("Stage 1 (AUC = ", auc1, ")"),
                                     paste0("Stage 2 (AUC = ", auc2, ")")),
           col = c("skyblue", "orange"), lwd = 2)

    #Feature nodes
    features <- rbind(evaltab[evaltab[,"keep"] == "1",])
    features <- features[order(as.numeric(features[,"stage2"]),
                               decreasing = TRUE),]
    featuretable <- data.frame(row.names = 1:nkeep,
                               node = features[1:nkeep,"node"],
                               type = features[1:nkeep,"type"],
                               prob = features[1:nkeep,"stage2"],
                               stringsAsFactors = FALSE)

    #Gene nodes
    genes <- rbind(evaltab[evaltab[,"type"] == "-1",])
    genes <- genes[order(as.numeric(genes[,"stage2"]), decreasing = TRUE),]
    genetable <- data.frame(row.names = 1:nrow(genes), node = genes[, "node"],
                            prob = genes[,"stage2"],
                            stringsAsFactors = FALSE)

    return(list(features = featuretable, genes = genetable))

}

