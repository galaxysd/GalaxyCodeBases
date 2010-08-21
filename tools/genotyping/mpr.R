# MPR-internal.R
`.loopMPR` <-
function (genoData, allele.matrix, numStep = 1, ...) 
{
    rawGenoSum <- NumRecomEvents(genoData = genoData)
    for (x in 1:(nrow(allele.matrix) - numStep + 1)) {
        ids <- x:(x + numStep - 1)
        new.genoData <- genoData
        new.genoData[ids, ] <- abs(new.genoData[ids, ] - 1)
        newGenoSum <- NumRecomEvents(genoData = new.genoData)
        if (newGenoSum < rawGenoSum) {
            allele.matrix[ids, ] <- allele.matrix[ids, c(2, 1)]
            genoData <- new.genoData
            rawGenoSum <- newGenoSum
        }
    }
    list(rawGenoSum, allele.matrix, genoData)
}
`.MPR_hetero_` <-
0.5
`.MPR_p0_` <-
0
`.MPR_p1_` <-
1
`._rice_phy2get_factor_` <-
244000

# NumRecomEvents.R
`NumRecomEvents` <-
function (baseData, allele.matrix, genoData = NULL) 
{
    if (is.null(genoData)) {
        genoData <- base2Geno(baseData, allele.matrix)
    }
    y <- !is.na(genoData)
    idsBorder <- cumsum(colSums(y))
    idsBorder <- idsBorder[idsBorder < sum(y)]
    sum(diff(as.numeric(genoData[y]))[-idsBorder] != 0, na.rm = TRUE)
}

#R A or G 
#Y C or T 
#S G or C 
#W A or T 
#K G or T 
#M A or C 
#B C or G or T 
#D A or G or T 
#H A or C or T 
#V A or C or G 
#N any base 
#        if (length(x[x == 'R'])) {
#            x[x == 'R']=NA;
#            x=c(x,'A','G');
#        }

# base2Allele.R
.IUPAC=c('R','Y','S','W','K','M','B','D','H','V','N');
.Decode=list(R=c('A','G'),Y=c('C','T'),S=c('G','C'),W=c('A','T'),K=c('G','T'),M=c('A','C'),B=c('C','G','T'),D=c('A','G','T'),H=c('A','C','T'),V=c('A','C','G'),N=NA);

`base2Allele` <-
function (baseData = NULL) 
{
    allele.matrix <- t(apply(baseData, 1, function(x) {
        x <- unique(x[!is.na(x)])

        for(i in 1:length(.IUPAC)) {
            if ( length(grep(.IUPAC[i],x)) ) {
                  x=x[-(grep(.IUPAC[i],x))]
                  x=c(x,.Decode[[i]])
            }
        }

        x <- unique(x[!is.na(x)])
        if (length(x) == 2) 
            return(x)
        else warning('SNP sit is not biallelic!')
        c(NA, NA)
    }))
    colnames(allele.matrix) <- c('P1', 'P2')
    na.exclude(allele.matrix)
}

# base2Geno.R
`base2Geno` <-
function (baseData = NULL, allele.matrix = NULL) 
{
    if (nrow(baseData) == ncol(allele.matrix)) 
        allele.matrix <- t(allele.matrix)
    if (nrow(baseData) != nrow(allele.matrix)) 
        stop('nrow(baseData)!=nrow(allele.matrix), allele.matrix error!!!')
    genoData <- baseData
    genoData[baseData == allele.matrix[, 1]] <- 0
    genoData[baseData == allele.matrix[, 2]] <- 1
    genoData[genoData != 1 & genoData != 0] <- NA
    genoData <- matrix(as.numeric(genoData), ncol = ncol(genoData))
    dimnames(genoData) <- dimnames(baseData)
    genoData
}

# calculateLODScore.R
`calculateLODScore` <-
function (geno.data, pheno.data, mapInfo) 
{
    geno.data.tmp <- geno.data > 0.5
    geno.data.tmp[geno.data.tmp == FALSE] <- NA
    ph1 <- matrix(t(matrix(pheno.data, nrow = length(pheno.data), 
        ncol = nrow(geno.data)))[geno.data.tmp], ncol = ncol(geno.data))
    geno.data.tmp <- geno.data < 0.5
    geno.data.tmp[geno.data.tmp == FALSE] <- NA
    ph0 <- matrix(t(matrix(pheno.data, nrow = length(pheno.data), 
        ncol = nrow(geno.data)))[geno.data.tmp], ncol = ncol(geno.data))
    g1 <- rowMeans(ph1, na.rm = T)
    n1 <- rowSums(!is.na(ph1))
    s1 <- rowVars(ph1, na.rm = T)
    g0 <- rowMeans(ph0, na.rm = T)
    n0 <- rowSums(!is.na(ph0))
    s0 <- rowVars(ph0, na.rm = T)
    n <- rowSums(!is.na(geno.data))
    ss <- (s1 * (n1 - 1) + s0 * (n0 - 1))/(n1 + n0 - 2)
    t <- (g1 - g0)/sqrt(ss/n1 + ss/n0)
    res <- data.frame(mapInfo, '-log10 P' = -log10(pt(abs(t), 
        df = pmax(1, n1 + n0 - 2), lower.tail = F, log.p = FALSE) * 
        2), t = t, d = g1 - g0, n0 = n0, n1 = n1, check.names = FALSE)
    class(res) <- c('scanone', 'data.frame')
    res
}

# correctFUNHMM.R
`correctFUNHMM` <-
function (x, base.position = as.numeric(sub('[^0-9]*([0-9]*)[^0-9]*', 
    '\\1', names(x))), hmmFUN = hmm.vitFUN.rils, geno.probability = c(0.495, 
    0.495, 0.01), transitionFUN = phy2get.haldane.rils, emissionFUN = makeEmissionFUN(errorRate = 0.01), 
    ...) 
{
    x.cr <- hmmFUN(geno = x + 1, position = base.position, geno.probability = geno.probability, 
        transitionFUN = transitionFUN, emissionFUN = emissionFUN, 
        ...) - 1
    x.cr[x.cr == 2] <- .MPR_hetero_
    x.cr
}

# correctGeno.R
`correctGeno` <-
function (geno.data, base.position = as.numeric(sub('[^0-9]*([0-9]*)[^0-9]*', 
    '\\1', rownames(geno.data))), correct.FUN = correctFUNHMM, 
    minInterval = 1, verbose = TRUE, ...) 
{
    i <- 0
    geno.data.cr <- apply(geno.data, 2, function(x, ...) {
        if (verbose) 
            cat('\r', i <<- i + 1)
        x.nna <- which(!is.na(x))
        ids <- sort(unique(x.nna[c(1, which(diff(x[x.nna]) != 
            0 | diff(base.position[x.nna]) >= minInterval) + 
            1)]))
        x.cr <- correct.FUN(x[ids], base.position = base.position[ids], 
            ...)
        x[x.nna] <- NA
        x[ids] <- x.cr
        x
    }, ...)
    if (verbose) 
        cat('\tDone.\n')
    geno.data.cr
}

# findBlockAndFilter.R
`findBlockAndFilter` <-
function (x, base.position = 1:length(x), size = sum(blocks[, 
    'size'], na.rm = TRUE)/100, num = sum(blocks[, 'num'], na.rm = TRUE)/100, 
    extendNA = TRUE, fillSmallNA = FALSE) 
{
    x <- as.numeric(x)
    t <- length(x)
    block.num <- 1
    block_start <- base.position[1]
    block_end <- base.position[t]
    blockNumItem <- t
    blockType <- x[1]
    x.val <- x
    x.val[is.na(x)] <- max(x, na.rm = TRUE) + 1
    x.diff <- sort(unique(c(1, which(x.val[-1] != x.val[-t]) + 
        1)))
    x.diff <- x.diff[x.diff <= t]
    block.num <- length(x.diff)
    if (block.num <= 1) 
        return(blocks <- cbind(start = block_start, end = block_end, 
            size = block_end - block_start + 1, num = blockNumItem, 
            type = blockType))
    block_start <- base.position[x.diff]
    block_end <- base.position[c(x.diff[-1] - 1, t)]
    blockNumItem <- c(diff(x.diff), t - x.diff[block.num] + 1)
    blockType <- x[x.diff]
    blocks <- cbind(start = block_start, end = block_end, size = block_end - 
        block_start + 1, num = blockNumItem, type = blockType)
    if (extendNA) {
        block.ids <- which(is.na(blocks[, 'type']) & blocks[, 
            'size'] > -1)
        block.ids.next <- block.ids[block.ids < nrow(blocks)]
        block.ids.pre <- block.ids[block.ids > 1]
        if (length(block.ids.next) > 0) 
            blocks[block.ids.next, 'end'] <- blocks[block.ids.next + 
                1, 'start'] - 1
        if (length(block.ids.pre) > 0) 
            blocks[block.ids.pre, 'start'] <- blocks[block.ids.pre - 
                1, 'end'] + 1
        blocks[, 'size'] <- blocks[, 'end'] - blocks[, 'start'] + 
            1
    }
    blocks.margin <- block_start[-1] - block_end[-length(block_end)] - 
        1
    blocks.extentSize <- as.numeric(blocks[, 'size']) + c(blocks.margin, 
        0) + c(0, blocks.margin)
    filter.block <- (blocks.extentSize < size & as.numeric(blocks[, 
        'num']) < num)
    if (sum(filter.block, na.rm = T) > 0) {
        blocks[filter.block, 'type'] <- NA
    }
    blocks <- mergeBlocks(blocks)
    if (nrow(blocks) == 1) 
        return(blocks)
    if (fillSmallNA) {
        blocks.margin <- blocks[-1, 'start'] - blocks[-nrow(blocks), 
            'end'] - 1
        blocks.extentSize <- as.numeric(blocks[, 'size']) + c(blocks.margin, 
            0) + c(0, blocks.margin)
        filter.block <- blocks.extentSize < size
        block.na.ids <- which(is.na(blocks[, 'type']) & filter.block)
        block.na.num <- length(block.na.ids)
        if (block.na.num > 0) {
            blocks[block.na.ids, 'type'] <- sapply(block.na.ids, 
                function(i) {
                  if (i == 1) 
                    return(blocks[2, 'type'])
                  if (i == nrow(blocks)) 
                    return(blocks[nrow(blocks) - 1, 'type'])
                  if (blocks[i - 1, 'type'] == blocks[i + 1, 
                    'type']) 
                    return(blocks[i - 1, 'type'])
                  blocks[i, 'type']
                })
        }
    }
    blocks <- mergeBlocks(blocks)
    if (nrow(blocks) == 1) 
        return(blocks)
    block.ids <- which(is.na(blocks[, 'type']))
    block.ids.next <- block.ids[block.ids < nrow(blocks)]
    block.ids.pre <- block.ids[block.ids > 1]
    if (length(block.ids.next) > 0) 
        blocks[block.ids.next, 'end'] <- blocks[block.ids.next + 
            1, 'start'] - 1
    if (length(block.ids.pre) > 0) 
        blocks[block.ids.pre, 'start'] <- blocks[block.ids.pre - 
            1, 'end'] + 1
    tmp.blocks <- blocks
    for (i in 2:nrow(blocks)) if (blocks[i, 1] - blocks[i - 1, 
        2] > 1) 
        tmp.blocks <- rbind(tmp.blocks, c(blocks[i - 1, 'end'] + 
            1, blocks[i, 'start'] - 1, blocks[i, 'start'] - blocks[i - 
            1, 'end'] - 1, 0, NA))
    blocks <- tmp.blocks[order(tmp.blocks[, 'end']), ]
    mergeBlocks(blocks)
}

# geno2Cross.R
`geno2Cross` <-
function (geno.data, pheno.data, cm_unit = 250000) 
{
    myGeno.data <- geno.data
    geno.value <- sort(unique(na.omit(c(myGeno.data))))
    if (identical(as.numeric(geno.value), c(0, 1))) {
        myGeno.data <- myGeno.data + 1
        geno.value <- sort(unique(na.omit(c(myGeno.data))))
    }
    if (!identical(as.numeric(geno.value), c(1, 2))) 
        stop('the value of geno.data should be 1, 2, or NA.')
    ids.RILs <- match(colnames(myGeno.data), rownames(pheno.data))
    myGeno.site <- data.frame(Chr = as.numeric(substr(rownames(myGeno.data), 
        1, 2)), Position = as.numeric(substr(rownames(myGeno.data), 
        3, 10)))
    myCrossData <- list()
    myCrossData$geno <- lapply(split(1:nrow(myGeno.data), myGeno.site$Chr), 
        function(ids) {
            myMarkers <- myGeno.site$Position[ids]/cm_unit
            names(myMarkers) <- rownames(myGeno.data)[ids]
            myMap <- list(data = t(myGeno.data[ids, ]), map = myMarkers)
            class(myMap) <- 'A'
            myMap
        })
    myCrossData$pheno <- as.data.frame(pheno.data[ids.RILs, ])
    class(myCrossData) <- c('riself', 'cross')
    attr(myCrossData, 'alleles') <- c('A', 'B')
    myCrossData
}

# genoToBin.R
`genoToBin` <-
function (genoData, base.position = 1:nrow(genoData), corrected = FALSE, 
    correct.FUN = correctFUNHMM, size = 250000, num = 5, fillSmallNA = TRUE, 
    minBinsize = 0, seqERR = 0.01, heterozygote = FALSE, ...) 
{
    SNPbyChr <- split(1:nrow(genoData), 1)
    i.chr <- 0
    res <- lapply(SNPbyChr, function(ids) {
        i.line <- 0
        i.chr <<- i.chr + 1
        blocks <- apply(genoData[ids, ], 2, function(x) {
            cat('chr: ', i.chr, '\tline: ', i.line <<- i.line + 
                1, '\t', colnames(genoData)[i.line], '\r')
            x.nna <- !is.na(x)
            if (corrected == FALSE) 
                x.correct <- correct.FUN(x[x.nna], base.position[ids][x.nna], 
                  ...)
            else x.correct <- x[x.nna]
            blocks.mat <- findBlockAndFilter(x.correct, base.position = base.position[ids][x.nna], 
                size = size, num = num, fillSmallNA = fillSmallNA)
            blocks.mat[nrow(blocks.mat), 'end'] <- max(base.position[ids])
            blocks.mat[nrow(blocks.mat), 'size'] <- blocks.mat[nrow(blocks.mat), 
                'end'] - blocks.mat[nrow(blocks.mat), 'start'] + 
                1
            if (heterozygote == FALSE) 
                blocks.mat[blocks.mat[, 'type'] == .MPR_hetero_, 
                  'type'] <- NA
            blocks.mat <- mergeBlocks(blocks.mat)
            t(blocks.mat)
        })
        blocks.mat <- matrix(unlist(blocks, recursive = TRUE, 
            use.names = FALSE), ncol = 5, byrow = TRUE)
        bin.border <- sort(unique(as.numeric(blocks.mat[, 2])))
        cat('\rchr: ', i.chr, '\tTotal', length(bin.border), 
            'borders. ')
        bin.border <- sort(unique(bin.border[c(1, which(diff(bin.border) >= 
            minBinsize) + 1)]))
        cat(length(bin.border), 'borders after filtering out bins less than', 
            round(minBinsize/1000, 1), 'kb.\n')
        geno.bin <- sapply(blocks, function(blocks.line) {
            blocks.end <- rbind(blocks.line)[2, ]
            ids <- match(bin.border, blocks.end)
            ids[is.na(ids)] <- findInterval(bin.border[is.na(ids)], 
                blocks.end, rightmost.closed = FALSE) + 1
            filter <- ids > ncol(blocks.line)
            ids[filter] <- ncol(blocks.line)
            rbind(blocks.line)[5, ids]
        })
        rownames(geno.bin) <- sprintf('%010d', bin.border)
        list(block = blocks, bin = geno.bin, border = bin.border)
    })
    geno.bin <- NULL
    for (i in 1:length(res)) geno.bin <- rbind(geno.bin, res[[i]][[2]])
    cat('Done.\n')
    list(block = res, bin = geno.bin, border = as.numeric(rownames(geno.bin)))
}

# genotypeCallsBayes.R
`genotypeCallsBayes` <-
function (ALLELE.num, errorRate = 5e-04, eps = 1e-10, maxIterate = 100, 
    verbose = FALSE) 
{
    P0.n <- 0.5
    P0.p1 <- P0.p2 <- (1 - P0.n)/2
    ALLELE.prob <- NULL
    P.n <- 1/3
    P.p1 <- P.p2 <- (1 - P.n)/2
    numIterate <- 1
    E <- errorRate
    while (abs(P0.n - P.n) > eps && numIterate <= maxIterate) {
        P0.p1 <- P.p1
        P0.p2 <- P.p2
        P0.n <- P.n
        n <- rowSums(ALLELE.num)
        k <- ALLELE.num[, 1]
        ALLELE.prob <- cbind((1 - E)^k * E^(n - k), 0.5^n, (1 - 
            E)^(n - k) * E^k) %*% diag(c(P0.p1, P0.n, P0.p2))
        ALLELE.type <- rep(NA, nrow(ALLELE.prob))
        ALLELE.maxprob <- rowMax(ALLELE.prob)
        ALLELE.type[ALLELE.prob[, 1] == ALLELE.maxprob] <- 1
        ALLELE.type[ALLELE.prob[, 3] == ALLELE.maxprob] <- 3
        ALLELE.type[ALLELE.prob[, 2] == ALLELE.maxprob | ALLELE.prob[, 
            1] == ALLELE.prob[, 3]] <- 2
        P.n <- mean(ALLELE.type == 2)
        P.p1 <- P.p2 <- (1 - P.n)/2
        if (verbose) 
            cat(P.p1, P.n, P.p2, table(filter.markers <- ALLELE.type != 
                2), '\n', sep = '\t')
        numIterate <- numIterate + 1
    }
    ALLELE.prob <- ALLELE.prob/rowSums(ALLELE.prob, na.rm = TRUE)
    list(prop = c(P.p1, P.n, P.p2), prob = ALLELE.prob, type = ALLELE.type)
}

# globalMPRByMarkers.R
`globalMPRByMarkers` <-
function (baseData, markers = NULL, alleleA = NULL, numTry = 3, 
    numBaseStep = 50, numBaseCandidateStep = numBaseStep * 2, 
    numKnownStep = pmax(numBaseStep/5, 10), numKnownCandidateStep = numKnownStep * 
        1.5, useMedianToFindKnown = TRUE, maxIterate = 150, maxNStep = 3, 
    scoreMin = 0.8, verbose = FALSE, strSTART = '\r', strEND = '', 
    ...) 
{
    if (is.null(alleleA)) 
        alleleA <- markers
    ALLELE.mat <- matrix(NA, nrow = nrow(baseData), ncol = 2)
    rownames(ALLELE.mat) <- rownames(baseData)
    alleleA.base <- alleleA[match(rownames(ALLELE.mat), names(alleleA))]
    alleleA.ids <- which(!is.na(alleleA.base))
    if (numKnownCandidateStep > length(alleleA.ids)) 
        numKnownCandidateStep <- length(alleleA.ids)
    if (numKnownStep > numKnownCandidateStep) 
        numKnownStep <- numKnownCandidateStep
    j <- 0
    ids.RILrows <- which(rowSums(!is.na(cbind(baseData))) > 0)
    rowN <- length(ids.RILrows)
    ids.times <- rep(0, rowN)
    ids.ok <- rep(0, rowN)
    ids.candidate <- na.omit(which(ids.times < numTry & ids.ok == 
        0)[1:numBaseCandidateStep])
    n <- length(ids.candidate)
    while (n > 1) {
        if (length(ids.candidate) > numBaseStep) {
            filter.dis <- ids.candidate < (median(ids.candidate) - 
                numBaseCandidateStep)
            ids.dis <- ids.candidate[filter.dis]
            if (length(ids.dis) < numBaseStep) 
                ids.candidate <- c(ids.dis, sample(ids.candidate[!filter.dis], 
                  numBaseStep - length(ids.dis)))
            else ids.candidate <- ids.dis
        }
        ids.times[ids.candidate] <- ids.times[ids.candidate] + 
            1
        ids <- ids.RILrows[ids.candidate]
        ids.point <- ifelse(useMedianToFindKnown == TRUE, median(ids), 
            ids[1])
        is.known <- !is.na(alleleA.base[ids])
        if (sum(is.known) < numKnownStep) {
            ids.known <- na.omit(alleleA.ids[order(abs(alleleA.ids - 
                ids.point))[sample(numKnownCandidateStep, numKnownStep)]])
            if (length(ids.known) > (numKnownStep - sum(is.known))) 
                ids.know <- ids.known[1:(numKnownStep - sum(is.known))]
            ids <- unique(c(ids, ids.known))
        }
        ids <- sort(ids)
        is.known <- !is.na(alleleA.base[ids])
        iResult <- localMPR(baseData[ids, ], maxIterate = maxIterate, 
            maxNStep = maxNStep, returnNumIterate = TRUE, verbose = 0)
        allele.matrix <- iResult[['allele']]
        a <- allele.matrix[is.known, ]
        b <- alleleA.base[ids[is.known]]
        a1 <- colSums(a == b, na.rm = T)/length(b)
        if (sum(a1 >= scoreMin) > 0) {
            if (a1[1] < a1[2]) 
                allele.matrix <- allele.matrix[, c(2, 1)]
            j <- j + 1
            a <- ALLELE.mat[ids, ]
            ids.na <- rowSums(is.na(a)) > 0
            if (sum(ids.na) > 0) 
                ALLELE.mat[ids[ids.na], ] <- allele.matrix[ids.na, 
                  ]
            ids.ok[ids.candidate] <- 1
        }
        ids.all <- which(ids.times < numTry & ids.ok == 0)
        ids.candidate <- na.omit(ids.all[1:numBaseCandidateStep])
        n <- length(ids.candidate)
        if (verbose) 
            cat(strSTART, length(ids.all), j, strEND, sep = '\t')
    }
    invisible(ALLELE.mat)
}

# globalMPRByMarkersPlus.R
`globalMPRByMarkersPlus` <-
function (baseData, markers = NULL, alleleA = NULL, numGroup = 1, 
    groupSort = FALSE, numPerm = 1, numTry = 3, numBaseStep = 50, 
    numBaseCandidateStep = numBaseStep * 2, numKnownStep = pmax(numBaseStep/5, 
        10), numKnownCandidateStep = numKnownStep * 2, useMedianToFindKnown = TRUE, 
    maxIterate = 150, maxNStep = 3, scoreMin = 0.8, saveMidData = FALSE, 
    verbose = FALSE, strSTART = '\r', strEND = '', ...) 
{
    if (is.null(alleleA)) 
        alleleA <- markers
    ALLELE.mat <- base2Allele(baseData)
    ALLELE.num <- matrix(0, nrow = nrow(ALLELE.mat), ncol = ncol(ALLELE.mat))
    dimnames(ALLELE.num) <- dimnames(ALLELE.mat)
    alleleA.base <- alleleA[match(rownames(ALLELE.mat), names(alleleA))]
    alleleA.ids <- which(!is.na(alleleA.base))
    if (numKnownCandidateStep > length(alleleA.ids)) 
        numKnownCandidateStep <- length(alleleA.ids)
    if (numKnownStep > numKnownCandidateStep) 
        numKnownStep <- numKnownCandidateStep
    j <- 0
    midData <- NULL
    for (perm in 1:numPerm) {
        if (numGroup >= 2) {
            if (groupSort) 
                groupRIL <- split(order(colSums(!is.na(baseData)), decreasing = TRUE), cut(1:ncol(baseData), numGroup))
            else groupRIL <- split(sample(ncol(baseData)), cut(1:ncol(baseData), numGroup))
            names(groupRIL) <- 1:numGroup
        }
        else {
            groupRIL <- list('1' = 1:ncol(baseData))
        }
        groupRIL <- groupRIL[sapply(groupRIL, length) > 0]
        for (ids.cols in groupRIL) {
            ids.RILrows <- which(rowSums(!is.na(cbind(baseData[, ids.cols]))) > 0)
            rowN <- length(ids.RILrows)
            ids.times <- rep(0, rowN)
            ids.ok <- rep(0, rowN)
            ids.candidate <- na.omit(which(ids.times < numTry & 
                ids.ok == 0)[1:numBaseCandidateStep])
            n <- length(ids.candidate)
            while (n > 1) {
                if (length(ids.candidate) > numBaseStep) {
                  filter.dis <- ids.candidate < (median(ids.candidate) - numBaseCandidateStep)
                  ids.dis <- ids.candidate[filter.dis]
                  if (length(ids.dis) < numBaseStep) 
                    ids.candidate <- c(ids.dis, sample(ids.candidate[!filter.dis], numBaseStep - length(ids.dis)))
                  else ids.candidate <- ids.dis
                }
                ids.times[ids.candidate] <- ids.times[ids.candidate] + 1
                ids <- ids.RILrows[ids.candidate]
                ids.point <- ifelse(useMedianToFindKnown == TRUE, median(ids), ids[1])
                is.known <- !is.na(alleleA.base[ids])
                if (sum(is.known) < numKnownStep) {
                  ids.known <- na.omit(alleleA.ids[order(abs(alleleA.ids - ids.point))[sample(numKnownCandidateStep, numKnownStep)]])
                  if (length(ids.known) > (numKnownStep - sum(is.known))) 
                    ids.know <- ids.known[1:(numKnownStep - sum(is.known))]
                  ids <- unique(c(ids, ids.known))
                }
                ids <- sort(ids)
                is.known <- !is.na(alleleA.base[ids])
                iResult <- localMPR(baseData[ids, ], allele.matrix = ALLELE.mat[ids, 
                  ], maxIterate = maxIterate, maxNStep = maxNStep, 
                  returnNumIterate = TRUE, verbose = 0)
                allele.matrix <- iResult[['allele']]
                a <- allele.matrix[is.known, ]
                b <- alleleA.base[ids[is.known]]
                a1 <- colSums(a == b, na.rm = T)/length(b)
                if (sum(a1 >= scoreMin) > 0) {
                  if (a1[1] < a1[2]) 
                    allele.matrix <- allele.matrix[, c(2, 1)]
                  j <- j + 1
                  a <- ALLELE.mat[ids, ]
                  a1 <- rowSums(a == allele.matrix) == 2
                  ALLELE.num[ids, ] <- ALLELE.num[ids, ] + cbind(a1, 
                    !a1)
                  ids.ok[ids.candidate] <- 1
                }
                ids.all <- which(ids.times < numTry & ids.ok == 
                  0)
                ids.candidate <- na.omit(ids.all[1:numBaseCandidateStep])
                n <- length(ids.candidate)
                if (verbose) 
                  cat(strSTART, 'perm', perm, 'unprocessed', 
                    length(ids.all), j, strEND, sep = '\t')
            }
        }
        filter.exchange <- ALLELE.num[, 1] < ALLELE.num[, 2]
        tmp <- ALLELE.mat
        tmp[filter.exchange, ] <- tmp[filter.exchange, c(2, 1)]
        ALLELE.mat <- tmp
        tmp <- ALLELE.num
        tmp[filter.exchange, ] <- tmp[filter.exchange, c(2, 1)]
        ALLELE.num <- tmp
        if (saveMidData) 
            midData[[perm]] <- list(allele = ALLELE.mat, call = ALLELE.num)
    }
    invisible(list(allele = ALLELE.mat, call = ALLELE.num, midData = midData))
}

# globalMPRRefine.R
`globalMPRRefine` <-
function (baseData, markers = NULL, alleleA = NULL, numGroup = ncol(baseData), 
    groupSort = FALSE, numPerm = 10, numTry = 3, numBaseStep = 50, 
    numBaseCandidateStep = numBaseStep * 2, numKnownStep = numBaseStep/2, 
    numKnownCandidateStep = numBaseStep * 2, useMedianToFindKnown = TRUE, 
    maxIterate = 150, maxNStep = 3, scoreMin = 0.8, useOnlyKnownToType = FALSE, 
    useBayes = FALSE, errorRate = 5e-04, saveMidData = FALSE, 
    verbose = FALSE, strSTART = '\r', strEND = '', ...) 
{
    if (is.null(alleleA)) 
        alleleA <- markers
    ALLELE.mat <- base2Allele(baseData)
    ALLELE.num <- matrix(0, nrow = nrow(ALLELE.mat), ncol = ncol(ALLELE.mat))
    dimnames(ALLELE.num) <- dimnames(ALLELE.mat)
    alleleA.base <- alleleA[match(rownames(ALLELE.mat), names(alleleA))]
    alleleA.ids <- which(!is.na(alleleA.base))
    ALLELE.num[ALLELE.mat == alleleA.base] <- 1
    tmpALLELE.num <- ALLELE.num
    if (numKnownCandidateStep > length(alleleA.ids)) 
        numKnownCandidateStep <- length(alleleA.ids)
    if (numKnownStep > numKnownCandidateStep) 
        numKnownStep <- numKnownCandidateStep
    j <- 0
    midData <- NULL
    for (perm in 1:numPerm) {
        if (numGroup >= 2) {
            if (groupSort) 
                groupRIL <- split(order(colSums(!is.na(baseData)), 
                  decreasing = TRUE), cut(1:ncol(baseData), numGroup))
            else groupRIL <- split(sample(ncol(baseData)), cut(1:ncol(baseData), 
                numGroup))
            names(groupRIL) <- 1:numGroup
        }
        else {
            groupRIL <- list('1' = 1:ncol(baseData))
        }
        groupRIL <- groupRIL[sapply(groupRIL, length) > 0]
        for (g in 1:length(groupRIL)) {
            ids.cols <- groupRIL[[g]]
            ids.RILrows <- which(rowSums(!is.na(cbind(baseData[, 
                ids.cols]))) > 0)
            rowN <- length(ids.RILrows)
            ids.times <- rep(0, rowN)
            ids.ok <- rep(0, rowN)
            ids.candidate <- na.omit(which(ids.times < numTry & 
                ids.ok == 0)[1:numBaseCandidateStep])
            n <- length(ids.candidate)
            while (n > 1) {
                if (length(ids.candidate) > numBaseStep) {
                  filter.dis <- ids.candidate < (median(ids.candidate) - 
                    numBaseCandidateStep)
                  ids.dis <- ids.candidate[filter.dis]
                  if (length(ids.dis) < numBaseStep) 
                    ids.candidate <- c(ids.dis, sample(ids.candidate[!filter.dis], 
                      numBaseStep - length(ids.dis)))
                  else ids.candidate <- ids.dis
                }
                ids.times[ids.candidate] <- ids.times[ids.candidate] + 
                  1
                ids <- ids.RILrows[ids.candidate]
                ids.point <- ifelse(useMedianToFindKnown == TRUE, 
                  median(ids), ids[1])
                is.known <- !is.na(alleleA.base[ids])
                if (sum(is.known) < numKnownStep) {
                  ids.known <- na.omit(alleleA.ids[order(abs(alleleA.ids - 
                    ids.point))[sample(numKnownCandidateStep, 
                    numKnownStep)]])
                  if (length(ids.known) > (numKnownStep - sum(is.known))) 
                    ids.know <- ids.known[1:(numKnownStep - sum(is.known))]
                  ids <- unique(c(ids, ids.known))
                }
                ids <- sort(as.numeric(ids))
                is.known <- !is.na(alleleA.base[ids])
                iResult <- localMPR(baseData[ids, ], allele.matrix = ALLELE.mat[ids, 
                  ], maxIterate = maxIterate, maxNStep = maxNStep, 
                  returnNumIterate = TRUE, verbose = 0)
                allele.matrix <- iResult[['allele']]
                if (useOnlyKnownToType) {
                  a <- allele.matrix[is.known, ]
                  b <- alleleA.base[ids[is.known]]
                  a1 <- colSums(a == b, na.rm = T)/length(b)
                }
                else {
                  a <- ALLELE.mat[ids, ]
                  b <- ALLELE.num[ids, ]
                  a.x <- a
                  b.f <- b[, 1] < b[, 2]
                  a.x[b.f, ] <- a.x[b.f, c(2, 1)]
                  a1 <- colMeans(allele.matrix == a.x[, 1], na.rm = T)
                }
                if (sum(a1 >= scoreMin) > 0) {
                  if (a1[1] < a1[2]) 
                    allele.matrix <- allele.matrix[, c(2, 1)]
                  j <- j + 1
                  a <- ALLELE.mat[ids, ]
                  a1 <- rowSums(a == allele.matrix) == 2
                  ALLELE.num[ids, ] <- ALLELE.num[ids, ] + cbind(a1, 
                    !a1)
                  ids.ok[ids.candidate] <- 1
                }
                ids.all <- which(ids.times < numTry & ids.ok == 
                  0)
                ids.candidate <- na.omit(ids.all[1:numBaseCandidateStep])
                n <- length(ids.candidate)
                if (verbose) 
                  cat(strSTART, 'perm', perm, 'group', g, 'unprocessed', 
                    length(ids.all), 'MPR', j, strEND, sep = '\t')
            }
            if (useBayes) {
                geno.res <- genotypeCallsBayes(ALLELE.num, errorRate = errorRate, 
                  eps = 1e-10, maxIterate = 100, verbose = FALSE)$type
                tmp <- ALLELE.mat
                tmp[geno.res == 2, ] <- NA
                tmp[geno.res == 3, ] <- tmp[geno.res == 3, c(2, 
                  1)]
            }
            else {
                tmp <- ALLELE.mat
                filter <- ALLELE.num[, 1] < ALLELE.num[, 2]
                tmp[filter, ] <- tmp[filter, c(2, 1)]
            }
            alleleA.base <- tmp[, 1]
            alleleA.ids <- which(!is.na(alleleA.base))
        }
        if (saveMidData) 
            midData[[perm]] <- list(allele = ALLELE.mat, call = ALLELE.num - 
                tmpALLELE.num, base = alleleA.base)
    }
    invisible(list(allele = ALLELE.mat, call = ALLELE.num - tmpALLELE.num, 
        base = alleleA.base, midData = midData))
}

# hmm.vitFUN.rils.R
`hmm.vitFUN.rils` <-
function (geno, position, geno.probability, transitionFUN = phy2get.haldane.rils, 
    emissionFUN = makeEmissionFUN(errorRate = 0.01), ...) 
{
    n.obs <- length(geno)
    n.state <- length(geno.probability)
    psi <- delta <- matrix(0, nrow = n.state, ncol = n.obs)
    n.con <- geno.cr <- numeric(n.obs)
    geno.dis <- abs(diff(position))
    n.con[1] <- 1
    g <- geno[1]
    for (i in 2:n.obs) {
        n.con[i] <- ifelse(geno[i] == g, n.con[i - 1] + 1, 1)
        g <- geno[i]
    }
    for (i in 1:n.state) delta[i, 1] <- log(geno.probability[i]) + 
        emissionFUN(i, geno[1], n.con[1])
    preProb <- numeric(n.state)
    for (t in 2:n.obs) {
        for (j in 1:n.state) {
            for (i in 1:n.state) preProb[i] <- delta[i, t - 1] + 
                log(transitionFUN(i, j, geno.dis[t - 1]))
            psi[j, t] <- which.max(preProb)
            delta[j, t] <- max(preProb) + emissionFUN(j, geno[t], 
                n.con[t])
        }
    }
    geno.cr[n.obs] <- which.max(delta[, n.obs])
    for (t in seq(n.obs - 1, 1, by = -1)) geno.cr[t] <- psi[geno.cr[t + 
        1], t + 1]
    geno.cr
}

# localMPR.R
`localMPR` <-
function (baseData, allele.matrix = NULL, maxIterate = 50, returnNumIterate = FALSE, 
    maxNStep = 5, verbose = FALSE, strEND = '\n', ...) 
{
    if (is.null(allele.matrix)) 
        allele.matrix <- base2Allele(baseData)
    if (nrow(baseData) == ncol(allele.matrix)) 
        allele.matrix <- t(allele.matrix)
    if (nrow(baseData) != nrow(allele.matrix)) 
        stop('nrow(baseData)!=nrow(allele.matrix), allele.matrix error!!!')
    newALLELE <- ALLELE <- allele.matrix
    genoData <- base2Geno(baseData, allele.matrix)
    newRSum <- oriRSum <- length(baseData)
    numStepSize <- 1
    numIterate <- 0
    while (numStepSize <= maxNStep) {
        numIterate <- numIterate + 1
        if (verbose) 
            cat('\r', numIterate, '\t', newRSum, '\t', numStepSize, 
                '\n')
        oriRSum <- newRSum
        loopRes <- .loopMPR(genoData, allele.matrix = ALLELE, 
            numStep = numStepSize)
        newRSum <- loopRes[[1]]
        if (oriRSum > newRSum) {
            ALLELE <- loopRes[[2]]
            genoData <- loopRes[[3]]
            numStepSize <- 1
        }
        else {
            numStepSize <- numStepSize + 1
        }
    }
    if (verbose) 
        cat('\tDone.', strEND)
    rownames(ALLELE) <- rownames(baseData)
    if (returnNumIterate) 
        list(allele = ALLELE, num_iterate = numIterate, numR = newRSum)
    else ALLELE
}

# makeEmissionFUN.R
`makeEmissionFUN` <-
function (errorRate = 0.01) 
{
    E <- log(errorRate)
    E2 <- log(1 - errorRate)
    E3 <- log(0.5)
    function(h, x, n) {
        if (h != 3) 
            return(ifelse(h == x, E2, E))
        else return(n * E3)
    }
}

# mergeBinMap.R
`mergeBinMap` <-
function (geno.fill) 
{
    pre <- NULL
    pre.nna <- NULL
    IsUniq <- rep(1, nrow(geno.fill))
    for (j in 1:nrow(geno.fill)) {
        cur.nna <- which(!is.na(geno.fill[j, ]))
        cur <- geno.fill[j, cur.nna]
        if (identical(pre.nna, cur.nna) && identical(pre, cur)) 
            IsUniq[j] <- 0
        else {
            pre <- cur
            pre.nna <- cur.nna
        }
    }
    geno.fill[unique(c(which(IsUniq == 1)[-1] - 1, nrow(geno.fill))), 
        ]
}

# mergeBlocks.R
`mergeBlocks` <-
function (blocks) 
{
    if (nrow(blocks) < 2) 
        return(blocks)
    type.val <- blocks[, 'type']
    t <- length(type.val)
    type.val[is.na(blocks[, 'type'])] <- as.numeric(max(blocks[, 
        'type'], na.rm = TRUE)) + 1
    block.start.ids <- sort(unique(c(1, which(type.val[-1] != 
        type.val[-t]) + 1)))
    block.start.ids <- block.start.ids[block.start.ids <= nrow(blocks)]
    block.block.num <- length(block.start.ids)
    blocks[-block.start.ids, 'start'] <- NA
    blocks[block.start.ids, c('end', 'num')] <- t(sapply(1:block.block.num, 
        function(i) {
            ids <- block.start.ids[i]:ifelse(i >= block.block.num, 
                nrow(blocks), block.start.ids[i + 1] - 1)
            c(blocks[ids[length(ids)], 'end'], sum(blocks[ids, 
                'num'], na.rm = TRUE))
        }))
    blocks <- rbind(blocks[!is.na(blocks[, 'start']), ])
    blocks[, 'size'] <- blocks[, 'end'] - blocks[, 'start'] + 
        1
    blocks
}

# phy2get.haldane.rils.R
`phy2get.haldane.rils` <-
function (a, b, l) 
{
    d <- l/(._rice_phy2get_factor_ * 100)
    p <- (1 - exp(-2 * d))
    p <- p/(1 + p)
    ifelse(a == b, 1 - p, p)
}

# Append for de-dependence rowQ from package:Biobase
`rowMax` <-
function (matrix)
{
	rowmax <-apply(matrix,1,function(x){
		return(max(x))})
	rowmax
}
`rowMin` <-
function (matrix)
{
	rowmin <-apply(matrix,1,function(x){
		return(min(x))})
	rowmin
}
