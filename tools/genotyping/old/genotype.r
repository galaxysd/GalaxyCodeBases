#!/ifs1/GAG/population/huxuesong/chroot/bin/r

self='./genotype.r';

if (is.null(argv) | length(argv)<2) {
  cat("Error: No options. Run",self,"<input> <output_prefix>.{bin,png,block}\n")
  q(status=1)
}


geno.data=read.table(argv[1],T,row.names=1)

`.MPR_hetero_` <-0.5
`.MPR_p0_` <-0
`.MPR_p1_` <-1
`._rice_phy2get_factor_` <-244000

`phy2get.haldane.rils` <-
function (a, b, l) 
{
    d <- l/(._rice_phy2get_factor_ * 100)
    p <- (1 - exp(-2 * d))
    p <- p/(1 + p)
    ifelse(a == b, 1 - p, p)
}

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

`correctFUNHMM` <-
function (x, base.position = as.numeric(sub("[^0-9]*([0-9]*)[^0-9]*", 
    "\\1", names(x))), hmmFUN = hmm.vitFUN.rils, geno.probability = c(0.495, 
    0.495, 0.01), transitionFUN = phy2get.haldane.rils, emissionFUN = makeEmissionFUN(errorRate = 0.01), 
    ...) 
{
    x.cr <- hmmFUN(geno = x + 1, position = base.position, geno.probability = geno.probability, 
        transitionFUN = transitionFUN, emissionFUN = emissionFUN, 
        ...) - 1
    x.cr[x.cr == 2] <- .MPR_hetero_
    x.cr
}

`mergeBlocks` <-
function (blocks)
{
    if (nrow(blocks) < 2)
        return(blocks)
    type.val <- blocks[, "type"]
    t <- length(type.val)
    type.val[is.na(blocks[, "type"])] <- as.numeric(max(blocks[,
        "type"], na.rm = TRUE)) + 1
    block.start.ids <- sort(unique(c(1, which(type.val[-1] !=
        type.val[-t]) + 1)))
    block.start.ids <- block.start.ids[block.start.ids <= nrow(blocks)]
    block.block.num <- length(block.start.ids)
    blocks[-block.start.ids, "start"] <- NA
    blocks[block.start.ids, c("end", "num")] <- t(sapply(1:block.block.num,
        function(i) {
            ids <- block.start.ids[i]:ifelse(i >= block.block.num,
                nrow(blocks), block.start.ids[i + 1] - 1)
            c(blocks[ids[length(ids)], "end"], sum(blocks[ids,
                "num"], na.rm = TRUE))
        }))
    blocks <- rbind(blocks[!is.na(blocks[, "start"]), ])
    blocks[, "size"] <- blocks[, "end"] - blocks[, "start"] +
        1
    blocks
}

`findBlockAndFilter` <-
function (x, base.position = 1:length(x), size = sum(blocks[, 
    "size"], na.rm = TRUE)/100, num = sum(blocks[, "num"], na.rm = TRUE)/100, 
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
        block.ids <- which(is.na(blocks[, "type"]) & blocks[, 
            "size"] > -1)
        block.ids.next <- block.ids[block.ids < nrow(blocks)]
        block.ids.pre <- block.ids[block.ids > 1]
        if (length(block.ids.next) > 0) 
            blocks[block.ids.next, "end"] <- blocks[block.ids.next + 
                1, "start"] - 1
        if (length(block.ids.pre) > 0) 
            blocks[block.ids.pre, "start"] <- blocks[block.ids.pre - 
                1, "end"] + 1
        blocks[, "size"] <- blocks[, "end"] - blocks[, "start"] + 
            1
    }
    blocks.margin <- block_start[-1] - block_end[-length(block_end)] - 
        1
    blocks.extentSize <- as.numeric(blocks[, "size"]) + c(blocks.margin, 
        0) + c(0, blocks.margin)
    filter.block <- (blocks.extentSize < size & as.numeric(blocks[, 
        "num"]) < num)
    if (sum(filter.block, na.rm = T) > 0) {
        blocks[filter.block, "type"] <- NA
    }
    blocks <- mergeBlocks(blocks)
    if (nrow(blocks) == 1) 
        return(blocks)
    if (fillSmallNA) {
        blocks.margin <- blocks[-1, "start"] - blocks[-nrow(blocks), 
            "end"] - 1
        blocks.extentSize <- as.numeric(blocks[, "size"]) + c(blocks.margin, 
            0) + c(0, blocks.margin)
        filter.block <- blocks.extentSize < size
        block.na.ids <- which(is.na(blocks[, "type"]) & filter.block)
        block.na.num <- length(block.na.ids)
        if (block.na.num > 0) {
            blocks[block.na.ids, "type"] <- sapply(block.na.ids, 
                function(i) {
                  if (i == 1) 
                    return(blocks[2, "type"])
                  if (i == nrow(blocks)) 
                    return(blocks[nrow(blocks) - 1, "type"])
                  if (blocks[i - 1, "type"] == blocks[i + 1, 
                    "type"]) 
                    return(blocks[i - 1, "type"])
                  blocks[i, "type"]
                })
        }
    }
    blocks <- mergeBlocks(blocks)
    if (nrow(blocks) == 1) 
        return(blocks)
    block.ids <- which(is.na(blocks[, "type"]))
    block.ids.next <- block.ids[block.ids < nrow(blocks)]
    block.ids.pre <- block.ids[block.ids > 1]
    if (length(block.ids.next) > 0) 
        blocks[block.ids.next, "end"] <- blocks[block.ids.next + 
            1, "start"] - 1
    if (length(block.ids.pre) > 0) 
        blocks[block.ids.pre, "start"] <- blocks[block.ids.pre - 
            1, "end"] + 1
    tmp.blocks <- blocks
    for (i in 2:nrow(blocks)) if (blocks[i, 1] - blocks[i - 1, 
        2] > 1) 
        tmp.blocks <- rbind(tmp.blocks, c(blocks[i - 1, "end"] + 
            1, blocks[i, "start"] - 1, blocks[i, "start"] - blocks[i - 
            1, "end"] - 1, 0, NA))
    blocks <- tmp.blocks[order(tmp.blocks[, "end"]), ]
    mergeBlocks(blocks)
}

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
            cat("chr: ", i.chr, "\tline: ", i.line <<- i.line + 
                1, "\t", colnames(genoData)[i.line], "\r")
            x.nna <- !is.na(x)
            if (corrected == FALSE) 
                x.correct <- correct.FUN(x[x.nna], base.position[ids][x.nna], 
                  ...)
            else x.correct <- x[x.nna]
            blocks.mat <- findBlockAndFilter(x.correct, base.position = base.position[ids][x.nna], 
                size = size, num = num, fillSmallNA = fillSmallNA)
            blocks.mat[nrow(blocks.mat), "end"] <- max(base.position[ids])
            blocks.mat[nrow(blocks.mat), "size"] <- blocks.mat[nrow(blocks.mat), 
                "end"] - blocks.mat[nrow(blocks.mat), "start"] + 1
            if (heterozygote == FALSE) 
                blocks.mat[blocks.mat[, "type"] == .MPR_hetero_, 
                  "type"] <- NA
            blocks.mat <- mergeBlocks(blocks.mat)
            t(blocks.mat)
        })
        blocks.mat <- matrix(unlist(blocks, recursive = TRUE, 
            use.names = FALSE), ncol = 5, byrow = TRUE)
        bin.border <- sort(unique(as.numeric(blocks.mat[, 2])))
        cat("\rchr: ", i.chr, "\tTotal", length(bin.border), 
            "borders. ")
        bin.border <- sort(unique(bin.border[c(1, which(diff(bin.border) >= 
            minBinsize) + 1)]))
        cat(length(bin.border), "borders after filtering out bins less than", 
            round(minBinsize/1000, 1), "kb.\n")
        geno.bin <- sapply(blocks, function(blocks.line) {
            blocks.end <- rbind(blocks.line)[2, ]
            ids <- match(bin.border, blocks.end)
            ids[is.na(ids)] <- findInterval(bin.border[is.na(ids)], 
                blocks.end, rightmost.closed = FALSE) + 1
            filter <- ids > ncol(blocks.line)
            ids[filter] <- ncol(blocks.line)
            rbind(blocks.line)[5, ids]
        })
        rownames(geno.bin) <- bin.border
        list(block = blocks, bin = geno.bin, border = bin.border)
    })
    geno.bin <- NULL
    for (i in 1:length(res)) geno.bin <- rbind(geno.bin, res[[i]][[2]])
    cat("Done.\n")
    list(block = res, bin = geno.bin, border = as.numeric(rownames(geno.bin)))
}

theBin=genoToBin(geno.data,as.numeric(rownames(geno.data)),corrected=F,heterozygote=T,minBinsize=250000)

#geno.bin=theBin$bin

#geno.colors <- geno.bin;geno.colors[is.na(geno.colors)] <- rgb(0,0,0)
#geno.colors[geno.colors==0]=rgb(0,0.7,0)
#geno.colors[geno.colors==1]=rgb(0,0,0.7)
#geno.colors[geno.colors==0.5]=rgb(0.6,0,0)
#rils=matrix(theBin$border,nrow=nrow(geno.bin),ncol=ncol(geno.bin))	# 1:nrow(geno.bin) or theBin$border
#poses=matrix(rep(1:ncol(geno.bin),each=nrow(geno.bin)),nrow=nrow(geno.bin),ncol=ncol(geno.bin))
#png(paste(argv[2],'.png',sep=''), 1600, 1200, pointsize=24,res=96,bg='white')
#plot(poses, rils, col=geno.colors,pch=20, ylab='',xlab="RIL index",mar=c(1.1,2.1,2.1,1.1))
#dev.off()

write.table(theBin$bin,paste(argv[2],'.bin',sep=''),col.names=T,quote=F,sep='\t')
write.table(theBin$border,paste(argv[2],'.border',sep=''),col.names=F,row.names=F,quote=F,sep='\t')
write.table(theBin$block$'1'$block,paste(argv[2],'.block',sep=''),col.names=F,quote=F,sep='\t')

