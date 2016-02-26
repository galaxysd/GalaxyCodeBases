# source("a1funs.r")

distance <- function (pa,pb) {
	differ <- pa-pb
	diss <- sum(differ*differ)
	sqrt(diss)
}

total_distance <- function (mpoints,lpairs) {
	dis_summer <- 0
	for (edge in lpairs) {
		curr_dis <- distance(mpoints[edge[1],],mpoints[edge[2],])
		dis_summer <- dis_summer + curr_dis
		#print(c(edge,curr_dis,dis_summer))
	}
	dis_summer
}

plot_pairs <- function (mpoints,lpairs,mytitle='',plot=FALSE, wait=TRUE) {
	if (plot == FALSE) return
	if (wait == TRUE) readline(prompt="Press [enter] to continue")
	plot(mpoints,type='p',xlab="x",ylab='y',pch = 19)
	for (edge in lpairs) {
		xyxy <- mpoints[edge[1:2],]
		lines(xyxy[,1],xyxy[,2])
	}
	tdis <- total_distance(mpoints,lpairs)
	title(paste('Total distance',tdis,mytitle))
}

random_pairings <- function (mpoints) {
	n <- nrow(mpoints)
	n <- as.integer(n/2)*2
	smp <- sample(n)
	foo <- matrix(smp,ncol=2)
	lapply(1:nrow(foo), function(i) foo[i,])
	# lapply(seq_len(ncol(x)), function(i) x[,i])
}

improve_pairings <- function (mpoints,inpairs,cnt=0,plot=FALSE) {
	maxcnt <- length(inpairs)
	lpairs <- inpairs[sample(maxcnt)]
	if (cnt==0 || cnt > maxcnt) cnt <- maxcnt
	i <- 1
	while (i < cnt) {
		edgeO <- lpairs[c(i,i+1)]
		abxy <- unlist(edgeO)
		#xxyy <- simplify2array(edgeO)
		axby <- list(c(abxy[c(1,3)]),abxy[c(2,4)])
		aybx <- list(c(abxy[c(1,4)]),abxy[c(2,3)])
		disO <- total_distance(mpoints,edgeO)
		disaxby <- total_distance(mpoints,axby)
		disaybx <- total_distance(mpoints,aybx)
		if (disaxby < disO) {
			lpairs[i] <- axby[1]
			lpairs[i+1] <- axby[2]
		}
		if (disaybx < disO) {
			lpairs[i] <- aybx[1]
			lpairs[i+1] <- aybx[2]
			if (disaxby < disaybx) {
				lpairs[i] <- axby[1]
				lpairs[i+1] <- axby[2]
			}
		}
		#plot_pairs(mpoints,lpairs,paste(',Time(s)',(i+1)/2),plot)
		i <- i + 2
	}
	lpairs
}

find_pairingsY <- function (mpoints, tries=1, plot=FALSE) {
	pairs0 <- random_pairings(mpoints)
	plot_pairs(mpoints,pairs0,'Init.',plot)
	cnt <- length(pairs0)
	if (tries > cnt) tries <- cnt
	if (tries < 1) tries <- 1
	new_pairs <- improve_pairings(mpoints,pairs0,tries,plot)
	plot_pairs(mpoints,new_pairs,paste('Time(s)',tries),plot)
}

find_pairings <- function (mpoints, tries=1, plot=FALSE) {
	pairs0 <- random_pairings(mpoints)
	plot_pairs(mpoints,pairs0,'Init.',plot)
	lastpairs <- pairs0
	for (i in 1:tries) {
		new_pairs <- improve_pairings(mpoints,lastpairs)
		plot_pairs(mpoints,new_pairs,paste('Time(s)',i),plot)
		lastpairs <- new_pairs
	}
}

find_pairingsX <- function (mpoints, tries=1, plot=FALSE) {
	pairs0 <- random_pairings(mpoints)
	plot_pairs(mpoints,pairs0,'Init.',plot)
	cnt <- length(pairs0)
	if (tries > cnt) tries <- cnt
	if (tries < 1) tries <- 1
	need_id <- sample(cnt)[1:tries]
	need_imp <- pairs0[need_id]
	had_imp <- improve_pairings(mpoints,need_imp)
	pairs1 <- pairs0
	for (i in 1:tries) {
		pairs1[need_id[i]] <- had_imp[i]
	}
	plot_pairs(mpoints,pairs1,paste('Tries(s)',tries),plot)
}

