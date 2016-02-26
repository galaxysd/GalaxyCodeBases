# source("a1funs.r")

find_pairings <- function () {
	;
}

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

plot_pairs <- function (mpoints,lpairs) {
	plot(mpoints,type='p',xlab="x",ylab='y',pch = 19)
	for (edge in lpairs) {
		xyxy <- mpoints[edge[1:2],]
		lines(xyxy[,1],xyxy[,2])
	}
}

random_pairings <- function (mpoints) {
	
}