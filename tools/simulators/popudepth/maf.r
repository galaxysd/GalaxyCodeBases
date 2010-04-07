# generate genotype for each individual
genotype=function(nsample,MAF){ # nsample为多少个个体，MAF为一维数组，存储所有SNP位点的MAF值
	nsite=length(MAF)
	a=array(runif(nsite*nsample*2,0,1),dim=c(nsite,2*nsample))
	b=a
	# generate genotype for each allele
	for (i in 1:nsite){
		for (j in 1:(2*nsample)){
			if(a[i,j]>MAF[i]){ b[i,j] = 0 }
			else{ b[i,j] = 1 }
		}
	}
	# generate genotype for each individual: 00->0;11->3;01->1,10->2
	c=array(dim=c(nsite,nsample+2))
	for (i in 1:nsite){
		c[i,1] = i # site ID
		c[i,2] = MAF[i] # MAF value
		for (j in seq(from=1,to=(2*nsample),by=2)){
			c[i,(j+1)/2+2]=b[i,j]*2+b[i,j+1]
			#if(b[i,j]==0 & b[i,j+1]==0){ c[i,(j+1)/2+2] = 0 }
			#else if(b[i,j]==1 & b[i,j+1]==1){ c[i,(j+1)/2+2] = 3 }
			#else if(b[i,j]==0 & b[i,j+1]==1){ c[i,(j+1)/2+2] = 1 }
			#else { c[i,(j+1)/2+2] = 2 }
		}
	}
	return (c)
}

genotype(5,c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))
