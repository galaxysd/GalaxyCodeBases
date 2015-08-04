#include "BWTSearch.h"

unsigned int REVBWTForwardSearch(const unsigned char *convertedkey, const unsigned int keylength, const BWT *rev_bwt, unsigned int *resultsaindexleft, unsigned int *resultsaindexright, unsigned int *rev_resultsaindexleft, unsigned int *rev_resultsaindexright) {

	unsigned int sacount=0;
	unsigned int rev_startsaindex, rev_endsaindex;
	unsigned int startsaindex, endsaindex;
	unsigned int pos = 1;
	int i;
	unsigned int c = convertedkey[0];
	unsigned int occcount_start[4];
	unsigned int occcount_end[4];
	unsigned int occcount[4];

	rev_startsaindex = rev_bwt->cumulativeFreq[convertedkey[0]]+1;
	rev_endsaindex = rev_bwt->cumulativeFreq[convertedkey[0]+1];
	startsaindex = rev_bwt->cumulativeFreq[convertedkey[0]]+1;
	endsaindex = rev_bwt->cumulativeFreq[convertedkey[0]+1];

	while (pos < keylength && startsaindex <= endsaindex) {
		c = convertedkey[pos];

		BWTAllOccValue(rev_bwt,rev_startsaindex,occcount_start);
		BWTAllOccValue(rev_bwt,rev_endsaindex + 1,occcount_end);

		rev_startsaindex = rev_bwt->cumulativeFreq[c] + occcount_start[c] + 1;
		rev_endsaindex = rev_bwt->cumulativeFreq[c] + occcount_end[c];

		occcount[3]=0;
		for (i=2;i>=0;i--) {
			occcount[i]=occcount[i+1]+occcount_end[i+1]-occcount_start[i+1];
		}

		endsaindex = endsaindex - occcount[c];
		startsaindex = endsaindex - (rev_endsaindex-rev_startsaindex);
		pos++;
	}

	*resultsaindexleft = startsaindex;
	*resultsaindexright = endsaindex;
	*rev_resultsaindexleft = rev_startsaindex;
	*rev_resultsaindexright = rev_endsaindex;

	sacount+=endsaindex-startsaindex+1;
	// number of occurrence = endsaindex - startsaindex + 1
	return sacount;

}


unsigned int REVBWTContForwardSearch(const unsigned char *convertedkey, const unsigned int start, const unsigned int len, const BWT *rev_bwt, unsigned int *sal, unsigned int *sar, unsigned int *rev_sal, unsigned int *rev_sar) {

	unsigned int sacount=0;
	unsigned int pos = start;
	unsigned char c;
	unsigned int occcount_start[4];
	unsigned int occcount_end[4];
	unsigned int occcount[4];
	int k;
	while (pos < start+len  && *sal <= *sar) {
		c = convertedkey[pos];

		BWTAllOccValue(rev_bwt,*rev_sal,occcount_start);
		BWTAllOccValue(rev_bwt,*rev_sar + 1,occcount_end);

		*rev_sal = rev_bwt->cumulativeFreq[c] + occcount_start[c] + 1;
		*rev_sar = rev_bwt->cumulativeFreq[c] + occcount_end[c];

		occcount[3]=0;
		for (k=2;k>=0;k--) {
			occcount[k]=occcount[k+1]+occcount_end[k+1]-occcount_start[k+1];
		}

		*sar = *sar - occcount[c];
		*sal = *sar - (*rev_sar-*rev_sal);

		pos++;
	}
	sacount+=*sar-*sal+1;
	return sacount;

}
unsigned int BWTContBackwardSearch(const unsigned char *convertedkey, const unsigned int start, const unsigned int len, const BWT *bwt, unsigned int *sal, unsigned int *sar) {

	unsigned int sacount=0;
	unsigned int pos = len;
	unsigned char c;

	if (*sal > *sar) {
		return 0;
	}

	while (pos > 0 && *sal <= *sar) {
		c = convertedkey[pos-1];
		*sal = bwt->cumulativeFreq[c] + BWTOccValue(bwt, *sal, c) + 1;
		*sar = bwt->cumulativeFreq[c] + BWTOccValue(bwt, *sar + 1, c);
		pos--;
	}
    sacount+=*sar-*sal+1;
	return sacount;
}

unsigned int BWTBackward1Error(char *querypattern, int chain, BWT *bwt, unsigned int start, unsigned int len, unsigned int pl, unsigned int pr, unsigned int allele2, HitInfo *hits,  unsigned int *numOfHits) {
		unsigned int mk_l=1,mk_r=0;
		unsigned int occcount_pstart[4];
		unsigned int occcount_pend[4];
		unsigned int occcount_start[4];
		unsigned int occcount_end[4];
		unsigned int occcount[4];
		unsigned int sacount=0;

		unsigned char c;
		unsigned char ec;
		int i;
		unsigned int allele1;
		//printf("bwtbackward1error %u %u\n",pl,pr);
		//                                  v--start       v----start+len
		// querypattern = xxxxxxxxxxxxxxxxx[xxxxxxxxxxxxxx]xxxxxx
		//                                      <------- search direction
		//                               querypattern[start+len-1], querypattern[start+len-2]...querypattern[start]
		//                           for i=0 to len-1,
		//                                append querypatter[start+len-1-i]!

		for (i=0;(i<len && pl<=pr);i++) {
			//call once only proceduressssss - great
			BWTAllOccValue(bwt,pl,occcount_pstart);
			BWTAllOccValue(bwt,pr + 1,occcount_pend);

			//backward manner
			for (ec=0;ec<4;ec++) {
				if (querypattern[start+len-1-i]==ec)
					continue;

				allele1= ((((start+len-1-i) & 0xfff)<<2) | (ec & 0x3));
				mk_l=pl;
				mk_r=pr;

				unsigned int pos = i+1;
				mk_l = bwt->cumulativeFreq[ec] + occcount_pstart[ec] + 1;
				mk_r = bwt->cumulativeFreq[ec] + occcount_pend[ec];

				if (BWTContBackwardSearch(querypattern,start,len-i-1,bwt,&mk_l,&mk_r)) {
//					printf("%d\t", start+len-1-i);
//					printf("%d\n", ec);
//					printf("%d\t%d\n",allele1, allele2);
					OCCProcess(mk_l,mk_r, chain, allele1, allele2, hits, numOfHits);
					sacount+=mk_r-mk_l+1;
				}
			}
			c = querypattern[start+len-1-i];
			pl = bwt->cumulativeFreq[c] + occcount_pstart[c] + 1;
			pr = bwt->cumulativeFreq[c] + occcount_pend[c];

		}
		return sacount;
}


unsigned int REVBWTForward1Error(char *querypattern, int chain, BWT * bwt,BWT * rev_bwt, unsigned int start,unsigned int len,unsigned int pl,unsigned int pr,unsigned int rev_pl,unsigned int rev_pr, unsigned int allele2, HitInfo *hits, unsigned int *numOfHits) {

		unsigned int mk_l=1,mk_r=0,rev_mk_l,rev_mk_r;
		unsigned int occcount_pstart[4];
		unsigned int occcount_pend[4];
		unsigned int occcountp[4];
		unsigned int occcount_start[4];
		unsigned int occcount_end[4];
		unsigned int occcount[4];
		unsigned int sacount=0;

		unsigned char c;
		unsigned char ec;
		unsigned int i;
		int k;
		unsigned int allele1;
		for (i=0;(i<len && pl<=pr);i++) {
			//call once only proceduressssss - great
			BWTAllOccValue(rev_bwt,rev_pl,occcount_pstart);
			BWTAllOccValue(rev_bwt,rev_pr + 1,occcount_pend);

			occcountp[3]=0;
			for (k=2;k>=0;k--) {
				occcountp[k]=occcountp[k+1]+occcount_pend[k+1]-occcount_pstart[k+1];
			}


			//forward manner
			for (ec=0;ec<4;ec++) {
				if (querypattern[start+i]==ec)
					continue;

				allele1= ((((start+i) & 0xfff)<<2) | (ec & 0x3));
				mk_l=pl;
				mk_r=pr;
				rev_mk_l=rev_pl;
				rev_mk_r=rev_pr;

				unsigned int pos = i+1;

				rev_mk_l = rev_bwt->cumulativeFreq[ec] + occcount_pstart[ec] + 1;
				rev_mk_r = rev_bwt->cumulativeFreq[ec] + occcount_pend[ec];

				mk_r = mk_r - occcountp[ec];
				mk_l = mk_r - (rev_mk_r-rev_mk_l);

				if (REVBWTContForwardSearch(querypattern,start+pos,len-i-1,rev_bwt,&mk_l,&mk_r,&rev_mk_l,&rev_mk_r)) {
//					printf("%d\t", start+i);
//					printf("%d\n", ec);
					OCCProcess(mk_l,mk_r, chain, allele1, allele2, hits, numOfHits);
					//return mk_l, mk_r
					sacount+=mk_r-mk_l+1;
				}
			}
			c = querypattern[start+i];

			rev_pl = rev_bwt->cumulativeFreq[c] + occcount_pstart[c] + 1;
			rev_pr = rev_bwt->cumulativeFreq[c] + occcount_pend[c];

			pr = pr - occcountp[c];
			pl = pr - (rev_pr-rev_pl);
		}
		return sacount;
}

