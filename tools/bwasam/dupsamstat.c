#include <stdlib.h> //EXIT_FAILURE
#include <stdio.h>
//#include <float.h>
#include <math.h>
#include <search.h> //hsearch
#include "bam/sam.h"
#include "bitarray.h"

/*
typedef struct {
	int beg, end;
	samfile_t *in;
} tmpstruct_t;

// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
	bam_plbuf_t *buf = (bam_plbuf_t*)data;
	bam_plbuf_push(b, buf);
	return 0;
}
// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
	tmpstruct_t *tmp = (tmpstruct_t*)data;
	if ((int)pos >= tmp->beg && (int)pos < tmp->end)
		printf("%s\t%d\t%d\n", tmp->in->header->target_name[tid], pos + 1, n);
	return 0;
}
*/
typedef struct {
    void *pBitArray;
    uint32_t ChrLen;
} hData_t;

int32_t ChrNum;

int main (int argc, char *argv[]) {
	if (argc!=3+1) {
		fprintf(stderr,"Usage: %s <in.sam> <snp.lst> <outprefix>\n",argv[0]);
		exit(EXIT_FAILURE);
	}
	samfile_t *samIn = samopen(argv[1], "r", NULL);
	if (samIn == 0) {
		fprintf(stderr, "Fail to open SAM file %s\n", argv[1]);
		exit(EXIT_FAILURE);
	}
	bam_header_t *samHeader=samIn->header;
	ChrNum = samHeader->n_targets;
	hcreate(2*ChrNum);
	ENTRY hItem, *phItem;
	hData_t *hData = malloc(ChrNum * sizeof(hData_t));
    if (hData == NULL) {
       fprintf(stderr, "[x]malloc failed\n");
       exit(EXIT_FAILURE);
    }
	for (int32_t i=0;i<ChrNum;++i) {
	    printf("[%0*d] %21s -> %9u bp\n",(int)(1+log10(ChrNum)),1+i,samHeader->target_name[i],samHeader->target_len[i]);
	    hItem.key  = samHeader->target_name[i];
	    hData[i].ChrLen = samHeader->target_len[i];
	    hData[i].pBitArray = calloc(1, BITNSLOTS(samHeader->target_len[i]) );
	    hItem.data = &hData[i];
        phItem = hsearch(hItem, ENTER);
        /* there should be no failures */
        if (phItem == NULL) {
           fprintf(stderr, "[x]hash entry failed\n");
           exit(EXIT_FAILURE);
        }

	}
/*
	if (argc == 2) { // if a region is not specified
		sampileup(samIn, -1, pileup_func, &tmp);
	} else {
		int ref;
		bam_index_t *idx;
		bam_plbuf_t *buf;
		idx = bam_index_load(argv[1]); // load BAM index
		if (idx == 0) {
			fprintf(stderr, "BAM indexing file is not available.\n");
			return 1;
		}
		bam_parse_region(samIn->header, argv[2], &ref,
		                 &tmp.beg, &tmp.end); // parse the region
		if (ref < 0) {
			fprintf(stderr, "Invalid region %s\n", argv[2]);
			return 1;
		}
		buf = bam_plbuf_init(pileup_func, &tmp); // initialize pileup
		bam_fetch(samIn->x.bam, idx, ref, tmp.beg, tmp.end, buf, fetch_func);
		bam_plbuf_push(0, buf); // finalize pileup
		bam_index_destroy(idx);
		bam_plbuf_destroy(buf);
	}
*/
	samclose(samIn);
	for (int32_t i=0;i<ChrNum;++i) {
	    free(hData[i].pBitArray);
	}
	free(hData);
	return 0;
}
