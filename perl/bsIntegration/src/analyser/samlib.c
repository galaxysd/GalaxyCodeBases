#include "functions.h"
#include <htslib/sam.h>

int getPairedSam(htsFile *fp, hts_idx_t *idx, bam1_t *b, bam1_t *d) {
	int r, ret = 1;
	static htsFile *fp2;
	if (idx==NULL) {
		ret = 3;
		r = sam_close(fp2);
		if (r < 0) {
			fprintf(stderr, "Error closing input.\n");
			ret = 1;
		}
		return ret;
	}
	if (fp2) {
		if (strcmp(fp->fn,fp2->fn) != 0) {	// https://github.com/samtools/htslib/blob/f689e1880bec15b768d098e9e324bbaaad9ab040/hts.c#L717
			r = sam_close(fp2);
			if (r < 0) {
				fprintf(stderr, "Error closing input.\n");
				ret = 1;
			}
			fp2 = sam_open(fp->fn, "r");
		}
	} else {
		fp2 = sam_open(fp->fn, "r");
	}
	const bam1_core_t *c = &b->core;
	if (c->mtid < 0) {
		return 1;
	}
	char *qname = bam_get_qname(b);
	uint16_t flag1 = c->flag & (BAM_FREAD1 | BAM_FREAD2);
	hts_itr_t *iter;
	if ((iter = sam_itr_queryi(idx, c->mtid, c->mpos, c->mpos+1)) == 0) {
		fprintf(stderr, "[E::%s] fail to parse region '(%d):%d'\n", __func__, c->mtid, c->mpos);
		ret = 1;
	}
	while ((r = sam_itr_next(fp2, iter, d)) >= 0) {	// 存在 第一个左端点符合的 && Read1 with Read2 only. Read12互补 && qname一致
		uint16_t flag2 = (d->core).flag & (BAM_FREAD1 | BAM_FREAD2);
		char *qname2 = bam_get_qname(d);
		if ( ((flag1 + flag2) == (flag1 | flag2)) && (c->mpos == (d->core).pos) && (strcmp(qname,qname2) == 0) ) {
			ret = 0;	// 理论上符合条件的只会有一个。不深入hts_itr_next的话，就这样了。
			goto final;
		}
	}
	final:
	hts_itr_destroy(iter);
	return ret;
}
