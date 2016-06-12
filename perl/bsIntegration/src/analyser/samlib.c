#include "functions.h"
#include <htslib/sam.h>

int getPairedSam(htsFile *fp, bam1_t *b, , bam1_t *d) {
	hts_idx_t *idx;
	htsFile *fp2 = sam_open(fp->fn, "r");	// https://github.com/samtools/htslib/blob/f689e1880bec15b768d098e9e324bbaaad9ab040/hts.c#L717
	if ((idx = sam_index_load(in, fp->fn)) == 0) {
		fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
		return 1;
	}
	const bam1_core_t *c = &b->core;
	char *qname = bam_get_qname(b);
	uint16_t flag1 = c->flag & (BAM_FREAD1 | BAM_FREAD2);
	hts_itr_t *iter;
	if ((iter = sam_itr_queryi(idx, c->mtid, c->mpos, c->mpos+1)) == 0) {
		fprintf(stderr, "[E::%s] fail to parse region '%s(%d):%d'\n", __func__, h->target_name[c->mtid], c->mtid, c->mpos);
		continue;
	}
	///ARRAY
	while ((r = sam_itr_next(fp2, iter, d)) >= 0) {	// 存在 第一个左端点符合的
		uint16_t flag2 = (d->core).flag & (BAM_FREAD1 | BAM_FREAD2);
		if ( (flag1 + flag2) != (flag1 | flag2) ) continue;	// Read1 with Read2 only. Read12互补
		if (c->mpos != (d->core).pos) continue;
		char *qname2 = bam_get_qname(d);
		if (strcmp(qname,qname2) != 0) continue;	// qname一致
		///PUSH
	}
	hts_itr_destroy(iter);
	///Choose ONE
	return 0;
}