#include "functions.h"
#include <htslib/sam.h>

// from sam.c
/*******************
 *** Memory pool ***
 *******************/
/*
typedef struct {
    int k, x, y, end;
} cstate_t;

static cstate_t g_cstate_null = { -1, 0, 0, 0 };

typedef struct __linkbuf_t {
    bam1_t b;
    int32_t beg, end;
    cstate_t s;
    struct __linkbuf_t *next;
} lbnode_t;

typedef struct {
    int cnt, n, max;
    lbnode_t **buf;
} mempool_t;

static mempool_t *mp_init(void)
{
    mempool_t *mp;
    mp = (mempool_t*)calloc(1, sizeof(mempool_t));
    return mp;
}
static void mp_destroy(mempool_t *mp)
{
    int k;
    for (k = 0; k < mp->n; ++k) {
        free(mp->buf[k]->b.data);
        free(mp->buf[k]);
    }
    free(mp->buf);
    free(mp);
}
static inline lbnode_t *mp_alloc(mempool_t *mp)
{
    ++mp->cnt;
    if (mp->n == 0) return (lbnode_t*)calloc(1, sizeof(lbnode_t));
    else return mp->buf[--mp->n];
}
static inline void mp_free(mempool_t *mp, lbnode_t *p)
{
    --mp->cnt; p->next = 0; // clear lbnode_t::next here
    if (mp->n == mp->max) {
        mp->max = mp->max? mp->max<<1 : 256;
        mp->buf = (lbnode_t**)realloc(mp->buf, sizeof(lbnode_t*) * mp->max);
    }
    mp->buf[mp->n++] = p;
}
*/
/***********************
 *** Pileup iterator ***
 ***********************/
/*
// Dictionary of overlapping reads
KHASH_MAP_INIT_STR(olap_hash, lbnode_t *)
typedef khash_t(olap_hash) olap_hash_t;

struct __bam_plp_t {
    mempool_t *mp;
    lbnode_t *head, *tail;
    int32_t tid, pos, max_tid, max_pos;
    int is_eof, max_plp, error, maxcnt;
    uint64_t id;
    bam_pileup1_t *plp;
    // for the "auto" interface only
    bam1_t *b;
    bam_plp_auto_f func;
    void *data;
    olap_hash_t *overlaps;
};

// from https://github.com/samtools/samtools/blob/develop/bam_sort.c
typedef bam1_t *bam1_p;
// Function to compare reads and determine which one is < the other
static inline int bam1_lt(const bam1_p a, const bam1_p b) {
   else return (((uint64_t)a->core.tid<<32|(a->core.pos+1)<<1|bam_is_rev(a)) < ((uint64_t)b->core.tid<<32|(b->core.pos+1)<<1|bam_is_rev(b)));
}
KSORT_INIT(sort, bam1_p, bam1_lt)
*/

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
