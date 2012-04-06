#include <stdlib.h>
#include "chrtable.h"
#include <err.h>

//extern struct ChrData_hash_struct *ChrData;
/*
void add_chr(char *id, uint32_t len) {
    struct ChrData_hash_struct *s;

    s = malloc(sizeof(struct ChrData_hash_struct));
    s->ChrID = malloc(strlen(id)+1);
    strcpy(s->ChrID, id);
    s->ChrLen = len;
    s->Data = malloc(len*sizeof(uint8_t));
    HASH_ADD_STR( ChrData, ChrID, s );  // id: name of key field
}
*/

void inc_depth(int32_t left, int32_t right, uint8_t *ThisDat) {
    //if (right <= left) errx(3,"Position Error:[%i,%i].",left,right);
    for (int32_t i=left; i<=right; ++i) {
        int tmp = ThisDat[i];
        if (tmp < UINT8_MAX) {
            ++ThisDat[i];
        }
    }
}

void do_stat(bam1_t *b, uint16_t overlap, uint8_t **ChrDat) {
    uint16_t k = overlap + 1;
    //bam1_t *b = balignd;
    //bam1_core_t *core = &(b->core);
    if (b->core.tid < 0) return;
    uint8_t *ThisDat = ChrDat[b->core.tid];
    int32_t left = b->core.pos;
    int32_t right = bam_calend(&b->core, bam1_cigar(b));
    if (left >= right) return;
//printf("%i %i %i\n",b->core.tid,left,right);
    inc_depth(left, right - k+1, ThisDat);
}

// b->core.pos, bam_calend(&b->core, bam1_cigar(b))