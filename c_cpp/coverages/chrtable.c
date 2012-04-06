#include <stdlib.h>
#include "chrtable.h"
#include <err.h>
#include <zlib.h>

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

char *strlinker(const char * const main, const char * const suffix) {
    char *outstr = malloc(strlen(main)+strlen(suffix)+1);
    strcpy(outstr,main);
    strcat(outstr,suffix);
    return outstr;
}

void inc_depth(int32_t left, int32_t right, uint8_t *ThisDat) {
    //if (right <= left) errx(3,"Position Error:[%i,%i].",left,right);
    for (int32_t i=left; i<=right; ++i) {
        int tmp = ThisDat[i];
        if (tmp < UINT8_MAX) {
            ++ThisDat[i];
        }
    }
}

void do_stat(bam1_t *b, const uint16_t overlap, uint8_t **ChrDat) {
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

char tmpChr[255];
int do_contig(const uint8_t mindep, const struct myData * Data, const char * const outfile) {
    sprintf(tmpChr,"_%u.contig.gz",mindep);
    gzFile *fpd = gzopen(strlinker(outfile,tmpChr), "wb9");
    gzprintf(fpd,"#minDepth = %u\n",mindep);
    int_fast8_t inContig=0;
    for (int32_t i=0; i < Data->n_targets; ++i) {
        for (int32_t p=0; p < Data->target_len[i]; ++p) {
            if ( (Data->ChrDat[i][p] >= mindep) ) {
                if (!inContig) {
                    gzprintf(fpd,"%s\t%i\t",Data->target_name[i],p);
                    inContig=1;
                }
            } else {
                if (inContig) {
                    gzprintf(fpd,"%i\n",p);
                    inContig=0;
                }
            }
        }
    }
    gzclose(fpd);
    return 0;
}