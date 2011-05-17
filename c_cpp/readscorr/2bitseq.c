#include "2bitseq.h"
#include <string.h> // strcpy

char *unit2basechr(uint64_t seq32mer){
    char *basechar=malloc(33);
    char *tmpbc;
    tmpbc=basechar;
    int i;
    for(i=0;i<64;i+=2) {
        *tmpbc++ = DBIT2BASE((seq32mer&(3LLU<<i))>>i);
        //*tmpbc++ = '0' + ((seq32mer&(3u<<i))>>i);
        //printf("=%lx,%lx,%lx,%lx=\n",seq32mer,3LLU<<i,seq32mer&(3LLU<<i),(seq32mer&(3LLU<<i))>>i);
    }
    *tmpbc++=0;
    return basechar;
}

char *dib2basechr(uint64_t *diBseq, size_t len){
    char *basechar=malloc(len+1);
    char *tmpstr;
    size_t i;
    for (i=0;i<len/32;i++) {
        tmpstr=unit2basechr(diBseq[i]);
        strcpy(basechar+i*32,tmpstr);
        free(tmpstr);
    }
    tmpstr=unit2basechr(diBseq[i]);
    strncpy(basechar+i*32,tmpstr,len % 32);
    free(tmpstr);
    *(basechar+len)='\0';
    return basechar;
}

/*
|--------->|    |------>nnn|, lastBits=7, remainBits=3 (%10)
|   -------|--> |   ------>|. Well, left is right here ...
*/
uint64_t *dibrevcomp(uint64_t const *const inseq, size_t len){
    size_t needtomallocQQW = (len+31u)>>5;
    uint64_t *outseq = malloc(needtomallocQQW*8);
    char lastBits = len % 64;   // well "% 64" will also be "andl $63,%esi"
    char remainBits = (-lastBits) % 64;
    uint64_t tmpstr=0;
    uint64_t highmask = ~((1LLU<<remainBits)-1U);
    for (size_t i=0;i<needtomallocQQW;i++) {    // from 0
        *(outseq+needtomallocQQW-1-i)=unitReverseComp((*(inseq+i)<<remainBits) | tmpstr);
        if (lastBits) { // this if will be optimized outsides of for, so I can have a rest here.
            tmpstr = (*(inseq+i) & highmask)>>lastBits;
        }
    }
    return outseq;
}

