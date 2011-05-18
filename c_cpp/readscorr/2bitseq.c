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

uint64_t *dibmalloc(size_t len){
    size_t needtomallocQQW = (len+31u)>>5;
    uint64_t *outseq = malloc(needtomallocQQW*8);
    return outseq;
}

/*
|000000000>||111111111>||2222222>...|, lastBits=7, nullBits=3 (%10)
|...0000000||00>1111111||11>2222222>|. Well, left is right here ...
|11>2222222>||00>1111111||...0000000|
|<2222222<11||1111111<00||0000000...|
*/
uint64_t *dibrevcomp(uint64_t const *const inseq, size_t len){
    size_t needtomallocQQW = (len+31u)>>5;
    uint64_t *outseq = malloc(needtomallocQQW*8);
    char lastBits = 2*(len % 64);   // well "% 64" will also be "andl $63,%esi"
    char nullBits = 64-lastBits;  //(-lastBits) % 64; only use if lastBits.
    uint64_t tmpstr=0;
    uint64_t highmask = ~((1LLU<<lastBits)-1U);
//printf("QQW:%zu,L:%d,N:%d,M:%lx\n",needtomallocQQW,lastBits,nullBits,highmask);    
    for (size_t i=0;i<needtomallocQQW;i++) {    // from 0
        *(outseq+needtomallocQQW-1-i)=unitReverseComp((*(inseq+i)<<nullBits) | tmpstr);
        //*(outseq+needtomallocQQW-1-i)=(*(inseq+i)<<nullBits) | tmpstr;
        if (lastBits) { // this if will be optimized outsides of for, so I can have a rest here.
            tmpstr = (*(inseq+i) & highmask)>>lastBits;
        }
    }
    return outseq;
}

