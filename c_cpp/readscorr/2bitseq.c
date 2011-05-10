#include "2bitseq.h"

unsigned char *unit2basechr(uint64_t seq32mer){
    unsigned char *basechar=malloc(33);
    unsigned char *tmpbc;
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

