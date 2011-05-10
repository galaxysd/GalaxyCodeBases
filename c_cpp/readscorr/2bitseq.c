#include "2bitseq.h"

char *unit2basechr(uint64_t seq32mer){
    char *basechar=malloc(33);
    char *tmpbc;
    tmpbc=basechar;
    int i;
    for(i=62;i>=0;i-=2) {
        *tmpbc++ = DBIT2BASE((seq32mer&(3<<i))>>i);
        //*tmpbc++ = '0' + ((seq32mer&(3<<i))>>i);
    }
    *tmpbc++=0;
    return basechar;
}

