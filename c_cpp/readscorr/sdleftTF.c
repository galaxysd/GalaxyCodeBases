#include <stdint.h> //uint64_t
#include <stdlib.h> //calloc
#include "sdleft.h"
#include "sdleftTF.h"
#include <stdio.h>

#define TOKENPASTE(x, y) x ## y
#define TOKENPASTE2(x, y) TOKENPASTE(x, y)
#define CAT0(x) #x
#define CAT(x) CAT0(x)

#ifdef USEUINT16
    #define THETYPE uint16_t
    #define DILTYPE uint32_t
    #define FLOTYPE double
#elif defined USEUINT32
    #define THETYPE uint32_t
    #define DILTYPE uint64_t
    #define FLOTYPE double
#elif defined USEUINT64
    #define THETYPE uint64_t
    #define DILTYPE uint128_t
    #define FLOTYPE float128
#elif defined PUBLIC
    // NA
#else
    #error Must define USEUINT{16,32,64} or PUBLIC before compilering !
#endif

#define ADDSUFFIX(x) TOKENPASTE2(x ## _,THETYPE)

SDLeftStat_t *dleft_stat_uint16_t(SDLeftArray_t * const dleftobj, FILE *);
SDLeftStat_t *dleft_stat_uint32_t(SDLeftArray_t * const dleftobj, FILE *);
SDLeftStat_t *dleft_stat_uint64_t(SDLeftArray_t * const dleftobj, FILE *);

#ifndef PUBLIC  // USEUINT{16,32,64}
SDLeftStat_t * ADDSUFFIX(dleft_stat) (SDLeftArray_t * const dleftobj, FILE *stream) {
fprintf(stderr,"[!]Count with[%s]\n",CAT(THETYPE));
    SDLeftStat_t *pSDLeftStat = malloc(sizeof(SDLeftStat_t));
    THETYPE * const pCountHistArray = calloc(sizeof(THETYPE),1+dleftobj->maxCountSeen);
    size_t totalDLAsize = SDLA_ITEMARRAY*dleftobj->itemByte*dleftobj->ArraySize;
    //size_t firstlevelDLAitemsize = SDLA_ITEMARRAY*dleftobj->itemByte;
    const unsigned char * const pDLA = dleftobj->pDLA;
    THETYPE Item_CountBits=0;  // set value in case SDLA_ITEMARRAY*dleftobj->itemByte == 0 (EP ?)
    uint128_t theItem;
    for (size_t i=0;i<totalDLAsize;i+=dleftobj->itemByte) {
        //const unsigned char * pChunk = pDLA + i;
        theItem = 0;
        for (uint_fast8_t j=0;j<dleftobj->itemByte;j++) {
            theItem |= ((uint128_t)*(pDLA + i + j)) << (j*8u);
        }
        Item_CountBits = theItem & dleftobj->Item_CountBitMask;
        ++pCountHistArray[Item_CountBits];
        //++CountSum;
    }
    THETYPE CountSum=0;
    DILTYPE CountSumSquared=0;
    double SStd;    // We need to return in a fixed type for printf
    for (size_t p=1;p<=dleftobj->maxCountSeen;p+=sizeof(size_t)) {
        CountSum += pCountHistArray[p];
        CountSumSquared += pCountHistArray[p] * pCountHistArray[p];
    }
    //http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    SStd = ( (FLOTYPE)CountSumSquared - ((FLOTYPE)CountSum*(FLOTYPE)CountSum/(FLOTYPE)dleftobj->ItemInsideAll) )
            / (dleftobj->ItemInsideAll -1);
    pSDLeftStat->SStd = SStd;
    pSDLeftStat->Mean = (double)CountSum / (double)dleftobj->ItemInsideAll;
    free(pCountHistArray);
    // deal(*pextree);
    return pSDLeftStat;
}
#else   // PUBLIC
//G_SDLeftArray_IN *pf = dleft_stat_uint16_t;
SDLeftStat_t * dleft_stat(SDLeftArray_t * const dleftobj, FILE *stream) {
    G_SDLeftArray_IN *pf;
    if (dleftobj->ItemInsideAll <= UINT16_MAX) {
        pf = dleft_stat_uint16_t;
    } else if (dleftobj->ItemInsideAll <= UINT32_MAX) {
        pf = dleft_stat_uint32_t;
    } else {    // uint64_t ItemInsideAll
        pf = dleft_stat_uint64_t;
    }
    return (*pf)(dleftobj,stream);
}
#endif
