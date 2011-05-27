#include <stdint.h> //uint64_t
#include <stdlib.h> //calloc
#include "sdleft.h"
#include "sdleftTF.h"

#define TOKENPASTE(x, y) x ## y
#define TOKENPASTE2(x, y) TOKENPASTE(x, y)

#ifdef USEUINT16
    #define THETYPE uint16_t
    #define DILTYPE uint32_t
#elif defined USEUINT32
    #define THETYPE uint32_t
    #define DILTYPE uint64_t
#elif defined USEUINT64
    #define THETYPE uint64_t
    #define DILTYPE uint128_t
#elif defined PUBLIC
    // NA
#else
    #error Must define USEUINT{16,32,64} or PUBLIC before compilering !
#endif

#define ADDSUFFIX(x) TOKENPASTE2(x ## _,THETYPE)

SDLeftStat_t *dleft_stat_uint16_t(SDLeftArray_t * const dleftobj);
SDLeftStat_t *dleft_stat_uint32_t(SDLeftArray_t * const dleftobj);
SDLeftStat_t *dleft_stat_uint64_t(SDLeftArray_t * const dleftobj);

#ifndef PUBLIC  // USEUINT{16,32,64}
SDLeftStat_t * ADDSUFFIX(dleft_stat) (SDLeftArray_t * const dleftobj) {
    SDLeftStat_t *pSDLeftStat = malloc(sizeof(SDLeftStat_t));
    size_t * const pCountHistArray = calloc(sizeof(size_t),1+dleftobj->maxCountSeen);
    size_t totalDLAsize = SDLA_ITEMARRAY*dleftobj->itemByte*dleftobj->ArraySize;
    //size_t firstlevelDLAitemsize = SDLA_ITEMARRAY*dleftobj->itemByte;
    const unsigned char * const pDLA = dleftobj->pDLA;
    uint64_t Item_CountBits=0;  // set value in case SDLA_ITEMARRAY*dleftobj->itemByte == 0 (EP ?)
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
    uint64_t CountSum=0;
    uint128_t CountSumSquared=0;
    for (size_t p=1;p<=dleftobj->maxCountSeen;p+=sizeof(size_t)) {
        CountSum += pCountHistArray[p];
        CountSumSquared += pCountHistArray[p] * pCountHistArray[p];
    }
    free(pCountHistArray);
    // deal(*pextree);
    return pSDLeftStat;
}
#else   // PUBLIC
G_SDLeftArray_IN pf = dleft_stat_uint16_t;
#endif
