#include <stdint.h> //uint64_t
#include "sdleft.h"
#include "sdleftTF.h"

#define TOKENPASTE(x, y) x ## y
#define TOKENPASTE2(x, y) TOKENPASTE(x, y)

#ifdef USEUINT16
    #define THETYPE uint16_t
#elif defined USEUINT32
    #define THETYPE uint32_t
#elif defined USEUINT64
    #define THETYPE uint64_t
#elif defined PUBLIC
    // NA
#else
    #error Must define USEUINT{16,32,64} or PUBLIC before compilering !
#endif

#define ADDSUFFIX(x) TOKENPASTE2(x ## _,THETYPE)

uint16_t *dleft_stat_uint16_t(SDLeftArray_t *dleftobj);
uint32_t *dleft_stat_uint32_t(SDLeftArray_t *dleftobj);
uint64_t *dleft_stat_uint64_t(SDLeftArray_t *dleftobj);

#ifndef PUBLIC  // USEUINT{16,32,64}
THETYPE * ADDSUFFIX(dleft_stat) (SDLeftArray_t *dleftobj) {
        return 0;
}
#else   // PUBLIC
G_SDLeftArray_IN pf = (G_SDLeftArray_IN) dleft_stat_uint16_t;
#endif
