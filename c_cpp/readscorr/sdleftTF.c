#include <stdint.h> //uint64_t
#include "sdleft.h"
#include "sdleftTF.h"

#ifndef PUBLIC
    #ifdef USEUINT16
    uint16_t  *dleft_statu16
    #elif defined USEUINT32
    uint32_t  *dleft_statu32
    #elif defined USEUINT64
    uint64_t  *dleft_statu64
    #else
        #error Must define USEUINT{16,32,64} or PUBLIC before compilering !
    #endif
     (SDLeftArray_t *dleftobj) {
        return 0;
    }
#else

void *pf = (G_SDLeftArray_IN) dleft_statu16;
#endif
