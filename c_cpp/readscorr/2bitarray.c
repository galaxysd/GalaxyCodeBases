#include "2bitarray.h"

FORCE_INLINE void TBITSetValue ( char arr[], size_t index, uint_fast8_t value ) {
    if (value>3) value=3;
    unsigned char theByte=(arr[TBITSLOT(index)]) & ~TBITMASK(index);
    theByte |= value << TBITSHIFT(index);
    arr[TBITSLOT(index)] = theByte;
}

FORCE_INLINE uint_fast8_t TBITSaturatedADD ( char arr[], size_t index, int_fast8_t value ) {
    unsigned char theByte = (arr[TBITSLOT(index)]);
    int_fast8_t theValue = (theByte >> TBITSHIFT(index)) & 3 ;
    if (value==0) {
        return (uint_fast8_t) theValue;
    } else if (value > 0) {
        if (theValue > 3-value) theValue=3;
          else theValue += value;
    } else {
        if (theValue < -value) theValue =0;
          else theValue += value;
    }
    theByte &= ~TBITMASK(index);
    theByte |= theValue << TBITSHIFT(index);
    arr[TBITSLOT(index)] = theByte;
    return (uint_fast8_t) theValue;
}

FORCE_INLINE uint_fast8_t TBITSaturatedINC ( char arr[], size_t index ) {
    unsigned char theByte = (arr[TBITSLOT(index)]);
    uint_fast8_t theValue = (theByte >> TBITSHIFT(index)) & 3 ;
    if (theValue==3) {
        return 3;
    } else {
        ++theValue;
        theByte &= ~TBITMASK(index);
        theByte |= theValue << TBITSHIFT(index);
        arr[TBITSLOT(index)] = theByte;
        return theValue;
    }
}

FORCE_INLINE uint_fast8_t TBITSaturatedDEC ( char arr[], size_t index ) {
    unsigned char theByte = (arr[TBITSLOT(index)]);
    uint_fast8_t theValue = (theByte >> TBITSHIFT(index)) & 3 ;
    if (theValue==0) {
        return 0;
    } else {
        --theValue;
        theByte &= ~TBITMASK(index);
        theByte |= theValue << TBITSHIFT(index);
        arr[TBITSLOT(index)] = theByte;
        return theValue;
    }
}

