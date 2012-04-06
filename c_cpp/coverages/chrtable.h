#ifndef CHR_TABLE_H
#define CHR_TABLE_H
/*
#include "uthash/uthash.h"

struct ChrData_hash_struct {
    char *ChrID;
    uint32_t ChrLen;
    uint8_t *Data;
    UT_hash_handle hh;
};

struct ChrData_hash_struct *ChrData;

void add_chr(char *id, uint32_t len);
*/
#include <bam/sam.h>

void do_stat(bam1_t *balignd, uint16_t overlap, uint8_t **ChrDat);











#endif /* CHR_TABLE_H */
