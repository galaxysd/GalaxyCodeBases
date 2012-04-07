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
#include <stdlib.h>
#include <string.h>
#include <bam/sam.h>

struct myData {
    int32_t n_targets;
    char **target_name;
    uint32_t *target_len;
    uint8_t **ChrDat;
} Data;

void do_stat(bam1_t *balignd, const uint16_t overlap, const struct myData *Data);

int do_contig(const uint8_t mindep, const struct myData *, const char * const outfile);

char *strlinker(const char * const main, const char * const suffix);







#endif /* CHR_TABLE_H */
