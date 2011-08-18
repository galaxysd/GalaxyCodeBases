#ifndef PRB_TREE_H
#define PRB_TREE_H 1

#include <stdint.h>
#include "prb.h"

#define IDARRAY_INIT (6)
#define IDARRAY_INC  (2)

typedef struct __prb_kmer_pos_data {
    char* pKmer;    // without tailing '\0', for KmerLength == *param
    uint16_t pos;
    uint32_t IDcnt;
    uint32_t* pIDary;
} kmerpos_t;

int compare_ints (const void *pa, const void *pb, void *param);
int compare_strings (const void *pa, const void *pb, void *param);
int compare_fixed_strings (const void *pa, const void *pb, void *param);
void print_whole_tree (const struct prb_table *tree, const char *title);

#endif /* prb-tree.h */
