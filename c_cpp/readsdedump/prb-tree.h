#ifndef PRB_TREE_H
#define PRB_TREE_H 1

#include "prb.h"

char *pgm_name; // used for die.perl, thus can be overwritten.

int compare_ints (const void *pa, const void *pb, void *param);
int compare_strings (const void *pa, const void *pb, void *param);
int compare_fixed_strings (const void *pa, const void *pb, void *param);
void print_whole_tree (const struct prb_table *tree, const char *title);
void pgm_fail (const char *message, ...);

#endif /* prb-tree.h */
