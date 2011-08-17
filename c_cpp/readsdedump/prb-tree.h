#ifndef PRB_TREE_H
#define PRB_TREE_H 1

#include "prb.h"

int compare_ints (const void *pa, const void *pb, void *param);
int compare_strings (const void *pa, const void *pb, void *param);
int compare_fixed_strings (const void *pa, const void *pb, void *param);
void print_whole_tree (const struct prb_table *tree, const char *title);

#endif /* prb-tree.h */
