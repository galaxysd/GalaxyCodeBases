#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "prb.h"
#include "prb-tree.h"
//#include <stdarg.h>
//#include <stdio.h>

char * pinsert[10] = {"ABCDE","FGHIJ","ABCDES","FGHIJ","ABSOPQ","1BSOPQ","2BSOPQ",NULL};

//KmerPosID_t data[10];

int
main (int argc, char *argv[])
{
  int seedlen=5;
  int success=1;                  /* Everything okay so far? */

  struct prb_table *tree;
  tree = prb_create (compare_fixed_strings, &seedlen, NULL);
  if (tree == NULL)
    {
      printf ("  Out of memory creating tree.\n");
      return 1;
    }

  int verbosity=3;
  for (int i = 0; pinsert[i]; i++)
    {
      if (verbosity >= 2)
        printf ("  Inserting %d:[%s]...\n", i, pinsert[i]);

      /* Add the |i|th element to the tree. */
      {
        void **p = prb_probe (tree, pinsert[i]);
        if (p == NULL)
          {
            if (verbosity >= 0)
              printf ("    Out of memory in insertion.\n");
            prb_destroy (tree, NULL);
            return 1;
          }
        if (*p != pinsert[i])
          printf ("    Duplicate item in tree!\n");
      }

      if (verbosity >= 3)
        print_whole_tree (tree, "    Afterward");

    }


  //pgm_fail("%s","test");
  prb_destroy (tree, NULL);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

