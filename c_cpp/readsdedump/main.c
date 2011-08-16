#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include "prb.h"
#include "prb-tree.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

int
main (int argc, char *argv[])
{
  int success=1;                  /* Everything okay so far? */

  /* Initialize |pgm_name|, using |argv[0]| if sensible. */
  pgm_name = argv[0] != NULL && argv[0][0] != '\0' ? argv[0] : "readsdedump";

  //pgm_fail("%s","test");
  
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

