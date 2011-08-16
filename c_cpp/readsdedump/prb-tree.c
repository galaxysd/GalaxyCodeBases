/* Produced by texiweb from libavl.w. */

/* libavl - library for manipulation of binary trees.
   Copyright (C) 1998-2002, 2004 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
   02111-1307, USA.

   The author may be contacted at <blp@gnu.org> on the Internet, or
   write to Ben Pfaff, Stanford University, Computer Science Dept., 353
   Serra Mall, Stanford CA 94305, USA.
*/

#include <assert.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prb-tree.h"

/* Utility functions. */

/* http://www.stanford.edu/~blp/avl/libavl.html/Answers-for-Chapter-2.html#2-3%232
   Comparison function for pointers to |int|s.
   |param| is not used. */
int 
compare_ints (const void *pa, const void *pb, void *param) 
{
  const int *a = pa;
  const int *b = pb;

  return (*a > *b) - (*a < *b);
}
/* Comparison function for strings. 
   param is not used. */
int 
compare_strings (const void *pa, const void *pb, void *param) 
{
  return strcmp (pa, pb);
}
int 
compare_fixed_strings (const void *pa, const void *pb, void *param) 
{
  return memcmp (pa, pb, *(size_t *) param);
}

/* Prints |message| on |stderr|, which is formatted as for |printf()|,
   and terminates the program unsuccessfully. */
void
pgm_fail (const char *message, ...)
{
  va_list args;

  fprintf (stderr, "%s: ", pgm_name);

  va_start (args, message);
  vfprintf (stderr, message, args);
  va_end (args);

  fputs("\n", stderr);

  exit (EXIT_FAILURE);
}
/*
   Allocates and returns a pointer to |size| bytes of memory.
   Aborts if allocation fails.
static void *
xmalloc (size_t size)
{
  void *block = malloc (size);
  if (block == NULL && size != 0)
    pgm_name="prb-tree.c";
    pgm_fail ("out of memory");
  return block;
}
*/

/* Prints the structure of |node|,
   which is |level| levels from the top of the tree. */
static void
print_tree_structure (const struct prb_node *node, int level)
{
  /* You can set the maximum level as high as you like.
     Most of the time, you'll want to debug code using small trees,
     so that a large |level| indicates a ``loop'', which is a bug. */
  if (level > 16)
    {
      printf ("[...]");
      return;
    }

  if (node == NULL)
    return;

  printf ("%d", *(int *) node->prb_data);
  if (node->prb_link[0] != NULL || node->prb_link[1] != NULL)
    {
      putchar ('(');

      print_tree_structure (node->prb_link[0], level + 1);
      if (node->prb_link[1] != NULL)
        {
          putchar (',');
          print_tree_structure (node->prb_link[1], level + 1);
        }

      putchar (')');
    }
}

/* Prints the entire structure of |tree| with the given |title|. */
void
print_whole_tree (const struct prb_table *tree, const char *title)
{
  printf ("%s: ", title);
  print_tree_structure (tree->prb_root, 0);
  putchar ('\n');
}

