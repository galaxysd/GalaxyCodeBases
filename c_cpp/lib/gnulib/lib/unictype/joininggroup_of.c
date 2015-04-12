/* Arabic joining group of Unicode characters.
   Copyright (C) 2011-2015 Free Software Foundation, Inc.
   Written by Bruno Haible <bruno@clisp.org>, 2011.

   This program is free software: you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

#include <config.h>

/* Specification.  */
#include "unictype.h"

/* Define u_joining_group table.  */
#include "joininggroup_of.h"

int
uc_joining_group (ucs4_t uc)
{
  unsigned int index1 = uc >> joining_group_header_0;
  if (index1 < joining_group_header_1)
    {
      int lookup1 = u_joining_group.level1[index1];
      if (lookup1 >= 0)
        {
          unsigned int index2 = (uc >> joining_group_header_2) & joining_group_header_3;
          int lookup2 = u_joining_group.level2[lookup1 + index2];
          if (lookup2 >= 0)
            {
              unsigned int index3 = ((uc & joining_group_header_4) + lookup2) * 7;
              /* level3 contains 7-bit values, packed into 16-bit words.  */
              unsigned int lookup3 =
                ((u_joining_group.level3[index3>>4]
                  | ((unsigned int) u_joining_group.level3[(index3>>4)+1] << 16))
                 >> (index3 % 16))
                & 0x7f;

              return lookup3;
            }
        }
    }
  return UC_JOINING_GROUP_NONE;
}
