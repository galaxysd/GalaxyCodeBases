/*

   TypeNLimit.h		Miscellaneous Constants

   Copyright (C) 2004, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/
   
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>

typedef long long int int64;
typedef unsigned int word;
typedef unsigned long long int uint64;
typedef unsigned char boolean;

struct xtimer_t {
  struct timeval tv1, tv2;
  int64 time_elapsed;
  void reset() { time_elapsed = 0; }
  void start() { gettimeofday(&tv1, NULL); }
  void stop() {
    gettimeofday(&tv2, NULL);
    time_elapsed += (int64)(tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
  }
  double elapsed() { return time_elapsed / 1000000.0; }
};
