# Copyright 2007 Google Inc.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

CXXFLAGS = -O3
CPPFLAGS = -DNDEBUG

all: dcsbwt dcsunbwt

dcsbwt: dcsbwt.cc difference_cover.o bw_transform.o bwt_compress.o \
	rl_compress.o \
	stream.o stream_compressor.o inverse_bwt.o bbb_compress.o
dcsunbwt: dcsunbwt.cc bwt_compress.o inverse_bwt.o stream.o \
	difference_cover.o \
	stream_compressor.o  bw_transform.o rl_compress.o bbb_compress.o

clean:
	rm *.o *~ dcsbwt dcsunbwt
