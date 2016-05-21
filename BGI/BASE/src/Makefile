#
#    Makefile
#
#    BASE
#
#    Copyright (C) 2015, The University of Hong Kong.
#
#    This program is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License
#    as published by the Free Software Foundation; either version 3
#    of the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  03110-1301, USA.
#
####################################################################
#
#   Modification History
#
#   Date   : 3th Oct. 2014
#   Author : Binghang LIU
#
####################################################################

CC = g++ 
CFLAGS = -O3 -funroll-loops -march=nocona -maccumulate-outgoing-args -Wno-unused-result -Wno-write-strings -lm -lz -lrt -g -fomit-frame-pointer -pthread

LIBS = -lm -lz -lrt -fomit-frame-pointer
THREADLIBS = -lpthread

BWTLIB = ./2bwt
BWTOBJLIBS = $(BWTLIB)/2BWT-Interface.o $(BWTLIB)/BWTConstruct.o $(BWTLIB)/BWT.o $(BWTLIB)/dictionary.o $(BWTLIB)/DNACount.o $(BWTLIB)/HSP.o $(BWTLIB)/HSPstatistic.o $(BWTLIB)/iniparser.o $(BWTLIB)/inistrlib.o $(BWTLIB)/karlin.o $(BWTLIB)/MemManager.o $(BWTLIB)/MiscUtilities.o $(BWTLIB)/overlap-api.o $(BWTLIB)/QSufSort.o $(BWTLIB)/r250.o $(BWTLIB)/Socket.o $(BWTLIB)/TextConverter.o $(BWTLIB)/Timing.o

all:	base 

clean:
	rm base *.o $(BWTLIB)/*.o $(BWTLIB)/2bwt_builder_cpu $(BWTLIB)/AsiicBWT2BWT 

$(BWTLIB) : force_look
	cd $(BWTLIB); $(MAKE)


main.o:	main.c Makefile
contigs.o: contigs.c Makefile
extend.o: extend.c Makefile
bwt_api.o: bwt_api.c Makefile

base: main.o contigs.o extend.o bwt_api.o $(BWTOBJLIBS) Makefile
	$(CC) $(CFLAGS) main.o extend.o contigs.o bwt_api.o $(BWTOBJLIBS) -o base

