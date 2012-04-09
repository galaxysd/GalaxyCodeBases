CXX = g++
CFLAGS = -Wall -pipe -march=core2 -mtune=generic
LDFLAGS = -lz
OPT = -O3

MAKEARG = $(CFLAGS) $(OPT)

exefiles = filter_lowqual

all: $(exefiles)

filter_lowqual: filter_lowqual.cpp
	$(CXX) $(CFLAGS) $(OPT) $(LDFLAGS) -o filter_lowqual filter_lowqual.cpp gzstream.cpp

debug: override MAKEARG := $(CFLAGS) -O -D DEBUG -g
debug: all

test: override MAKEARG := $(CFLAGS) -O -D DEBUG -g
test: filter_adapter
	valgrind --leak-check=yes ./filter_lowqual t.bam

.PHONY : clean
clean:
	-rm $(exefiles)
