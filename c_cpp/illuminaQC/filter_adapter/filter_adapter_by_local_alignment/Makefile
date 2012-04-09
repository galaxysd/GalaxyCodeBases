CXX = g++
CFLAGS = -Wall -pipe -march=core2 -mtune=generic
LDFLAGS = -lz
OPT = -O3

MAKEARG = $(CFLAGS) $(OPT)

exefiles = filter_adapter

all: $(exefiles)

filter_adapter:
	$(CXX) $(CFLAGS) $(OPT) $(LDFLAGS) -o filter_adapter filter_adapter.cpp gzstream.cpp

debug: override MAKEARG := $(CFLAGS) -O -D DEBUG -g
debug: all

test: override MAKEARG := $(CFLAGS) -O -D DEBUG -g
test: filter_adapter
	valgrind --leak-check=yes ./filter_adapter t.bam

.PHONY : clean
clean:
	-rm $(exefiles)
