CC = gcc
CPP = g++
CFLAGS = -Wall -pipe -march=core2 -mtune=generic
DEBUG = -g
LDFLAGS = -lz -Wl,-O1 -Wl,--as-needed -Wl,--sort-common -Wl,--enable-new-dtags -Wl,--hash-style=both
OPT = -O3
ASMO = -Wa,-a,-ad

MAXKMERLENGTH=31
DEF = -D MAXKMERLENGTH=$(MAXKMERLENGTH) 

MAKEARG = $(CFLAGS) $(DEF) $(OPT) $(ASMO)

# OpenMP 
ifdef OPENMP
override CFLAGS := $(CFLAGS) -fopenmp
override DEF := $(DEF) -D OPENMP
endif
 
cobjects = main.o
cppobjects = 
c99objects = MurmurHash3.o
objects = $(cobjects) $(cppobjects) $(c99objects)

all: clean $(objects) readscorr

readscorr : $(objects)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o readscorr $(objects)

$(cobjects): %.o: %.c 
	$(CC) $(MAKEARG) -c $< -o $@ > $@.asm
$(cppobjects): %.o: %.cpp
	$(CPP) $(MAKEARG) -c $< -o $@ > $@.asm
$(c99objects): %.o: %.c
	$(CC) -std=c99 $(MAKEARG) -c $< -o $@ > $@.asm

debug: override MAKEARG := $(CFLAGS) $(DEF) -O -D DEBUG -g $(ASMO)
debug: clean $(objects) 
	$(CC) $(CFLAGS) $(LDFLAGS) $(DEBUG) -o readscorr $(objects)

.PHONY : clean
clean :
	-rm readscorr $(objects) $(cppobjects)
