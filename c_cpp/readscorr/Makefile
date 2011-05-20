CC = gcc
CPP = g++
CFLAGS = -Wall -pipe -march=core2 -mtune=generic -pedantic\
#-ftree-vectorize -ftree-vectorizer-verbose=5
LDFLAGS = -lz -lm -Wl,-O1 -Wl,--as-needed -Wl,--sort-common -Wl,--enable-new-dtags -Wl,--hash-style=both
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
 
cobjects =
customobjects =
c99objects = main.o MurmurHash3.o gFileIO.o sdleft.o 2bitseq.o
objects = $(cobjects) $(c99objects)

all: clean $(objects) readscorr

readscorr : $(objects) $(customobjects)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o readscorr $(objects)

$(cobjects): %.o: %.c 
	$(CC) $(MAKEARG) -c $< -o $@ > $@.asm
$(cppobjects): %.o: %.cpp
	$(CPP) $(MAKEARG) -c $< -o $@ > $@.asm
$(c99objects): %.o: %.c
	$(CC) -std=gnu99 $(MAKEARG) -c $< -o $@ > $@.asm

debug: override MAKEARG := $(CFLAGS) $(DEF) -O -D DEBUG -g $(ASMO)
debug: clean $(objects) 
	$(CC) $(CFLAGS) $(LDFLAGS) -g -o readscorr $(objects)

.PHONY : clean
clean :
	-rm readscorr $(objects) $(customobjects)
