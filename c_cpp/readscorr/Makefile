CC = gcc
CPP = g++
CFLAGS = -Wall
DEBUG = -g
LDFLAGS = -lm
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
	$(CPP) $(CFLAGS) $(OPT) $(LDFLAGS) -o readscorr $(objects)

$(cobjects): %.o: %.c 
	$(CC) $(MAKEARG) -c $< -o $@ > $@.asm
$(cppobjects): %.o: %.cpp
	$(CPP) $(MAKEARG) -c $< -o $@ > $@.asm
$(c99objects): %.o: %.c
	$(CC) $(MAKEARG) -std=c99 -c $< -o $@ > $@.asm

debug: override MAKEARG := $(CFLAGS) $(DEF) -D DEBUG -g $(ASMO)
debug: clean $(objects) 
	$(CPP) $(CFLAGS) $(LDFLAGS) $(DEBUG) -o readscorr $(objects)

.PHONY : clean
clean :
	-rm readscorr $(objects) $(cppobjects)
