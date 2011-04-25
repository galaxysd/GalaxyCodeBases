CC = gcc
CPP = g++
CFLAGS = -Wall
DEBUG = -g
LDFLAGS = -lm
OPT = -O3
ASMO = -Wa,-a,-ad

MAXKMERLENGTH=31
DEF = -D MAXKMERLENGTH=$(MAXKMERLENGTH)

# OpenMP 
ifdef OPENMP
override CFLAGS := $(CFLAGS) -fopenmp
override DEF := $(DEF) -D OPENMP
endif
 
cobjects = main.o
cppobjects = MurmurHash3.o
objects = $(cobjects) $(cppobjects)

all: clean $(objects) readscorr

readscorr : $(objects)
	$(CPP) $(CFLAGS) $(OPT) $(LDFLAGS) -o readscorr $(objects)

$(cobjects): %.o: %.c 
	$(CC) $(CFLAGS) $(OPT) $(DEF) -c $< -o $@ > $@.asm
$(cppobjects): %.o: %.cpp
	$(CPP) $(CFLAGS) $(OPT) $(DEF) -c $< -o $@ > $@.asm

debug: override DEF := $(DEF) -D DEBUG -g $(ASMO)
debug: clean $(objects) 
	$(CPP) $(CFLAGS) $(LDFLAGS) $(DEBUG) -o readscorr $(objects)

.PHONY : clean
clean :
	-rm readscorr $(objects) $(cppobjects)
