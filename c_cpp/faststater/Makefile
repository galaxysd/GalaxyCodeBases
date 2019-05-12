#CC = gcc
CFLAGS = -Wall -pipe -march=core2 -mtune=generic
LDFLAGS = -lz -lm
OPT = -O3

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CC=gcc-8
	CFLAGS += -I/usr/local/include
	LDFLAGS += -L/usr/local/lib -largp
endif

OBJDIR = ./objects/

MAKEARG = $(CFLAGS) $(OPT)

cobjects =
c99objects = ./objects/gFileIO.o
mainobjects = ./objects/fcount.o ./objects/gcount.o
objects = $(cobjects) $(c99objects)
exefiles = readsCounter homopolymerCounter

all: clean $(objects) $(exefiles)

readsCounter: $(objects) $(mainobjects)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o readsCounter $(objects) ./objects/fcount.o
homopolymerCounter: $(objects) $(mainobjects)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o homopolymerCounter $(objects) ./objects/gcount.o

$(cobjects): $(OBJDIR)%.o: %.c 
	$(CC) $(MAKEARG) -c $< -o $@

$(c99objects) $(mainobjects): $(OBJDIR)%.o: %.c
	$(CC) -std=gnu99 $(MAKEARG) -c $< -o $@

debug: override MAKEARG := $(CFLAGS) -O -D DEBUG -g
debug: all

.PHONY : clean
clean:
	-rm $(exefiles) $(mainobjects) $(objects)
	@mkdir -p $(OBJDIR)

