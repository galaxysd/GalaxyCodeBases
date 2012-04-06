CC = gcc
CFLAGS = -Wall -pipe -march=core2 -mtune=generic
LDFLAGS = -lz -lm -lbam
OPT = -O3

OBJDIR = ./objects/

MAKEARG = $(CFLAGS) $(OPT)

cobjects =
c99objects = ./objects/chrtable.o
mainobjects = ./objects/constator.o
objects = $(cobjects) $(c99objects)
exefiles = constator

all: $(objects) $(exefiles)

constator: $(objects) $(mainobjects)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o constator $(objects) ./objects/constator.o

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

