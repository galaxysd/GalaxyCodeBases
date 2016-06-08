CC = gcc
CFLAGS = -Wall -pipe -march=core2 -mtune=generic
LDFLAGS = -lz -lm -lhts -lbam
OPT = -O3

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	LDFLAGS += -largp
endif

OBJDIR = ./objects/

MAKEARG = $(CFLAGS) $(OPT)

cobjects =
c99objects = 
mainobjects = ./objects/main.o
objects = $(cobjects) $(c99objects)
exefiles = constator

all: $(objects) $(exefiles)

constator: $(objects) $(mainobjects)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o bsanalyser $(mainobjects) $(objects)

$(cobjects): $(OBJDIR)%.o: %.c tmpdir
	$(CC) $(MAKEARG) -c $< -o $@

$(c99objects) $(mainobjects): $(OBJDIR)%.o: %.c tmpdir
	$(CC) $(MAKEARG) -c $< -o $@

debug: override MAKEARG := $(CFLAGS) -O -D DEBUGa -g
debug: all

tmpdir:
	-mkdir ./objects

test: override MAKEARG := $(CFLAGS) -O -D DEBUGa -g
test: constator
	valgrind --leak-check=yes ./bsanalyser t.bam

.PHONY : clean
clean:
	-rm $(exefiles) $(mainobjects) $(objects)
	@mkdir -p $(OBJDIR)

