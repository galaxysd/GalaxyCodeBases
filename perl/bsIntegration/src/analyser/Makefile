CC ?= gcc
CFLAGS = -std=gnu11 -Wall -pipe -march=core2 -mtune=generic
LDFLAGS = -lz -lm -lhts -lbam
OPT = -O3

exefiles = bsanalyser

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	LDFLAGS += -largp
endif

OBJDIR = ./objects/

MAKEARG = $(CFLAGS) $(OPT)

cobjects = ./objects/ini.o
c99objects = 
mainobjects = ./objects/main.o
objects = $(cobjects) $(c99objects)

all: $(objects) $(exefiles)

$(exefiles): $(objects) $(mainobjects)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o $(exefiles) $(mainobjects) $(objects)

$(cobjects) $(mainobjects): $(OBJDIR)%.o: %.c tmpdir
	$(CC) $(MAKEARG) -c $< -o $@

$(c99objects): $(OBJDIR)%.o: %.c tmpdir
	$(CC) $(MAKEARG) -c $< -o $@

debug: override MAKEARG := $(CFLAGS) -O -D DEBUGa -g
debug: all

tmpdir:
	-mkdir ./objects

test: override MAKEARG := $(CFLAGS) -O -D DEBUGa -g
test: $(exefiles)
	valgrind --leak-check=yes ./$(exefiles) -p g ToGrep.ini

pp:
	$(CC) -E main.c | indent > mm.c
	-mate mm.c 

.PHONY : clean
clean:
	-rm $(exefiles) $(mainobjects) $(objects)
	@mkdir -p $(OBJDIR)

