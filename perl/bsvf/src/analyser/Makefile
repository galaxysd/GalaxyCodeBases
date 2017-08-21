CC = gcc
CXX = g++
CFLAGS = -std=gnu11 -Wall -pedantic -pipe -march=core2 -mtune=generic -I/usr/local/include
CXXFLAGS = -std=gnu++14 -Wall -pedantic -pipe -march=core2 -mtune=generic -I/usr/local/include
LDFLAGS = -L/usr/local/lib -lhts #-lpthread -pthread
OPT = -O3

exefiles = bsanalyser

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CC = gcc-7
	CXX = g++-7
	LDFLAGS += -largp -D_THREAD_SAFE #-static-libgcc -static-libstdc++
endif
ifeq ($(UNAME_S),Linux)
	LDFLAGS += -static -lbz2 -lm -lpthread -lz -llzma -pthread
endif

OBJDIR = ./objects/

#MAKEARG = $(CFLAGS) $(OPT)

cobjects = ./objects/ini.o ./objects/samio.o ./objects/samlib.o ./objects/analyse.o
cxxobjects = 
c99objects = 
mainobjects = ./objects/main.o
objects = $(cobjects) $(c99objects)

all: $(objects) $(exefiles)

$(exefiles): $(objects) $(cxxobjects) $(mainobjects)
	#$(CXX) $(CXXFLAGS) $(OPT) $(LDFLAGS) -o $(exefiles) $(mainobjects) $(objects)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o $(exefiles) $(mainobjects) $(objects)

$(mainobjects): $(OBJDIR)%.o: %.c tmpdir
	#$(CXX) $(CXXFLAGS) $(OPT) -c $< -o $@
	$(CC) $(CFLAGS) $(OPT) -c $< -o $@

$(cobjects): $(OBJDIR)%.o: %.c tmpdir
	$(CC) $(CFLAGS) $(OPT) -c $< -o $@

$(cxxobjects): $(OBJDIR)%.o: %.c tmpdir
	$(CXX) $(CXXFLAGS) $(OPT) -c $< -o $@

$(c99objects): $(OBJDIR)%.o: %.c tmpdir
	$(CC) $(CFLAGS) $(OPT) -c $< -o $@

debug: override MAKEARG := $(CFLAGS) -O -D DEBUGa -g
debug: all

tmpdir:
	-mkdir ./objects

test: override MAKEARG := $(CFLAGS) -O -D DEBUGa -g
test: $(exefiles)
	valgrind --leak-check=yes ./$(exefiles) -p g ToGrep.ini
run: $(exefiles)
	valgrind --leak-check=yes ./$(exefiles) -p g ToGrep.ini

pp:
	$(CC) -E main.c | indent > mm.c
	-mate mm.c 
sign: $(exefiles)
	codesign -s "Mac Developer" -v $(exefiles)
	codesign -d -v $(exefiles)

.PHONY : clean
clean:
	-rm $(exefiles) $(mainobjects) $(objects)
	@mkdir -p $(OBJDIR)

