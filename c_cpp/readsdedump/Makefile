CC = gcc
#CPP = g++
CFLAGS = -Wall -pipe -march=core2 -mtune=generic -pedantic -Wfloat-equal\
# -Wno-packed-bitfield-compat\
# -Wpadded -Wextra
# -DMKL_ILP64 -m64 -I$(MKLROOT)/include \
# -ftree-vectorize -ftree-vectorizer-verbose=5
LDFLAGS = -lz -lm -Wl,-O1 -Wl,--as-needed -Wl,--sort-common -Wl,--enable-new-dtags -Wl,--hash-style=both \
# $(MKLROOT)/lib/em64t/libmkl_solver_ilp64_sequential.a -Wl,--start-group $(MKLROOT)/lib/em64t/libmkl_intel_ilp64.a $(MKLROOT)/lib/em64t/libmkl_sequential.a $(MKLROOT)/lib/em64t/libmkl_core.a -Wl,--end-group -lpthread
OPT = -O3
ASMO = -Wa,-a,-ad

OBJDIR = ./objects/

DEF = #-D NEW

MAKEARG = $(CFLAGS) $(DEF) $(OPT) $(ASMO)

# OpenMP 
ifdef OPENMP
override CFLAGS := $(CFLAGS) -fopenmp
override DEF := $(DEF) -D OPENMP
endif

# pthread
ifdef PTHREAD
override LDFLAGS := $(LDFLAGS) -lpthread
override DEF := $(DEF) -D PTHREAD
endif

cobjects =
customobjects = $(addprefix $(OBJDIR),)# sdleftTFuint16.o sdleftTFuint32.o sdleftTFuint64.o sdleftTFpublic.o )
c99objects = $(addprefix $(OBJDIR), prb.o)
noshareobjects = $(addprefix $(OBJDIR), prbtest.o main.o simple.o prb-tree.o prb-test.o)
objects = $(cobjects) $(c99objects) $(customobjects)
exefiles = prb-test readsdedump simplededump

all: clean $(objects) $(exefiles) stat_dup
	@echo [$(CC) $(CFLAGS) $(OPT) $(LDFLAGS)][$(MAKEARG)].

readsdedump: $(objects) $(noshareobjects)
	@$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o readsdedump $(objects) $(addprefix $(OBJDIR), prb-tree.o main.o)
simplededump: $(objects) $(noshareobjects)
	@$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o simplededump $(objects) $(addprefix $(OBJDIR), simple.o)
prb-test: $(objects) $(noshareobjects)
	@$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o prb-test $(objects) $(addprefix $(OBJDIR), prb-test.o prbtest.o)
stat_dup: 
	g++ $(CFLAGS) $(OPT) $(LDFLAGS) -o stat_dup -I. duplication.cpp gzstream.C -lz

$(cobjects): $(OBJDIR)%.o: %.c 
	@$(CC) $(MAKEARG) -c $< -o $@ > $(@:.o=.asm)
#$(cppobjects): %.o: %.cpp
#	$(CPP) $(MAKEARG) -c $< -o $@ > $@.asm
$(c99objects) $(noshareobjects): $(OBJDIR)%.o: %.c
	@$(CC) -std=gnu99 $(MAKEARG) -c $< -o $@ > $(@:.o=.asm)

debug: override MAKEARG := $(CFLAGS) $(DEF) -O -D DEBUG -g $(ASMO)
debug: all

alpha: override MAKEARG := $(CFLAGS) $(DEF) $(OPT) $(ASMO) -D TEST
alpha: all

alphadbg: override MAKEARG := $(CFLAGS) $(DEF) $(OPT) $(ASMO) -D TEST -D DEBUG -g
alphadbg: all

alphapf: override MAKEARG := $(CFLAGS) $(DEF) $(OPT) $(ASMO) -D TEST -pg -g
alphapf: override LDFLAGS += -pg
alphapf: gprofmain

gprofmain: clean readscorr
	./readscorr -o gprof test.lst 2> gprof.log
	gprof -l readscorr > gprof.report

test: readscorr
	valgrind --leak-check=yes ./readscorr test.lst -k 17

.PHONY : clean
clean:
	-rm $(exefiles) $(noshareobjects) $(objects)
	@mkdir -p $(OBJDIR)

