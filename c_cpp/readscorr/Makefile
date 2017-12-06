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
c99objects = $(addprefix $(OBJDIR), MurmurHash3.o gFileIO.o 2bitseq.o cfgparser.o gtools.o bihash.o)
mainobjects = $(addprefix $(OBJDIR), main.o tester.o fcount.o gcount.o )
objects = $(cobjects) $(c99objects) $(customobjects)
exefiles = readscorr tester fcount gcount

all: clean $(objects) $(exefiles)
	@echo [$(CC) $(CFLAGS) $(OPT) $(LDFLAGS)][$(MAKEARG)].

readscorr: $(objects) $(mainobjects)
	@$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o readscorr $(objects) $(addprefix $(OBJDIR), main.o)
tester: $(objects) $(mainobjects)
	@$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o tester $(objects) $(addprefix $(OBJDIR), tester.o)
fcount: $(objects) $(mainobjects)
	@$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o fcount $(objects) $(addprefix $(OBJDIR), fcount.o)
gcount: $(objects) $(mainobjects)
	@$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o gcount $(objects) $(addprefix $(OBJDIR), gcount.o)

#$(OBJDIR)sdleftTFuint16.o:
#	$(CC) -std=gnu99 $(MAKEARG) -D USEUINT16 -c sdleftTF.c -o $@ > $(@:.o=.asm)
#$(OBJDIR)sdleftTFuint32.o:
#	$(CC) -std=gnu99 $(MAKEARG) -D USEUINT32 -c sdleftTF.c -o $@ > $(@:.o=.asm)
#$(OBJDIR)sdleftTFuint64.o:
#	$(CC) -std=gnu99 $(MAKEARG) -D USEUINT64 -c sdleftTF.c -o $@ > $(@:.o=.asm)
#$(OBJDIR)sdleftTFpublic.o:
#	$(CC) -std=gnu99 $(MAKEARG) -D PUBLIC -c sdleftTF.c -o $@ > $(@:.o=.asm)

$(cobjects): $(OBJDIR)%.o: %.c 
	@$(CC) $(MAKEARG) -c $< -o $@ > $(@:.o=.asm)
#$(cppobjects): %.o: %.cpp
#	$(CPP) $(MAKEARG) -c $< -o $@ > $@.asm
$(c99objects) $(mainobjects): $(OBJDIR)%.o: %.c
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
	-rm $(exefiles) $(mainobjects) $(objects)
	@mkdir -p $(OBJDIR)

