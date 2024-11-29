CC ?= gcc
CXX ?= g++
PKGCONFIG ?= pkg-config

ZLIBNAME ?= libisal
ZLIB_ID ?= 2

MKDIR_P = mkdir -p

PKGLIBS = $(ZLIBNAME)

PKGCPPFLAGS := $(shell $(PKGCONFIG) --cflags $(PKGLIBS))
PKGLDFLAGS := $(shell $(PKGCONFIG) --libs $(PKGLIBS))
PKGLIBPATH := $(shell $(PKGCONFIG) --variable=libdir $(PKGLIBS))
WARNINGFLAGS := -Wall -Wextra -Wpedantic -Werror=odr -Werror=strict-aliasing # -Wconversion

CFLAGS := -pipe -march=x86-64-v3 -mtune=generic -funroll-loops -flto -fPIE
# https://en.wikipedia.org/wiki/X86-64#Microarchitecture_levels
LDFLAGS := $(PKGLDFLAGS) -Wl,-pie -lc -lm -pthread
OPT := -Os
CFLAGS += -fstack-protector-strong -fcf-protection

ifeq ($(OS),Windows_NT) 
    detected_OS := Windows
else
    detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')
endif
OBJECTOOL := ldd
ifeq ($(detected_OS),Windows)
    CFLAGS += -D WIN32
endif
ifeq ($(detected_OS),Darwin)
	CC = /usr/local/opt/llvm/bin/clang
	CXX = /usr/local/opt/llvm/bin/clang++
	OBJECTOOL := otool -L
	LDFLAGS += -L/usr/local/opt/llvm/lib/c++ -L/usr/local/opt/llvm/lib -lunwind
	#PKGLDFLAGS := -Wl,-search_paths_first $(PKGLDFLAGS)
endif
ifeq ($(detected_OS),Linux)
	OBJECTOOL := ldd
	CFLAGS += -fstack-clash-protection
	#PKGLDFLAGS := -Wl,-Bstatic $(PKGLDFLAGS) -Wl,-Bdynamic -Wl,--as-needed
	LDFLAGS += -Wl,-z,defs -Wl,-z,now -Wl,-z,relro
	WARNINGFLAGS += -Werror=lto-type-mismatch
endif

ifeq ($(shell $(CC) -v 2>&1 | grep -c "clang version"),1)
	#clang
	LDFLAGS += -Wl,-rpath,$(PKGLIBPATH) -fuse-ld=lld
else
	#gcc
endif

SRCDIR := src1
SRCEXTS := .c
BUILDIR := build1
BUILT_PROGRAMS = salustsFstCoord1

SOURCES := $(foreach d,$(SRCDIR),$(wildcard $(addprefix $(d)/*,$(SRCEXTS))))
OBJFILES := $(SOURCES:$(SRCDIR)/%.c=$(BUILDIR)/%.o)

VPATH := $(SRCDIR) klib
vpath %.o $(BUILDIR)

CPPFLAGS := $(addprefix -I,$(VPATH)) $(PKGCPPFLAGS) -DZLIB_ID=$(ZLIB_ID) -D_FORTIFY_SOURCE=2
CSTANDARD := -std=gnu2x# -std=c23

all: $(BUILT_PROGRAMS)
	$(OBJECTOOL) $^

debug: override CPPFLAGS += -DDEBUG -g -fsanitize=address
debug: override OPT := -Og
debug: all

release: override CPPFLAGS += -DNDEBUG -DRELEASE
release: override OPT := -Ofast -finline-functions -funswitch-loops
release: all

z1: clean
	$(MAKE) ZLIBNAME=zlib-ng ZLIB_ID=1

z0: clean
	$(MAKE) ZLIBNAME=zlib ZLIB_ID=0

clean:
	@: # Clean objects and executables.
	rm -f $(BUILDIR)/*.o
	rm -f $(BUILT_PROGRAMS)

$(BUILDIR):
	mkdir -p $(BUILDIR)

$(BUILDIR)/%.o: $(SRCDIR)/%.c $(BUILDIR)
	$(CC) -c $(CPPFLAGS) $(WARNINGFLAGS) $(CFLAGS) $(OPT) $< -o $@

%: $(BUILDIR)/%.o	# sth. like `make uv_callback.test`
	$(CC) $(CFLAGS) $(WARNINGFLAGS) $(CPPFLAGS) $^ $(LDFLAGS) $(TARGET_ARCH) $(LOADLIBES) $(LDLIBS) -o $@

$(BUILT_PROGRAMS): $(OBJFILES)
	$(CC) $(CFLAGS) $(WARNINGFLAGS) $(CPPFLAGS) $^ $(LDFLAGS) $(TARGET_ARCH) $(LOADLIBES) $(LDLIBS) -o $@

.SECONDARY: $(OBJFILES)

#$(VERBOSE).SILENT:
echo:
	@echo $(LOCALPCPATH)
	@echo $(THISDIR) $(MAKEFILE_LIST)
	@echo PKGCPPFLAGS:[$(PKGCPPFLAGS)]
	@echo PKGLDFLAGS:[$(PKGLDFLAGS)]
	@echo SOURCES:[$(SOURCES)]
	@echo OBJFILES:[$(OBJFILES)]

test:
	UV_THREADPOOL_SIZE=1 MallocNanoZone=0 ./$(BUILT_PROGRAMS) n.fq xx
	UV_THREADPOOL_SIZE=2 valgrind -q --leak-check=yes ./$(BUILT_PROGRAMS) n.fq xx

fmt:
	find $(SRCDIR) -type f -exec clang-format -i {} \;

.PHONY: all clean debug release test echo z0 z1
