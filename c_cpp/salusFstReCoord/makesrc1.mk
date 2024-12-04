CC ?= gcc
CXX ?= g++
PKGCONFIG ?= pkg-config
ZLIBNAME ?= libisal

MKDIR_P = mkdir -p
ZCAT = zcat

PKGLIBS = $(ZLIBNAME)

PKGCPPFLAGS := $(shell $(PKGCONFIG) --cflags $(PKGLIBS))
PKGLDFLAGS := $(shell $(PKGCONFIG) --libs $(PKGLIBS))
PKGLIBPATH := $(shell $(PKGCONFIG) --variable=libdir $(PKGLIBS))
WARNINGFLAGS := -Wall -Wextra -Wpedantic -Werror=odr -Werror=strict-aliasing # -Wconversion

CFLAGS := -pipe -march=x86-64-v3 -mtune=generic -funroll-loops -flto -fPIE # -pthread
# https://en.wikipedia.org/wiki/X86-64#Microarchitecture_levels
LDFLAGS := $(PKGLDFLAGS) -Wl,-pie -lc -lm # -pthread
OPT := -Os
CFLAGS += -fstack-protector-strong -fcf-protection

lc = $(subst A,a,$(subst B,b,$(subst C,c,$(subst D,d,$(subst E,e,$(subst F,f,$(subst G,g,$(subst H,h,$(subst I,i,$(subst J,j,$(subst K,k,$(subst L,l,$(subst M,m,$(subst N,n,$(subst O,o,$(subst P,p,$(subst Q,q,$(subst R,r,$(subst S,s,$(subst T,t,$(subst U,u,$(subst V,v,$(subst W,w,$(subst X,x,$(subst Y,y,$(subst Z,z,$1))))))))))))))))))))))))))
uc = $(subst a,A,$(subst b,B,$(subst c,C,$(subst d,D,$(subst e,E,$(subst f,F,$(subst g,G,$(subst h,H,$(subst i,I,$(subst j,J,$(subst k,K,$(subst l,L,$(subst m,M,$(subst n,N,$(subst o,O,$(subst p,P,$(subst q,Q,$(subst r,R,$(subst s,S,$(subst t,T,$(subst u,U,$(subst v,V,$(subst w,W,$(subst x,X,$(subst y,Y,$(subst z,Z,$1))))))))))))))))))))))))))

ZLIBID := $(call uc,$(subst -,,$(ZLIBNAME)))

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
	ZCAT = gzcat
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

CPPFLAGS := $(addprefix -I,$(VPATH)) $(PKGCPPFLAGS) -DUSE_${ZLIBID} -D_FORTIFY_SOURCE=2
CSTANDARD := -std=gnu2x# -std=c23

all: $(BUILT_PROGRAMS)
	$(OBJECTOOL) $^

debug: override CPPFLAGS += -DDEBUG -g -fsanitize=address
debug: override OPT := -Og
debug: clean all

release: override CPPFLAGS += -DNDEBUG -DRELEASE
release: override OPT := -Ofast -finline-functions -funswitch-loops
release: clean all

z1: clean
	$(MAKE) release ZLIBNAME=zlib-ng

z0: clean
	$(MAKE) release ZLIBNAME=zlib

clean:
	@: # Clean objects and executables.
	rm -f $(BUILDIR)/*.o
	rm -f $(BUILT_PROGRAMS)

$(BUILDIR):
	$(MKDIR_P) $(BUILDIR)

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
	#UV_THREADPOOL_SIZE=1 MallocNanoZone=0 ./$(BUILT_PROGRAMS) n.fq >t.fq
	UV_THREADPOOL_SIZE=2 valgrind -q --leak-check=yes ./$(BUILT_PROGRAMS) n.fq >t1.fq

fmt:
	find $(SRCDIR) -type f -exec clang-format -i {} \;

cmp:
	$(MKDIR_P) tmp
	python3 SpatialOmicsCoord.py n.fq tmp 1.25
	$(ZCAT) tmp/newCoord_n.gz >t2.fq
	./$(BUILT_PROGRAMS) n.fq 1.25 >t1.fq
	diff -Nau t1.fq t2.fq
	shasum t1.fq t2.fq

pkg:
	tar -chz -f salusFstReCoord.tgz $(VPATH) SpatialOmicsCoord.py n.fq Makefile

.PHONY: all clean debug release test echo z0 z1 cmp pkg
