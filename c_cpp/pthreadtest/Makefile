# The options used in linking as well as in any direct use of ld.
LDLIB     = -lpthread

LDFLAGS   = -Wl,-O1 -Wl,--sort-common -Wl,--enable-new-dtags -Wl,--hash-style=both $(LDLIB)
# The directories in which source files reside.
# If not specified, only the current directory will be serached.
SRCDIRS   =

## Implicit Section: change the following only when necessary.
##==========================================================================

# The source file types (headers excluded).
# .c indicates C source files, and others C++ ones.
SRCEXTS = .c .C .cc .cpp .CPP .c++ .cxx .cp

CFLAGS  = -pipe -march=core2 -mtune=generic -Wfloat-equal \
#-Wall -pedantic

ifeq ($(SRCDIRS),)
  SRCDIRS = .
endif
SOURCES = $(foreach d,$(SRCDIRS),$(wildcard $(addprefix $(d)/*,$(SRCEXTS))))
TARGET  = $(addprefix bin/,$(basename $(SOURCES)))

all: $(TARGET)
	ls -l $(TARGET)

bin/%: %.c dir
	gcc $(CFLAGS) $(LDFLAGS) -o $@ $<

dir:
	@-mkdir bin

.PHONY : clean
clean:
	-rm $(TARGET)
	-rmdir bin
