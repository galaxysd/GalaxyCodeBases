# learned from http://sourceforge.net/projects/gcmakefile/
# The options used in linking as well as in any direct use of ld.
LDLIB     = #-lpthread

LDFLAGS   = -Wl,-O1 -Wl,--sort-common -Wl,--enable-new-dtags -Wl,--hash-style=both $(LDLIB)
# The directories in which source files reside.
# If not specified, only the current directory will be serached.
SRCDIRS   =

## Implicit Section: change the following only when necessary.
##==========================================================================

# The source file types (headers excluded).
# .c indicates C source files, and others C++ ones.
SRCEXTS = .c# .C .cc .cpp .CPP .c++ .cxx .cp

OLevel= -O3
CFLAGS  = -pipe -march=core2 -mtune=generic -std=gnu99
WARNFLAGS =  -Wfloat-equal -Wall \
#-pedantic

ASMFLAGS = -S -fverbose-asm -g -masm=intel

ifeq ($(SRCDIRS),)
  SRCDIRS = .
endif

SOURCES = $(foreach d,$(SRCDIRS),$(wildcard $(addprefix $(d)/*,$(SRCEXTS))))

.PHONY: all clean

all: $(SOURCES:.c=)
#$(patsubst $.c,%,$(SOURCES))
	@echo [@gcc $(CFLAGS) $(OLevel) $(LDFLAGS)]
#@echo [$(SOURCES)] to [$(SOURCES:.c=)].

%: %.c
	@mkdir -p asm
	@echo  [$(OLevel)]:${<} -\> ${@}
	@gcc $(WARNFLAGS) $(CFLAGS) $(OLevel) $(LDFLAGS) -o $@ $<
	@gcc $(CFLAGS) $(OLevel) $(LDFLAGS) $(ASMFLAGS) -o ./asm/$@.s $<
	-as -alhnd -o /dev/null ./asm/$@.s > ./asm/$@.asm
#@gcc $(CFLAGS) $(OLevel) $(LDFLAGS) $(ASMFLAGS) -o - $< | as -alhnd > ./asm/$@.asm
	@rm ./asm/$@.s

duptc: override LDLIB += -lm

clean:
	-rm $(SOURCES:.c=)
	-rm -frv ./asm
