LDLIB     = -lz

LDFLAGS   = -Wl,-O1 -Wl,--sort-common -Wl,--enable-new-dtags -Wl,--hash-style=both $(LDLIB)

SRCDIRS   = .
SRCEXTS = .c# .C .cc .cpp .CPP .c++ .cxx .cp

OLevel= -O3
CFLAGS  = -pipe -march=core2 -mtune=generic -std=gnu99 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -pthread -I../../lib/htslib
WARNFLAGS =  -Wfloat-equal -Wall \
#-pedantic

ASMFLAGS = -S -fverbose-asm -g -masm=intel

SOURCES = $(foreach d,$(SRCDIRS),$(wildcard $(addprefix $(d)/*,$(SRCEXTS))))

.PHONY: all clean

all: $(SOURCES:.c=)
#$(patsubst $.c,%,$(SOURCES))
	@echo [@gcc $(CFLAGS) $(OLevel) $(LDFLAGS)]
#@echo [$(SOURCES)] to [$(SOURCES:.c=)].

%: %.c
	@mkdir -p asm
	@echo  [$(OLevel)]:${<} -\> ${@}
	gcc $(WARNFLAGS) $(CFLAGS) $(OLevel) $(LDFLAGS) -o $@ $< ../../lib/htslib/libhts.a
	@gcc $(CFLAGS) $(OLevel) $(LDFLAGS) $(ASMFLAGS) -o ./asm/$@.s $<
	-as -alhnd -o /dev/null ./asm/$@.s > ./asm/$@.asm
#@gcc $(CFLAGS) $(OLevel) $(LDFLAGS) $(ASMFLAGS) -o - $< | as -alhnd > ./asm/$@.asm
	@rm ./asm/$@.s

clean:
	-rm $(SOURCES:.c=)
	-rm -frv ./asm
