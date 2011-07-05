CC = gcc
CFLAGS = -Wall -pipe -march=core2 -mtune=generic -pedantic -Wfloat-equal
LDFLAGS = -lz -Wl,-O1 -Wl,--as-needed -Wl,--sort-common -Wl,--enable-new-dtags -Wl,--hash-style=both 
OPT = -O3

OBJDIR = ./objects/

MAKEARG = $(CFLAGS) $(OPT) 

objects = $(addprefix $(OBJDIR), gFileIO.o 2bitseq.o)
mainobjects = $(addprefix $(OBJDIR), example.o )

all: clean $(objects) example libgFileIO.a

example: libgFileIO.a $(mainobjects)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o example $(mainobjects) libgFileIO.a

$(objects) $(mainobjects): $(OBJDIR)%.o: %.c
	@$(CC) -std=gnu99 $(MAKEARG) -c $< -o $@

libgFileIO.a: $(objects)
	@touch libgFileIO.a
	@-rm libgFileIO.a 2> /dev/null
	@ar r libgFileIO.a $(objects) 

.PHONY : clean
clean:
	@-rm $(exefiles) $(mainobjects) $(objects) libgFileIO.a example 2> /dev/null
	@mkdir -p $(OBJDIR)

