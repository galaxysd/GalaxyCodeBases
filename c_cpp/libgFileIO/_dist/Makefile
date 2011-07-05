CC = gcc
CFLAGS = -Wall -pipe -march=core2 -mtune=generic -pedantic -Wfloat-equal
LDFLAGS = -lz -Wl,-O1 -Wl,--as-needed -Wl,--sort-common -Wl,--enable-new-dtags -Wl,--hash-style=both 
OPT = -O3 -std=gnu99

MAKEARG = $(CFLAGS) $(OPT) 

all: clean example

example: libgFileIO.a
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o example example.c libgFileIO.a

.PHONY : clean
clean:
	@-rm example 2> /dev/null
