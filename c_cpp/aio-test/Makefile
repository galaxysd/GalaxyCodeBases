CFLAGS=-O2 -pipe -march=core2 --param l1-cache-size=32 --param l1-cache-line-size=64 --param l2-cache-size=2048 -mtune=generic

LDFLAGS=-Wl,-O1 -Wl,--as-needed -Wl,--sort-common -Wl,--enable-new-dtags -Wl,--hash-style=both

BUILD=-Wall -laio -lpthread

all:
	gcc $(CFLAGS) $(BUILD) $(LDFLAGS) -o aio-stress aio-stress.c
	@echo See  https://github.com/repoforge/rpms/tree/master/specs/aio-stress for details.

clean:
	rm aio-stress

