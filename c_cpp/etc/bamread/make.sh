#!/bin/sh

rm bamread

gcc -Wfloat-equal -Wall  -pipe -march=core2 -mtune=generic -std=gnu99 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -pthread -I../../lib/htslib -O3 -Wl,-O1 -Wl,--sort-common -Wl,--enable-new-dtags -Wl,--hash-style=both -lz -o bamread bamread.c ../../lib/htslib/libhts.a

