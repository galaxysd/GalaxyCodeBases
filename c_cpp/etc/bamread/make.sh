#!/bin/sh

rm bamread

gcc -Wfloat-equal -Wall  -pipe -march=core2 -mtune=generic -std=gnu99 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -pthread -I./htslib -O3 -lz -o bamread bamread.c ./htslib/libhts.a
