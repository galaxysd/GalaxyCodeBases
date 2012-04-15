#include <stdlib.h> //EXIT_FAILURE
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#define handle_error(msg) \
   do { perror(msg); exit(EXIT_FAILURE); } while (0)

/*
[HTTP/1.1 404 File Not Found
Content-Type: text/octet
User-ReturnCode: -29104
Content-Length: 16

«Í˜vÿÿŽP]

Offset      0  1  2  3  4  5  6  7   8  9  A  B  C  D  E  F

00100000   48 54 54 50 2F 31 2E 31  20 34 30 34 20 46 69 6C   HTTP/1.1 404 Fil
00100010   65 20 4E 6F 74 20 46 6F  75 6E 64 0D 0A 43 6F 6E   e Not Found  Con
00100020   74 65 6E 74 2D 54 79 70  65 3A 20 74 65 78 74 2F   tent-Type: text/
00100030   6F 63 74 65 74 0D 0A 55  73 65 72 2D 52 65 74 75   octet  User-Retu
00100040   72 6E 43 6F 64 65 3A 20  2D 32 39 31 30 34 0D 0A   rnCode: -29104
00100050   43 6F 6E 74 65 6E 74 2D  4C 65 6E 67 74 68 3A 20   Content-Length:
00100060   31 36 0D 0A 0D 0A AB CD  98 76 FF FF 8E 50         16    «Í˜vÿÿŽP
*/
char theStr[]="HTTP/1.1 404 File Not Found\r\nContent-Type: text/octet\r\nUser-ReturnCode: -29104\r\nContent-Length: 16\r\n\r\n\xAB\xCD\x98\x76\xFF\xFF\x8E\x50";
int theStrLen = sizeof(theStr) - 1;

int main (int argc, char *argv[]) {
	if (argc!=4) {
		fprintf(stderr,"Usage: %s <main file> <2nd file> <out file>\n",argv[0]);
		exit(EXIT_FAILURE);
	}
    printf("%d %ld [%s]\n",theStrLen,strlen(theStr),theStr);
    int in1, in2, out;
    struct stat sb1, sb2;
    char *p1, *p2;
    in1 = open(argv[1], O_RDONLY);
    if (in1 == -1)
        handle_error("open1");
    in2 = open(argv[2], O_RDONLY);
    if (in2 == -1)
        handle_error("open2");

    if (fstat(in1, &sb1) == -1)           /* To obtain file size */
        handle_error("fstat1");
    if (fstat(in2, &sb2) == -1)           /* To obtain file size */
        handle_error("fstat2");
    if (sb1.st_size != sb2.st_size) {
		fprintf(stderr,"Input files with different size: %zd,%zd .\n",sb1.st_size,sb2.st_size);
		exit(EXIT_FAILURE);
    }
    if ( !(sb1.st_size && sb2.st_size) ) {
		fprintf(stderr,"Input files are empty.\n");
		exit(EXIT_FAILURE);
    }

    p1 = mmap(NULL, sb1.st_size, PROT_READ, MAP_PRIVATE, in1, 0);
    if (p1 == MAP_FAILED)
        handle_error("mmap1");
    p2 = mmap(NULL, sb2.st_size, PROT_READ, MAP_PRIVATE, in2, 0);
    if (p2 == MAP_FAILED)
        handle_error("mmap2");

    out = open(argv[3], O_WRONLY|O_CREAT|O_EXCL, S_IRUSR|S_IWUSR|S_IRGRP);
    if (out == -1)
        handle_error("open3");

    char *thisp, *lastp;
    thisp = lastp = p1;
    while ( thisp - p1 < sb1.st_size ) {
        if (*thisp == 'H') {
            if ( memcmp(thisp,theStr,theStrLen) == 0 ) {
                printf("1(%zx)   %zd\t%zx\n",(size_t)p1,thisp-p1,(size_t)thisp);
            }
        }
        ++thisp;
    }

    if (close(out) == -1)           /* To obtain file size */
        handle_error("close3");
    exit(EXIT_SUCCESS);
}
