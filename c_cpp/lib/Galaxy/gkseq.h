/* The MIT License
	From Heng Li <lh3@sanger.ac.uk>, as in bwa package.
	Macro translated by Hu Xuesong.

	Changelog:
	20110817  > Add kseq_t *kseq_open(char *const filename);
	             To open:
	              kseq_t* kseq = kseq_open("filename");
	              // gzFile fp=gzopen("filename", "r"); kseq_t* kseq=kseq_init(fp);
	             To close:
	              kseq_destroy(kseq); // gzclose(fp);
	          > Add RealSize,RealOffset as kseq->f->*
*/

#ifndef _G_KSEQ_H
#define _G_KSEQ_H

#include <ctype.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
//open
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
//lseek
#include <unistd.h>

#define __GKSEQ_FILETYPE gzFile	// Well, EP enough to use define here.
#define __GKSEQ_READFUNC gzread	// Maybe "lh3" just want to make debug harder for his code ?
#define __GKSEQ_OPENFUNC gzopen
#define __GKSEQ_CLOSEFUNC gzclose
//#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
//#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)
#define __GKSEQ_BUFSIZE (4096)

#ifndef KROUNDUP32
#define KROUNDUP32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif



// from gzip-1.4/gzip.h, added __gkseq_ / __GKSEQ_ prefix.
typedef unsigned char  __gkseq_uch;
typedef unsigned short __gkseq_ush;
typedef unsigned long  __gkseq_ulg;
/* Macros for getting two-byte and four-byte header values */
#define __GKSEQ_SH(p) ((__gkseq_ush)(__gkseq_uch)((p)[0]) | ((__gkseq_ush)(__gkseq_uch)((p)[1]) << 8))
#define __GKSEQ_LG(p) ((__gkseq_ulg)(__GKSEQ_SH(p)) | ((__gkseq_ulg)(__GKSEQ_SH((p)+2)) << 16))
// from gzip-1.4/gzip.h



//#74 "t.h"     // append `KSEQ_INIT(gzFile, gzread)` to kseq.h, remove all `#include`, you will get "t.h".
#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
    size_t l, m;
    char *s;
} kstring_t;
#endif

//#215 "t.h"    // gcc -E t.h > gkseq.h && indent -kr -brf -l120 -nut t.h -o gkseqn.h

/* RealSize: original (uncompressed) size of sequence file
   RealOffset: file offset in RealSize. Begin from 0, should pointing to '>' or '@'.
   ReadBuffer: current buffer count. The first buffer is 0. (init value is -1)
 */
typedef struct __kstream_t {
    char *buf;
    int begin, end, is_eof;
    unsigned long RealSize,RealOffset;
    long ReadBuffer;
    __GKSEQ_FILETYPE f;
} kstream_t;

typedef struct {
    kstring_t name, comment, seq, qual;
    int last_char;
    kstream_t *f;
} kseq_t;

static inline kstream_t *ks_init(__GKSEQ_FILETYPE f) {
    kstream_t *ks = (kstream_t *) calloc(1, sizeof(kstream_t));
    ks->f = f;
    ks->buf = (char *) malloc(__GKSEQ_BUFSIZE);
    ks->ReadBuffer=-1L;
    return ks;
}

static inline void ks_destroy(kstream_t * ks) {
    if (ks) {
        free(ks->buf);
        __GKSEQ_CLOSEFUNC(ks->f);
        free(ks);
    }
}

static inline int ks_getc(kstream_t * ks) {
    if (ks->is_eof && ks->begin >= ks->end)
        return -1;
    if (ks->begin >= ks->end) {
        ks->begin = 0;
        ks->end = __GKSEQ_READFUNC(ks->f, ks->buf, __GKSEQ_BUFSIZE);
        ++ks->ReadBuffer;
        if (ks->end < __GKSEQ_BUFSIZE)
            ks->is_eof = 1;
        if (ks->end == 0)
            return -1;
    }
    return (int) ks->buf[ks->begin++];
}

static int ks_getuntil(kstream_t * ks, int delimiter, kstring_t * str, int *dret) {
    if (dret)
        *dret = 0;
    str->l = 0;
    if (ks->begin >= ks->end && ks->is_eof)
        return -1;
    for (;;) {
        int i;
        if (ks->begin >= ks->end) {
            if (!ks->is_eof) {
                ks->begin = 0;
                ks->end = __GKSEQ_READFUNC(ks->f, ks->buf, __GKSEQ_BUFSIZE);
                ++ks->ReadBuffer;
                if (ks->end < __GKSEQ_BUFSIZE)
                    ks->is_eof = 1;
                if (ks->end == 0)
                    break;
            } else
                break;
        }
        if (delimiter) {
            for (i = ks->begin; i < ks->end; ++i)
                if (ks->buf[i] == delimiter)
                    break;
        } else {
            for (i = ks->begin; i < ks->end; ++i)
                if (isspace(ks->buf[i]))
                    break;
        }
        if (str->m - str->l < i - ks->begin + 1) {
            str->m = str->l + (i - ks->begin) + 1;
            KROUNDUP32(str->m);
            str->s = (char *) realloc(str->s, str->m);
        }
        memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin);
        str->l = str->l + (i - ks->begin);
        ks->begin = i + 1;
        if (i < ks->end) {
            if (dret)
                *dret = ks->buf[i];
            break;
        }
    }
    str->s[str->l] = '\0';
    return str->l;
}

static inline kseq_t *kseq_init(__GKSEQ_FILETYPE fd) {
    kseq_t *s = (kseq_t *) calloc(1, sizeof(kseq_t));
    s->f = ks_init(fd);
/*  realloc() may change the address, 
    thus, even we `malloc(32)` here and set `s->seq.m`, 
    `s->seq.s` still may be changed later on realloc().
*/
    return s;
}
static kseq_t *kseq_open(char *const filename) {
    unsigned long realsize=0;
	int ifd = open(filename, O_RDONLY);
	if (ifd==-1) return NULL;
	unsigned char buf[4];
	off_t bytes_in = lseek(ifd, (off_t)(-4), SEEK_END);
        if ( bytes_in != -1L) {
            bytes_in += 4L;
            if (read(ifd, (char*)buf, sizeof(buf)) == sizeof(buf)) {
                realsize = __GKSEQ_LG(buf);
            }
        }
	close(ifd);
    __GKSEQ_FILETYPE fp;
    fp = __GKSEQ_OPENFUNC(filename, "r");
    kseq_t *s = kseq_init(fp);
    s->f->RealSize = realsize;
    return s;
}

static inline void kseq_rewind(kseq_t * ks) {
    ks->last_char = 0;
    ks->f->is_eof = ks->f->begin = ks->f->end = 0;
}

static inline void kseq_destroy(kseq_t * ks) {
    if (!ks) return;
    free(ks->name.s); free(ks->comment.s); free(ks->seq.s);	free(ks->qual.s);
    ks_destroy(ks->f);
    free(ks);
}

/* Return value:
//   >=0  length of the sequence (normal)
    1   FASTA, just bases
    3   FASTQ, both bases and Q
   -1   end-of-file
   -2   truncated quality string
*/
static int_fast8_t kseq_read(kseq_t * seq) {	// better to be ssize_t
    int c;
    //uint_fast8_t with_last_chr = 0;
    kstream_t *ks = seq->f;
    if (seq->last_char == 0) {	// then jump to the next header line
        while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@');
        if (c == -1)
            return -1;	// end of file
        seq->last_char = c;
    }	// the first header char has been read
    seq->comment.l = seq->seq.l = seq->qual.l = 0;
    if (ks_getuntil(ks, 0, &seq->name, &c) < 0)
        return -1;
    if (c != '\n')
        ks_getuntil(ks, '\n', &seq->comment, 0);
    while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') {
        if (isgraph(c)) {	// printable non-space character
            if (seq->seq.l + 1 >= seq->seq.m) {	// double the memory
                seq->seq.m = seq->seq.l + 2;
                KROUNDUP32(seq->seq.m);	// rounded to next closest 2^k
                seq->seq.s = (char *) realloc(seq->seq.s, seq->seq.m);
            }
            seq->seq.s[seq->seq.l++] = (char) c;
        }
    }
    if (c == '>' || c == '@')
        seq->last_char = c;	// the first header char has been read
        //with_last_chr = 1;
    seq->seq.s[seq->seq.l] = 0;	// null terminated string
    if (c != '+') {
        //seq->qual.s[0] = 0;	// set Q to "". also bool, 0 for no Q. JUST TEMP, should change back later !!!
        seq->f->RealOffset = __GKSEQ_BUFSIZE*(seq->f->ReadBuffer) + seq->f->begin - 1;//with_last_chr;
        return 1;	// FASTA, only bases
    }
    if (seq->qual.m < seq->seq.m) {	// allocate enough memory
        seq->qual.m = seq->seq.m;
        seq->qual.s = (char *) realloc(seq->qual.s, seq->qual.m);
    }
    while ((c = ks_getc(ks)) != -1 && c != '\n');	// skip the rest of '+' line
    if (c == -1)
        return -2;	// we should not stop here
    while ((c = ks_getc(ks)) != -1 && seq->qual.l < seq->seq.l)
        if (c >= 33 && c <= 127)
            seq->qual.s[seq->qual.l++] = (unsigned char) c;
    seq->qual.s[seq->qual.l] = 0;	// null terminated string
    seq->last_char = 0;	// we have not come to the next header line
    //with_last_chr = 0;
    seq->f->RealOffset = __GKSEQ_BUFSIZE*(seq->f->ReadBuffer)+seq->f->begin;
    if (seq->seq.l != seq->qual.l)
        return -2;	// qual string is shorter than seq string
    //return seq->seq.l;
    return 3;   // both base and Q
}

#endif
