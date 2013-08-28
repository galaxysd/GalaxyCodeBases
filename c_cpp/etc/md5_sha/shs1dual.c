/*
 * shs1dual - multi digest code that works with NIST Secure Hash Standard-1
 *
 * @(#) $Revision: 13.2 $
 * @(#) $Id: shs1dual.c,v 13.2 2006/08/14 10:09:23 chongo Exp $
 * @(#) $Source: /usr/local/src/cmd/hash/RCS/shs1dual.c,v $
 *
 * Split our data into multiple byte index streams, digest them all
 * and output the digests on a line.
 *
 * This file was written by Landon Curt Noll.
 *
 * This code has been placed in the public domain.  Please do not
 * copyright this code.
 *
 * LANDON CURT NOLL DISCLAIMS ALL WARRANTIES WITH  REGARD  TO
 * THIS  SOFTWARE,  INCLUDING  ALL IMPLIED WARRANTIES OF MER-
 * CHANTABILITY AND FITNESS.  IN NO EVENT SHALL  LANDON  CURT
 * NOLL  BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
 * DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM  LOSS  OF
 * USE,  DATA  OR  PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 * NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR  IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * chongo (was here) /\oo/\
 * http://www.isthe.com/chongo/index.html
 *
 * Share and enjoy! :-)
 *
 * See shs1drvr.c for version and modification history.
 */

char *shs1dual_what="@(#)";	/* #(@) if checked in */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "shs1.h"

/* static declarations */
static void multiData(BYTE*, UINT, BYTE*, UINT, UINT, SHS1_INFO*);
static void multiStream(BYTE*, UINT, FILE*, UINT, SHS1_INFO*);
static void multiFile(BYTE*, UINT, char*, UINT, SHS1_INFO*);
static void multiOutput(char*, int, int, UINT, SHS1_INFO*);

/* dual test suite strings */
#define ENTRY(str) {(BYTE *)str, sizeof(str)-1}
struct dual_test {
    BYTE *data;		/* data or NULL to test */
    int len;		/* length of data */
} dual_test_data[] = {
    {NULL, 0},
    ENTRY(""),
    ENTRY("a"),
    ENTRY("aa"),
    ENTRY("aabbccddeeffgghhiijjkkllmmnnooppqqrrssttuuvvwwxxyyzz"),
    ENTRY("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"),
    ENTRY("chongo <Ich bin, du bist, aber ein Yit ist nicht!!! :-)> /\\../\\"),
    ENTRY("123456789 123456789 123456789 123456789 123456789 123456789 1234"),
    ENTRY("a123456789 123456789 123456789 123456789 123456789 123456789 1234")
};
#define MAX_DUAL_TEST_DATA (sizeof(dual_test_data)/sizeof(dual_test_data[0]))

/* dual test filenames */
char *dual_test_file[] = {
    "file1",
    "file2",
    "shs1.data",
    "/dev/null"
};
#define MAX_DUAL_TEST_FILE (sizeof(dual_test_file)/sizeof(dual_test_file[0]))

/* where the test files are located by default */
#if !defined(TLIB)
#define TLIB "."
#endif


/*
 * multiData - divide data into multiple sections
 *
 * given:
 *	pre_str		string prefix or NULL
 *	pre_len		length of pre_str
 *	inString	string to digest
 *	in_len		length of inString
 *	cnt		number of digests
 *	digs		array of digests, cnt elements long
 */
static void
multiData(BYTE *pre_str, UINT pre_len, BYTE *inString,
    UINT in_len, UINT cnt, SHS1_INFO *digs)
{
    BYTE **bufs;		/* byte arrays for digests */
    UINT *buflen;		/* bytes stored in bufs[i] */
    int len;			/* total length of pre_str and inString */
    UINT indx;			/* byte stream index */
    BYTE **n;			/* n[i] is next byte to use in bufs[i] */
    BYTE *p;
    UINT i;

    /*
     * determine lengths
     */
    len = (pre_str == NULL) ? 0 : pre_len;
    len += (inString == NULL) ? 0 : in_len;
    /* no strings, quick return */
    if (len == 0 || cnt <= 0) {
        return;
    }

    /*
     * malloc string arrays
     */
    bufs = (BYTE **)malloc(sizeof(BYTE *)*cnt);
    if (bufs == NULL) {
	fprintf(stderr, "%s: bad malloc #1\n", program);
	exit(51);
    }
    buflen = (UINT *)malloc(sizeof(UINT)*cnt);
    if (buflen == NULL) {
	fprintf(stderr, "%s: bad malloc #2\n", program);
	exit(52);
    }
    n = (BYTE **)malloc(sizeof(BYTE *)*cnt);
    if (n == NULL) {
	fprintf(stderr, "%s: bad malloc #3\n", program);
	exit(53);
    }
    for (i=0; i < cnt; ++i) {
	bufs[i] = (BYTE *)malloc(1+(len/cnt));
	if (bufs[i] == NULL) {
	    fprintf(stderr, "%s: bad malloc #4\n", program);
	    exit(54);
	}
	buflen[i] = 0;
	n[i] = bufs[i];
    }

    /*
     * divide the pre-string
     */
    indx = 0;
    if (pre_str != NULL) {
	for (p=pre_str, i=0; i < pre_len; ++i, ++p, indx++) {
	    *(n[(indx = ((indx >= cnt) ? 0 : indx))]++) = *p;
	    ++buflen[indx];
	}
    }

    /*
     * divide the string
     */
    if (inString != NULL) {
	for (p=inString, i=0; i < in_len; ++indx, ++i, ++p) {
	    *(n[(indx = ((indx >= cnt) ? 0 : indx))]++) = *p;
	    ++buflen[indx];
	}
    }

    /*
     * update arrays
     */
    for (i=0; i < cnt; ++i) {
	shs1Update(digs+i, bufs[i], buflen[i]);
    }

    /*
     * cleanup
     */
    free(buflen);
    free(n);
    for (i=0; i < cnt; ++i) {
	free(bufs[i]);
    }
    free(bufs);
}


/*
 * multiStream - divide a Stream into multiple sections
 *
 * given:
 *	pre_str		data prefix or NULL
 *	pre_len		length of pre_str
 *	stream		the stream to process
 *	cnt		number of digests
 *	digs		array of digests, cnt elements long
 */
static void
multiStream(BYTE *pre_str, UINT pre_len, FILE *stream,
    UINT cnt, SHS1_INFO *digs)
{
    BYTE data[SHS1_READSIZE];	/* our read buffer */
    int bytes;			/* bytes last read */
    BYTE **bufs;		/* byte arrays for digests */
    UINT *buflen;		/* bytes stored in bufs[i] */
    UINT indx;			/* byte stream index */
    BYTE *pbeyond;		/* beyond end of used data or pre_str */
    BYTE **n;			/* n[i] is next byte to use in bufs[i] */
    BYTE *p;
    UINT i;

    /*
     * no sets, quick return
     */
    if (cnt <= 0) {
        return;
    }

    /*
     * malloc string arrays
     */
    bufs = (BYTE **)malloc(sizeof(BYTE *)*cnt);
    if (bufs == NULL) {
	fprintf(stderr, "%s: bad malloc #5\n", program);
	exit(55);
    }
    buflen = (UINT *)malloc(sizeof(UINT)*cnt);
    if (buflen == NULL) {
	fprintf(stderr, "%s: bad malloc #6\n", program);
	exit(56);
    }
    n = (BYTE **)malloc(sizeof(BYTE *)*cnt);
    if (n == NULL) {
	fprintf(stderr, "%s: bad malloc #7\n", program);
	exit(57);
    }
    for (i=0; i < cnt; ++i) {
	bufs[i] = (BYTE *)malloc(SHS1_BLOCKSIZE);
	if (bufs[i] == NULL) {
	    fprintf(stderr, "%s: bad malloc #8\n", program);
	    exit(58);
	}
	buflen[i] = 0;
	n[i] = bufs[i];
    }

    /*
     * divide the pre-string
     */
    indx = 0;
    if (pre_str != NULL && pre_len > 0) {
	for (p=pre_str, pbeyond=pre_str+pre_len; p < pbeyond; ++p, ++indx) {
	    *(n[(indx = ((indx >= cnt) ? 0 : indx))]++) = *p;
	    if (++buflen[indx] >= SHS1_BLOCKSIZE) {
		shs1Update(digs+indx, bufs[indx], SHS1_BLOCKSIZE);
		buflen[indx] = 0;
		n[indx] = bufs[indx];
	    }
	}
    }

    /*
     * process the contents of the file
     */
    while ((bytes = fread((char *)data, 1, SHS1_READSIZE, stream)) > 0) {

	/*
	 * load the bytes into the bufs
	 */
	for (p=data, pbeyond=data+bytes; p < pbeyond; ++p, ++indx) {
	    *(n[(indx = ((indx >= cnt) ? 0 : indx))]++) = *p;
	    if (++buflen[indx] >= SHS1_BLOCKSIZE) {
		shs1Update(digs+indx, bufs[indx], SHS1_BLOCKSIZE);
		buflen[indx] = 0;
		n[indx] = bufs[indx];
	    }
	}
    }

    /*
     * process any partial buffers
     */
    for (i=0; i < cnt; ++i) {
	if (buflen[i] > 0) {
	    shs1Update(digs+i, bufs[i], buflen[i]);
	}
    }

    /*
     * cleanup
     */
    free(buflen);
    free(n);
    for (i=0; i < cnt; ++i) {
	free(bufs[i]);
    }
    free(bufs);
}


/*
 * multiFile - divide a file into alternating bytes and digest both halves
 *
 * given:
 *	pre_str		string prefix or NULL
 *	pre_len		length of pre_str
 *	filename	the filename to process
 *	cnt		number of digests
 *	digs		array of digests, cnt elements long
 */
static void
multiFile(BYTE *pre_str, UINT pre_len, char *filename,
    UINT cnt, SHS1_INFO *digs)
{
    FILE *inFile;		/* the open file stream */
    struct stat buf;		/* stat or lstat of file */
    struct shs1_stat hashbuf;	/* stat data to digest */
    struct shs1_stat hashlbuf;	/* lstat data to digest */
    UINT i;

    /*
     * firewall
     */
    if (cnt <= 0) {
	return;
    }

    /*
     * open the file
     */
    inFile = fopen(filename, "rb");
    if (inFile == NULL) {
	fprintf(stderr, "%s: cannot open %s: ", program, filename);
	perror("");
	return;
    }

    /*
     * pre-process prefix if needed
     */
    if (pre_str == NULL || pre_len <= 0) {
	if (i_flag) {
	    multiData(NULL, 0, (BYTE *)filename, strlen(filename), cnt, digs);
	}
    } else {
	if (i_flag) {
	    multiData(pre_str, pre_len, (BYTE *)filename, strlen(filename),
		     cnt, digs);
	} else {
	    multiData(pre_str, pre_len, NULL, 0, cnt, digs);
	}
    }

    /*
     * digest file stat and lstat
     */
    if (i_flag) {
	if (fstat(fileno(inFile), &buf) < 0) {
	    printf("%s can't be stated.\n", filename);
	    return;
	}
	hashbuf.stat_dev = buf.st_dev;
	hashbuf.stat_ino = buf.st_ino;
	hashbuf.stat_mode = buf.st_mode;
	hashbuf.stat_nlink = buf.st_nlink;
	hashbuf.stat_uid = buf.st_uid;
	hashbuf.stat_gid = buf.st_gid;
	hashbuf.stat_size = buf.st_size;
	hashbuf.stat_mtime = buf.st_mtime;
	hashbuf.stat_ctime = buf.st_ctime;
	if (lstat(filename, &buf) < 0) {
	    printf("%s can't be lstated.\n", filename);
	    return;
	}
	hashlbuf.stat_dev = buf.st_dev;
	hashlbuf.stat_ino = buf.st_ino;
	hashlbuf.stat_mode = buf.st_mode;
	hashlbuf.stat_nlink = buf.st_nlink;
	hashlbuf.stat_uid = buf.st_uid;
	hashlbuf.stat_gid = buf.st_gid;
	hashlbuf.stat_size = buf.st_size;
	hashlbuf.stat_mtime = buf.st_mtime;
	hashlbuf.stat_ctime = buf.st_ctime;
        multiData((BYTE *)&hashbuf, sizeof(hashbuf), (BYTE *)&hashlbuf,
        	 sizeof(hashlbuf), cnt, digs);

	/*
	 * pad sections with zeros to process file data faster
	 */
	for (i=0; i < cnt; ++i) {
	    if (digs[i].datalen > 0) {
		shs1Update(digs+i, (BYTE *)shs1_zero,
			   SHS1_CHUNKSIZE - digs[i].datalen);
	    }
	}
    }

    /*
     * process the data stream
     */
    multiStream(NULL, 0, inFile, cnt, digs);
    fclose(inFile);
}


/*
 * multiOutput - output the multiple digests
 *
 * given:
 *	str		print string after digest, NULL => none
 *	quot		1 => surround str with a double quotes
 *	nospace		1 => don't space seperate multi digests
 *	cnt		number of digests
 *	digs		array of digests, cnt elements long
 */
static void
multiOutput(char *str, int quot, int nospace, UINT cnt, SHS1_INFO *digs)
{
    UINT i;

    /*
     * firewall
     */
    if (cnt <= 0) {
	return;
    }

    /*
     * finalize the sets
     */
    for (i=0; i < cnt; ++i) {
	shs1Final(digs+i);
    }

    /*
     * print the digests
     */
    shs1Print(c_flag, 0, digs);
    for (i=1; i < cnt-1; ++i) {
	if (nospace == 0) {
	    putchar(' ');
	    shs1Print(c_flag, 0, digs+i);
	} else {
	    shs1Print(0, 0, digs+i);
	}
    }
    if (i < cnt) {
	if (nospace == 0) {
	    putchar(' ');
	    shs1Print(c_flag, i_flag, digs+cnt-1);
	} else {
	    shs1Print(0, i_flag, digs+cnt-1);
	}
    }
    if (str && !q_flag) {
	if (quot) {
	    printf(" \"%s\"\n", str);
	} else {
	    printf(" %s\n", str);
	}
    } else {
	putchar('\n');
    }
}


/*
 * multiTest - shs1 dual test suite
 */
void
multiTest(void)
{
    struct dual_test *t;	/* current dual test */
    struct dual_test *p;	/* current dual pre-string test */
    struct stat buf;		/* stat of a test file */
    SHS1_INFO digs[2];		/* even byte digest */
    char **f;			/* current file being tested */
    unsigned int i;
    unsigned int j;

    /*
     * find all of the test files
     */
    for (i=0, f=dual_test_file; i < MAX_DUAL_TEST_FILE; ++i, ++f) {
	if (stat(*f, &buf) < 0) {
	    /* no file1 in this directory, cd to the test suite directory */
	    if (chdir(TLIB) < 0) {
		fflush(stdout);
		fprintf(stderr,
		    "%s: cannot find %s or %s/%s\n", program, *f, TLIB, *f);
		return;
	    }
	}
    }

    /*
     * try all combinations of test strings as prefixes and data
     */
    for (i=0, t=dual_test_data; i < MAX_DUAL_TEST_DATA; ++i, ++t) {
	for (j=1, p=dual_test_data+1; j < MAX_DUAL_TEST_DATA; ++j, ++p) {
	    printf("pre:%u data:%u\n", i, j);
	    shs1Init(digs+0);
	    shs1Init(digs+1);
	    multiData(p->data, p->len, t->data, t->len, 2, digs);
	    multiOutput(NULL, 0, 0, 2, digs);
	}
    }

    /*
     * try the files with all test strings as prefixes
     */
    for (i=0, p=dual_test_data; i < MAX_DUAL_TEST_DATA; ++i, ++p) {
	for (j=0, f=dual_test_file; j < MAX_DUAL_TEST_FILE; ++j, ++f) {
	    printf("pre:%d file:%s\n", i, *f);
	    shs1Init(digs+0);
	    shs1Init(digs+1);
	    multiFile(p->data, p->len, *f, 2, digs);
	    multiOutput(NULL, 0, 0, 2, digs);
	}
    }
    exit(0);
}


/*
 * multiMain - main driver of shs1 dual routines
 *
 * given:
 *	argc		arg count left after getopt
 *	argv		args left after getopt
 *	pre_str		pre-process this data first
 *	pre_len		length of pre_str
 *	data_str	data is this string, not a file
 *	nospace		1 => don't space seperate multi digests
 *	cnt		number of digests to perform
 */
void
multiMain(int argc, char **argv, BYTE *pre_str, UINT pre_len,
    char *data_str, int nospace, UINT cnt)
{
    extern int optind;		/* option index */
    SHS1_INFO *digs;		/* multiple digest */
    unsigned int i;

    /*
     * firewall
     */
    if (cnt <= 0) {
	return;
    }

    /*
     * initialize multiple digests
     */
    digs = (SHS1_INFO *)malloc(sizeof(SHS1_INFO)*cnt);
    if (digs == NULL) {
	fprintf(stderr, "%s: bad malloc #1\n", program);
	exit(60);
    }
    for (i=0; i < cnt; ++i) {
	shs1Init(digs+i);
    }

    /*
     * digest a string
     */
    if (data_str != NULL) {
	multiData(pre_str, pre_len, (BYTE*)data_str, strlen(data_str),
		  cnt, digs);
	multiOutput(data_str, 1, nospace, cnt, digs);

    /*
     * case: digest stdin
     */
    } else if (optind == argc) {
	multiStream(pre_str, pre_len, stdin, cnt, digs);
	multiOutput(NULL, 0, nospace, cnt, digs);

    /*
     * case: digest files
     */
    } else {
	for (; optind < argc; optind++) {
	    multiFile(pre_str, pre_len, argv[optind], cnt, digs);
	    multiOutput(argv[optind], 0, nospace, cnt, digs);
	}
    }
}
