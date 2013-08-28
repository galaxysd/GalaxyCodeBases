/*
 * md5drvr - md5 driver code
 *
 * @(#) $Revision: 13.3 $
 * @(#) $Id: md5drvr.c,v 13.3 2006/08/14 10:21:59 chongo Exp $
 * @(#) $Source: /usr/local/src/cmd/hash/RCS/md5drvr.c,v $
 *
 * This file was written by RSA Data Security.
 *
 * This file was modified by Landon Curt Noll.
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
 ***
 *
 * NOTE: The version information below refers to all md5 code, not
 *	 just this file.  In particular, this file was created by
 *	 Landon Curt Noll.
 *
 * Version 1.1: 17 Feb 1990		RLR
 *	Original code written.
 *
 * Version 1.2: 27 Dec 1990		SRD,AJ,BSK,JT
 *	C reference version.
 *
 * Version 1.3: 27 Apr 1991		RLR
 *	G modified to have y&~z instead of y&z
 *	FF, GG, HH modified to add in last register done
 *	access pattern: round 2 works mod 5, round 3 works mod 3
 *	distinct additive constant for each step
 *	round 4 added, working mod 7
 *
 * Version 1.4: 10 Jul 1991		SRD,AJ,BSK,JT
 *	Constant correction.
 *
 * Version 2.1:	31 Dec 1993		Landon Curt Noll
 *     Modified/Re-wrote md5.c
 *
 * Version 2.2: 09 Jan 1994		Landon Curt Noll
 *     md5drvr.c and md5dual.c code cloned from shs version 2.5.8 94/01/09
 *     performance tuning
 *
 * Version 2.3: 10 Jan 1994		Landon Curt Noll
 *     added MUST_ALIGN for Sparc and other RISC architectures
 *     use must_align.c to automatically determine if MUST_ALIGN is required
 *     performance tuning
 *     increased test to 64 megabytes due to increased performance
 *
 * Version 2.4:	non-existent version
 *
 * Version 2.5:	non-existent version
 *
 * Version 2.6: 10 Jan 1994		Landon Curt Noll
 *     Merged the shs and md5 Makefiles to build both in the same directory
 *     Bumped version to 2.6 to match level to shs
 *     Test suite header now says md5 (not MD5)
 *     Minor performance improvements
 *
 * Version 2.7: 14 Jan 1994		Landon Curt Noll
 *     code cleanup
 *     chunk is now 64 bytes, block is determined by blocking factor
 *     magic 64 and 64 related values defined in terms of #defines
 *     fixed bit count carry bug
 *     fixed writable strings test bug
 *
 * Version 2.8: 22 Jan 1994		Landon Curt Noll
 *     code cleanup
 *     count bytes in driver, convert to 64 bit count in final transform
 *     handle read errors and EOF better
 *     prefix strings not multiple of 64 bytes in length do not slow down hash
 *     renumbered exit codes
 *     fixed dual digest split bug
 *
 * Version 2.9: 05 Feb 1994		Landon Curt Noll
 *     prep for beta release
 *
 * Version 2.10: 25 Mar 1994		Landon Curt Noll
 *     must_align catchs signal to detect misaligned access
 *     malloc type not declared
 *
 * Version 2.11: 09 Mar 1995		Landon Curt Noll
 *     Moved stream and file routines to md5io.c.
 *
 * Version 2.12: 17 Nov 1995		Landon Curt Noll
 *     Fixed help string.
 *
 * Version 3.1:	non-existent version
 *
 * Version 4.1: 13 Aug 2006		Landon Curt Noll
 *	Bump to version 4.1 to match the version of SHA and SHA-1 drvr files.
 *	Port to ANSI C.
 *	Fixed all known compile warnings.
 *	Allow access to the internal hash transform function (to help code)
 *	  that makes direct use of the cryprographic hash source.
 *	Improved the way -v prints version.  Now -v prints the RCS version,
 *	  not version listed in the above comment.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "md5.h"

/* size of test in megabytes */
#define TEST_MEG 64

/* number of chunks to process */
#define TEST_CHUNKS (TEST_MEG*1024*1024/MD5_READSIZE)

/* MD5 test suite strings */
#define ENTRY(str) {(BYTE *)str, NULL, sizeof(str)-1}
struct MD5_test {
    BYTE *ro_data;	/* read only string data or NULL to test */
    BYTE *data;		/* data or NULL to test */
    int len;		/* length of data */
} MD5_test_data[] = {
    ENTRY(""),
    ENTRY("a"),
    ENTRY("abc"),
    ENTRY("message digest"),
    ENTRY("abcdefghijklmnopqrstuvwxyz"),
    ENTRY("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"),
    ENTRY("12345678901234567890123456789012345678901234567890123456789012345678901234567890")
};
#define MAX_MD5_TEST_DATA ((int)(sizeof(MD5_test_data)/sizeof(MD5_test_data[0])))

/* MD5 test filenames */
char *MD5_test_file[] = {
    "file1",
    "file2",
};
#define MAX_MD5_TEST_FILE ((int)(sizeof(MD5_test_file)/sizeof(MD5_test_file[0])))

/* where the test files are located by default */
#if !defined(TLIB)
#define TLIB "."
#endif

/* Prototypes of the static functions */
static void MD5Output(char*, int, MD5_CTX*);
static int MD5PreFileRead(char*, BYTE**);
static void MD5TestSuite(void);
static void MD5Help(void);

/* global variables */
static int c_flag = 0;		/* 1 => print C style digest with leading 0x */
int i_flag = 0;			/* 1 => process inode & filename */
int q_flag = 0;			/* 1 => print the digest only */
int dot_zero = 0;		/* 1 => print .0 after the digest */
int debug = 0;			/* 1 => add debug */
char *program = "<program name unknown>";	/* our name */


/*
 * MD5Output - output the digest
 *
 * given:
 *	str	print string after digest, NULL => none
 *	quot	1 => surround str with a double quotes
 *	dig	current digest
 */
static void
MD5Output(char *str, int quot, MD5_CTX *dig)
{
    /*
     * finalize the digest
     */
    MD5Final(dig);

    /*
     * print the digest
     */
    MD5Print(dig);
    if (str && !q_flag) {
	if (quot) {
	    printf(" \"%s\"\n", str);
	} else {
	    printf(" %s\n", str);
	}
    } else {
	putchar('\n');
    }
    fflush(stdout);
}


/*
 * MD5Print - print a digest in hex
 *
 * Prints message digest buffer in MD5Info as 40 hexadecimal digits. Order is
 * from low-order byte to high-order byte of digest. Each byte is printed
 * with high-order hexadecimal digit first.
 *
 * If -c, then print a leading "0x".  If -i, then print a trailing ".0".
 */
void
MD5Print(MD5_CTX *MD5Info)
{
    int i;

    if (c_flag) {
	fputs("0x", stdout);
    }
    for (i = 0; i < 16; i++) {
	printf ("%02x", MD5Info->digest[i]);
    }
    if (dot_zero) {
	fputs(".0", stdout);
    }
}


/*
 * A time trial routine, to measure the speed of MD5.
 *
 * Measures user time required to digest TEST_MEG megabytes of characters.
 */
static void
MD5TimeTrial(void)
{
    BYTE data[MD5_READSIZE];	/* test buffer */
    MD5_CTX MD5Info;		/* hash state */
    struct rusage start;	/* test start time */
    struct rusage stop;		/* test end time */
    double usrsec;		/* duration of test in user seconds */
    unsigned int i;

    /* initialize test data */
    for (i = 0; i < MD5_READSIZE; i++)
	data[i] = (BYTE)(i & 0xFF);

    /* start timer */
    if (!q_flag) {
	printf("md5 time trial for %d megabytes of test data ...", TEST_MEG);
	fflush(stdout);
    }
    getrusage(RUSAGE_SELF, &start);

    /* digest data in MD5_READSIZE byte chunk */
    MD5Init(&MD5Info);
    for (i=0; i < TEST_CHUNKS; ++i) {
	MD5fullUpdate(&MD5Info, data, MD5_READSIZE);
    }
    MD5COUNT(&MD5Info, MD5_READSIZE*TEST_CHUNKS);
    MD5Final(&MD5Info);

    /* stop timer, get time difference */
    getrusage(RUSAGE_SELF, &stop);
    usrsec = (stop.ru_utime.tv_sec - start.ru_utime.tv_sec) +
    	   (double)(stop.ru_utime.tv_usec - start.ru_utime.tv_usec)/1000000.0;
    if (!q_flag) {
	putchar('\n');
    }
    MD5Print(&MD5Info);
    if (q_flag) {
	putchar('\n');
    } else {
	printf(" is digest of test data\n");
	printf("user seconds to process test data: %.2f\n", usrsec);
	printf("characters processed per user second: %d\n",
	    (int)((double)TEST_MEG*1024.0*1024.0/usrsec));
    }
}


/*
 * Runs a standard suite of test data.
 */
static void
MD5TestSuite(void)
{
    struct MD5_test *t;		/* current MD5 test */
    struct stat buf;		/* stat of a test file */
    MD5_CTX digest;		/* test digest */
    char **f;			/* current file being tested */
    int i;

    /*
     * copy our test strings into writable data
     */
    for (i=0, t=MD5_test_data; i < MAX_MD5_TEST_DATA; ++i, ++t) {
	if (t->ro_data != NULL) {
	    t->data = (BYTE *)malloc(t->len + 1);
	    if (t->data == NULL) {
		fprintf(stderr, "%s: malloc #4 failed\n", program);
		exit(4);
	    }
	    strcpy((char *)t->data, (char *)t->ro_data);
        }
    }

    /*
     * print test header
     */
    puts("md5 test suite results:");

    /*
     * find all of the test files
     */
    for (i=0, f=MD5_test_file; i < MAX_MD5_TEST_FILE; ++i, ++f) {
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
    for (i=0, t=MD5_test_data; i < MAX_MD5_TEST_DATA; ++i, ++t) {
	MD5Init(&digest);
	MD5Update(&digest, t->data, t->len);
	MD5COUNT(&digest, t->len);
	MD5Output((char *)(t->data), 1, &digest);
    }

    /*
     * try the files with all test strings as prefixes
     */
    for (i=0, f=MD5_test_file; i < MAX_MD5_TEST_FILE; ++i, ++f) {
	MD5Init(&digest);
	MD5File(NULL, 0, *f, 0, &digest);
	MD5Output(*f, 0, &digest);
    }
    exit(0);
}


/*
 * MD5PreFileRead - read and process a prepend file
 *
 * given:
 *	pre_file		form pre_str from file pre_file
 *	buf			pointer to pre_str pointer
 *
 * Returns the length of pre_str, and modifies pre_str to
 * point at the malloced prepend data.
 */
static int
MD5PreFileRead(char *pre_file, BYTE **buf)
{
    struct stat statbuf;	/* stat for pre_file */
    int pre_len;		/* length of pre_file to be used */
    int bytes;			/* bytes read from pre_file */
    FILE *pre;			/* pre_file descriptor */

    /* obtain the length that we will use */
    if (stat(pre_file, &statbuf) < 0) {
	fprintf(stderr, "%s: unpable to find prepend file %s\n",
	    program, pre_file);
	exit(5);
    }
    pre_len = statbuf.st_size;
    if (pre_len > MD5_MAX_PRE_FILE) {
	/* don't use beyond MD5_MAX_PRE_FILE in size */
	pre_len = MD5_MAX_PRE_FILE;
    }

    /* malloc our pre string */
    *buf = (BYTE *)malloc(pre_len+1);
    if (*buf == NULL) {
	fprintf(stderr, "%s: malloc #3 failed\n", program);
	exit(6);
    }

    /* open our pre_file */
    pre = fopen(pre_file, "rb");
    if (pre == NULL) {
	fprintf(stderr, "%s: unable to open prepend file %s\n",
	  program, pre_file);
	exit(7);
    }

    /* read our pre_file data */
    bytes = fread((char *)(*buf), 1, pre_len, pre);
    if (bytes != pre_len) {
	fprintf(stderr,
	  "%s: unable to read %d bytes from prepend file %s\n",
	  program, pre_len, pre_file);
	exit(8);
    }

    /* return our length */
    return (pre_len);
}


/*
 * MD5Help - print MD5 help message and exit
 */
static void
MD5Help(void)
{
    fprintf(stderr,
      "%s [-cdhiqtx][-p str][-P str][-s str] [file ...]\n",
      program);
    fprintf(stderr,
      "    -c          print C style digests with a leading 0x\n");
    fprintf(stderr,
      "    -d          dual digests of even and odd indexed bytes\n");
    fprintf(stderr,
      "    -h          prints this message\n");
    fprintf(stderr,
      "    -i          process inode and filename as well as file data\n");
    fprintf(stderr,
      "    -p str      prepend str to data before digesting\n");
    fprintf(stderr,
      "    -P str      prepend the file 'str' to data before digesting\n");
    fprintf(stderr,
      "    -q          print only the digest\n");
    fprintf(stderr,
      "    -s str      prints digest and contents of string\n");
    fprintf(stderr,
      "    -t          prints time statistics for %dM chars\n", TEST_MEG);
    fprintf(stderr,
      "    -v          print version\n");
    fprintf(stderr,
      "    -x          execute an extended standard suite of test data\n");
    fprintf(stderr,
      "    file        print digest and name of file\n");
    fprintf(stderr,
      "    (no args)   print digest of stdin\n");
    exit(0);
}


/*
 * main - MD5 main
 */
int
main(int argc, char *argv[])
{
    MD5_CTX digest;		/* our current digest */
    BYTE *pre_str = NULL;	/* pre-process this data first */
    char *pre_file = NULL;	/* pre-process this file first */
    char *data_str = NULL;	/* data is this string, not a file */
    UINT pre_str_len;		/* length of pre_str or pre_file */
    UINT data_str_len;		/* length of data_str */
    int d_flag = 0;		/* 1 => dual digest mode */
    int t_flag = 0;		/* 1 => -t was given */
    int x_flag = 0;		/* 1 => -x was given */
    extern char *optarg;	/* argument to option */
    extern int optind;		/* option index */
    int c;

    /*
     * parse args
     */
    program = argv[0];
    while ((c = getopt(argc, argv, "cdihp:P:qs:tvx")) != -1) {
        switch (c) {
        case 'c':
	    c_flag = 1;
	    break;
        case 'd':
	    d_flag = 1;
	    break;
	case 'h':
	    MD5Help();
	    /*NOTREACHED*/
	    break;
	case 'i':
            i_flag = 1;
            break;
	case 'p':
	    pre_str = (BYTE *)optarg;
	    break;
	case 'q':
	    q_flag = 1;
	    break;
	case 'P':
	    pre_file = optarg;
	    break;
        case 's':
            data_str = optarg;
            break;
	case 't':
	    t_flag = 1;
	    break;
	case 'v':
	    printf("%s: version %s\n", program, MD5_VERSION);
	    exit(0);
	case 'x':
	    x_flag = 1;
	    break;
	default:
	    MD5Help();
	    break;
        }
    }
    if (data_str && optind != argc) {
	fprintf(stderr, "%s: -s is not compatible with digesting files\n",
	    program);
	exit(9);
    }
    if (i_flag && optind == argc) {
	fprintf(stderr, "%s: -i works only on filenames\n", program);
	exit(10);
    }

    /*
     * process -x if needed
     */
    if (x_flag) {
        if (d_flag) {
            dualTest();
	} else {
	    MD5TestSuite();
	}
	exit(0);
    }

    /*
     * process -t if needed
     */
    if (t_flag) {
	MD5TimeTrial();
	exit(0);
    }

    /*
     * process -P or -p if needed
     */
    if (pre_str && pre_file) {
	fprintf(stderr, "%s: -p and -P conflict\n", program);
	exit(11);
    }
    if (pre_file) {
	pre_str_len = MD5PreFileRead(pre_file, &pre_str);
    } else if (pre_str) {
        pre_str_len = strlen((char *)pre_str);
    } else {
        pre_str_len = 0;
    }

    /*
     * if -d, perform dual digest processing instead
     */
    if (d_flag) {
	dualMain(argc, argv, pre_str, pre_str_len, data_str);

    /*
     * if no -d, process string, stdin or files
     */
    } else {

	/*
	 * case: digest a string
	 */
	if (data_str != NULL) {
	    data_str_len = strlen(data_str);
	    MD5Init(&digest);
	    MD5Update(&digest, pre_str, pre_str_len);
	    MD5Update(&digest, (BYTE *)data_str, data_str_len);
	    MD5COUNT(&digest, pre_str_len + data_str_len);
	    MD5Output(data_str, 1, &digest);

	/*
	 * case: digest stdin
	 */
	} else if (optind == argc) {
	    MD5Init(&digest);
	    MD5Stream(pre_str, pre_str_len, stdin, &digest);
	    MD5Output(NULL, 0, &digest);

	/*
	 * case: digest files
	 */
	} else {
	    if (i_flag) {
		dot_zero = 1;
	    }
	    for (; optind < argc; optind++) {
		MD5Init(&digest);
		MD5File(pre_str, pre_str_len, argv[optind], i_flag, &digest);
		MD5Output(argv[optind], 0, &digest);
	    }
	}
    }

    /* all done */
    /* exit(0); */
    return 0;
}
