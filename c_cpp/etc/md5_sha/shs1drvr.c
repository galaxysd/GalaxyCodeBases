/*
 * shs1drvr - NIST Secure Hash Standard-1 implimentation driver code
 *
 * @(#) $Revision: 13.4 $
 * @(#) $Id: shs1drvr.c,v 13.4 2010/10/12 21:09:35 chongo Exp $
 * @(#) $Source: /usr/local/src/cmd/hash/RCS/shs1drvr.c,v $
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
 ***
 *
 * NOTE: The version information below refers to all shs1 code, not
 *	 just this file.  In particular, this file was created by
 *	 Landon Curt Noll.
 *
 * Version 1.1: 02 Sep 1992?			original authors
 *     This code is based on code by Peter C. Gutmann.  Much thanks goes
 *     to Peter C. Gutman (pgut1@cs.aukuni.ac.nz) , Shawn A. Clifford
 *     (sac@eng.ufl.edu), Pat Myrto (pat@rwing.uucp) and others who wrote
 *     and/or worked on the original code.
 *
 * Version 2.1:	31 Dec 1993		Landon Curt Noll
 *     Reformatted, performance improvements and bug fixes
 *
 * Version 2.2:	02 Jan 1994		Landon Curt Noll
 *     fixed -p usage
 *     better error messages
 *     added -c help
 *     added -c 0	(concatenation)
 *     reordered -i stat buffer pre-pending
 *
 * Version 2.3:	03 Jan 1994		Landon Curt Noll
 *     added -c 1	(side by side)
 *     added -c 2	(even force to be odd)
 *     added -c x	(shs1 dual test suite)
 *     changed -c help to be -c h
 *     changed -c operand to type[,opt[,...]]
 *     prefix & string ABI now can take arbitrary binary data
 *     fixed memory leak
 *     fixed even/odd byte split bug
 *     added -P file
 *     added -q
 *     added UNROLL_LOOPS to control shs1.c loop unrolling
 *     major performance improvements
 *
 * Version 2.4: 05 Jan 1994		Landon Curt Noll
 *     renamed calc mode to dual mode
 *     removed all -c code
 *     added -d		(dual digests, space separated)
 *     rewrote most of the file, string and stream logic using shs1dual code
 *
 * Version 2.5: 08 Jan 1994		Landon Curt Noll
 *     added (new) -c	(print 0x in front of digests)
 *     removed st_blksize and st_blocks from -i preprocessing data
 *     only print .0 suffix if -i and digesting a file
 *     non-zero edit codes are now unique
 *     changed the dual test suite (shorter, added non alpha numeric chars)
 *     -i requires filenames
 *     fixed @(#) what string code
 *     boolean logic simplication by Rich Schroeppel (rcs@cs.arizona.edu)
 *     on the fly in a circular buffer by Colin Plumb (colin@nyx10.cs.du.edu)
 *
 * Version 2.6: 11 Jan 1994		Landon Curt Noll
 *     Merged the shs1 and md5 Makefiles to build both in the same directory
 *     alignment and byte order performance improvements
 *     eliminate wateful memory copies
 *     shs1 transform contains no function calls
 *     beta release
 *
 * Version 2.7: 13 Jan 1994		Landon Curt Noll
 *     code cleanup
 *     chunk is now 64 bytes, block is determined by blocking factor
 *     magic 64 and 64 related values defined in terms of #defines
 *     added blocking code (-b block_len)
 *     added xor feedback code (-f)
 *     added xor feedback and block options to performance test
 *     performance improvements
 *
 * Version 2.8: 16 Jan 1994		Landon Curt Noll
 *     code cleanup
 *     performance improvements
 *     removed blocking and feedback code
 *     count bytes in driver, convert to 64 bit count in final transform
 *     added debug mode
 *     handle read errors and EOF better
 *     prefix strings not multiple of 64 bytes in length do not slow down hash
 *     renumbered exit codes
 *     fixed dual digest split bug
 *     byte sex swapping is now controlled thru the SHS1_TRANSFORM macro
 *     shs1Transform() is now called via the SHS1_TRANSFORM macro
 *
 * Version 2.9: 12 Feb 1994		Landon Curt Noll
 *     prep for beta release
 *     removed all feedback code
 *
 * Version 2.10: 25 Mar 1994		Landon Curt Noll
 *     must_align catchs signal to detect misaligned access
 *
 * Version 3.1: 09 Mar 1995		Landon Curt Noll
 *     Changed to implement the new Secure Hash Standard-1 (SHS1).
 *     The Secure Hash Standard-1 (SHS1) is a United States Department
 *     of Commerce National Institute of Standards and Technology approved
 *     standard (FIPS Pub 180-1).
 *
 *     The only substantial change was made to the exor() macro in shs1.c
 *     Changed name to shs1.  Bumped version from 2.10 to 3.1.
 *
 *     Moved stream and file routines into shs1io.c.
 *
 * Version 3.2: 17 Nov 1995		Landon Curt Noll
 *     Fixed help string.
 *
 *     Added multiple digests capability instead of just dual.  Added
 *     -m num to denote 2 or more multiple digests.
 *
 *     Added -C to prevent spaces (and later 0x if -c) between multi digests.
 *
 * Version 3.3: 01 Sep 1996		Landon Curt Noll
 *     Provide sha1 as well as shs1.
 *
 * Version 4.0: 13 Aug 2006		Landon Curt Noll
 *	Port to ANSI C.
 *	Fixed all known compile warnings.
 *	Allow access to the internal hash transform function.
 *	The shs1Transform() function performs the byte swapping now.
 *	Eliminated use of SHS1_TRANSFORM macro.
 *	Improved the way -v prints version.  Now -v prints the RCS version,
 *	  not version listed in the above comment.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "shs1.h"

/* size of test in megabytes */
#define TEST_MEG 16

/* number of chunks to process */
#define TEST_CHUNKS (TEST_MEG*1024*1024/SHS1_READSIZE)

/* SHS1 test suite strings */
#define ENTRY(str) {(BYTE *)str, NULL, sizeof(str)-1}
struct shs1_test {
    BYTE *ro_data;	/* read only string data or NULL to test */
    BYTE *data;		/* data or NULL to test */
    int len;		/* length of data */
} shs1_test_data[] = {
    ENTRY(""),
    ENTRY("a"),
    ENTRY("abc"),
    ENTRY("message digest"),
    ENTRY("abcdefghijklmnopqrstuvwxyz"),
    ENTRY("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"),
    ENTRY("12345678901234567890123456789012345678901234567890123456789012345678901234567890"),
    ENTRY("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq")
};
#define MAX_SHS1_TEST_DATA ((int)(sizeof(shs1_test_data)/sizeof(shs1_test_data[0])))

/* shs1 test filenames */
char *shs1_test_file[] = {
    "file1",
    "file2",
};
#define MAX_SHS1_TEST_FILE ((int)(sizeof(shs1_test_file)/sizeof(shs1_test_file[0])))

/* where the test files are located by default */
#if !defined(TLIB)
#define TLIB "."
#endif

/* Prototypes of the static functions */
static void shs1Output(char*, int, int, int, SHS1_INFO*);
static int shs1PreFileRead(char*, BYTE**);
static void shs1TestSuite(void);
static void shs1Help(void);

/* global variables */
int c_flag = 0;			/* 1 => print C style digest with leading 0x */
int i_flag = 0;			/* 1 => process inode & filename */
int q_flag = 0;			/* 1 => print the digest only */
int debug = 0;			/* 1 => add debug */
char *program = "<program name unknown>";	/* our name */


/*
 * shs1Output - output the digest
 *
 * given:
 *	str		print string after digest, NULL => none
 *	quot		1 => surround str with a double quotes
 *	hexhdr		1 => print 0x before the hex value
 *	trail		1 => print a trailing ".0"
 *	dig		current digest
 */
static void
shs1Output(char *str, int quot, int hexhdr, int trail, SHS1_INFO *dig)
{
    /*
     * finalize the digest
     */
    shs1Final(dig);
#if defined(DEBUG)
    if (debug) {
	fprintf(stderr, "DEBUG: octet count: %lld\n", dig->octets);
	fprintf(stderr, "DEBUG: octet count: 0x%016llx\n", dig->octets);
	fprintf(stderr, "DEBUG: bit count: %lld\n", dig->octets<<3);
	fprintf(stderr, "DEBUG: bit count: 0x%016llx\n", dig->octets<<3);
    }
#endif

    /*
     * print the digest
     */
    shs1Print(hexhdr, trail, dig);
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
 * shs1Print - print a digest in hex
 *
 * given:
 *	hexhdr		1 => print 0x before the hex value
 *	trail		1 => print a trailing ".0"
 *	shs1Info	SHA-1 hash state
 *
 * Prints message digest buffer in shs1Info as 40 hexadecimal digits. Order is
 * from low-order byte to high-order byte of digest. Each byte is printed
 * with high-order hexadecimal digit first.
 *
 * If hexhdr, then print a leading "0x".  If trail, then print a trailing ".0".
 */
void
shs1Print(int hexhdr, int trail, SHS1_INFO *shs1Info)
{
    if (hexhdr) {
	fputs("0x", stdout);
    }
    printf("%08x%08x%08x%08x%08x",
	shs1Info->digest[0], shs1Info->digest[1], shs1Info->digest[2],
	shs1Info->digest[3], shs1Info->digest[4]);
    if (trail) {
	fputs(".0", stdout);
    }
}


/*
 * shs1TimeTrial - measure the speed of SHS1
 *
 * Measures user time required to digest TEST_MEG megabytes of characters.
 *
 * This function will time blocking and under xor feedback mode if they
 * are set.
 */
static void
shs1TimeTrial(void)
{
    ULONG data[SHS1_READWORDS];	/* test buffer */
    SHS1_INFO shs1Info;		/* hash state */
    struct rusage start;	/* test start time */
    struct rusage stop;		/* test end time */
    double usrsec;		/* duration of test in user seconds */
    unsigned int i;

    /*
     * initialize test data
     */
    for (i = 0; i < SHS1_READSIZE; i++) {
	((BYTE *)data)[i] = (BYTE)(i & 0xFF);
    }

    /*
     * announce test
     */
    if (!q_flag) {
	printf("shs1 time trial for %d megs of test data ...", TEST_MEG);
	fflush(stdout);
    }

    /*
     * digest data in SHS1_READSIZE byte chunks
     */
    getrusage(RUSAGE_SELF, &start);
    shs1Init(&shs1Info);
    for (i=0; i < TEST_CHUNKS; ++i) {
	shs1fullUpdate(&shs1Info, (BYTE *)data, SHS1_READSIZE);
    }
    shs1Final(&shs1Info);
    getrusage(RUSAGE_SELF, &stop);

    /*
     * announce the test results
     */
    usrsec = (stop.ru_utime.tv_sec - start.ru_utime.tv_sec) +
    	   (double)(stop.ru_utime.tv_usec - start.ru_utime.tv_usec)/1000000.0;
    if (!q_flag) {
	putchar('\n');
    }
    shs1Print(0, 0, &shs1Info);
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
 * shs1TestSuite - run a standard suite of test data
 */
static void
shs1TestSuite(void)
{
    struct shs1_test *t;		/* current shs1 test */
    struct stat buf;		/* stat of a test file */
    SHS1_INFO digest;		/* test digest */
    char **f;			/* current file being tested */
    int i;

    /*
     * copy our test strings into writable data
     */
    for (i=0, t=shs1_test_data; i < MAX_SHS1_TEST_DATA; ++i, ++t) {
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
    puts("shs1 test suite results:");

    /*
     * find all of the test files
     */
    for (i=0, f=shs1_test_file; i < MAX_SHS1_TEST_FILE; ++i, ++f) {
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
    for (i=0, t=shs1_test_data; i < MAX_SHS1_TEST_DATA; ++i, ++t) {
	shs1Init(&digest);
	shs1Update(&digest, t->data, t->len);
	shs1Output((char *)t->ro_data, 1, c_flag, i_flag, &digest);
    }

    /*
     * try the files with all test strings as prefixes
     */
    for (i=0, f=shs1_test_file; i < MAX_SHS1_TEST_FILE; ++i, ++f) {
	shs1Init(&digest);
	shs1File(NULL, 0, *f, 0, &digest);
	shs1Output(*f, 0, c_flag, i_flag, &digest);
    }
    exit(0);
}


/*
 * shs1PreFileRead - read and process a prepend file
 *
 * given:
 *	pre_file		form pre_str from file pre_file
 *	buf			pointer to pre_str pointer
 *
 * Returns the length of pre_str, and modifies pre_str to
 * point at the malloced prepend data.
 */
static int
shs1PreFileRead(char *pre_file, BYTE **buf)
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
    if (pre_len > SHS1_MAX_PRE_FILE) {
	/* don't use beyond SHS1_MAX_PRE_FILE in size */
	pre_len = SHS1_MAX_PRE_FILE;
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
 * shs1Help - print shs1 help message and exit
 */
static void
shs1Help(void)
{
    fprintf(stderr,
      "%s [-cCd%shimqtx][-p prefix][-P pfile][-s str] [file ...]\n",
      program,
#if defined(DEBUG)
      "D"
#else
      ""
#endif
      );
    fprintf(stderr,
      "    -c          print C style digests with a leading 0x\n");
    fprintf(stderr,
      "    -C          do not space seperate multi digests\n");
    fprintf(stderr,
      "    -d          same as -m 2\n");
#if defined(DEBUG)
    fprintf(stderr,
      "    -D          debug mode\n");
#endif
    fprintf(stderr,
      "    -h          prints this message\n");
    fprintf(stderr,
      "    -i          process inode and filename as well as file data\n");
    fprintf(stderr,
      "    -m num      process num digests on every num-th byte\n");
    fprintf(stderr,
      "    -p prefix   prepend str to data before digesting\n");
    fprintf(stderr,
      "    -P pfile    prepend the file 'str' to data before digesting\n");
    fprintf(stderr,
      "    -q          print only the digest\n");
    fprintf(stderr,
      "    -r          reverse feedback mode\n");
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
 * main - shs1 main control function
 */
int
main(int argc, char *argv[])
{
    SHS1_INFO digest;		/* our current digest */
    BYTE *pre_str = NULL;	/* pre-process this data first */
    char *pre_file = NULL;	/* pre-process this file first */
    char *data_str = NULL;	/* data is this string, not a file */
    UINT pre_str_len;		/* length of pre_str or pre_file */
    UINT data_str_len;		/* length of data_str */
    int C_flag = 0;		/* 1 => don't space seperate multi digests */
    int d_flag = 0;		/* 1 => dual digest mode */
    int m_flag = 0;		/* 1 => multi digest mode */
    int t_flag = 0;		/* 1 => -t was given */
    int x_flag = 0;		/* 1 => -x was given */
    int multi = 0;		/* number of digests to do in parallel */
    extern char *optarg;	/* argument to option */
    extern int optind;		/* option index */
    int c;

    /*
     * parse args
     */
    program = argv[0];
    while ((c = getopt(argc, argv, "cCdDihm:p:P:qs:tvx")) != -1) {
        switch (c) {
        case 'c':
	    c_flag = 1;
	    break;
        case 'C':
	    C_flag = 1;
	    break;
        case 'd':
	    d_flag = 1;
	    multi = 2;
	    break;
        case 'D':
#if defined(DEBUG)
	    debug = 1;
#else
	    fprintf(stderr, "%s: not compiled with -DDEBUG\n", program);
	    exit(9);
	    /*NOTREACHED*/
#endif
	    break;
	case 'h':
	    shs1Help();
	    /*NOTREACHED*/
	    break;
	case 'i':
            i_flag = 1;
            break;
	case 'm':
	    m_flag = 1;
	    multi = atoi(optarg);
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
	    printf("%s: version %s\n", program, SHS1_VERSION);
	    exit(0);
	case 'x':
	    x_flag = 1;
	    break;
	default:
	    shs1Help();
	    break;
        }
    }
    /* arg checking */
    if (data_str && optind != argc) {
	fprintf(stderr, "%s: -s is not compatible with digesting files\n",
	    program);
	exit(10);
    }
    if (i_flag && optind == argc) {
	fprintf(stderr, "%s: -i works only on filenames\n", program);
	exit(11);
    }
    if (d_flag && m_flag) {
	fprintf(stderr, "%s: -d imples -m 2, use one or the other\n", program);
	exit(12);
    }
    if ((d_flag || m_flag) && multi < 2) {
	fprintf(stderr, "%s: multi count must be > 1\n", program);
	exit(13);
    }
    if (C_flag && !d_flag && !m_flag) {
	fprintf(stderr, "%s: -C requires -d or -m num\n", program);
	exit(14);
    }

    /*
     * process -x if needed
     */
    if (x_flag) {
        if (d_flag) {
            multiTest();
	} else {
	    shs1TestSuite();
	}
	exit(0);
    }

    /*
     * process -t if needed
     */
    if (t_flag) {
	shs1TimeTrial();
	exit(0);
    }

    /*
     * process -P or -p if needed
     */
    if (pre_str && pre_file) {
	fprintf(stderr, "%s: -p and -P conflict\n", program);
	exit(15);
    }
    if (pre_file) {
	pre_str_len = shs1PreFileRead(pre_file, &pre_str);
    } else if (pre_str) {
        pre_str_len = strlen((char *)pre_str);
    } else {
        pre_str_len = 0;
    }
    if (pre_str_len > SHS1_MAX_PRE_FILE) {
	fprintf(stderr, "%s: prefix may not be longer than %d bytes\n",
	    program, SHS1_MAX_PRE_FILE);
    	exit(15);
    }

    /*
     * if -d of -m num, perform multi digest processing instead
     */
    if (d_flag || m_flag) {
	multiMain(argc, argv, pre_str, pre_str_len, data_str, C_flag,
		  (UINT)multi);

    /*
     * if no -d and no -m num, process string, stdin or files
     */
    } else {

	/*
	 * case: digest a string
	 */
	if (data_str != NULL) {
	    data_str_len = strlen(data_str);
	    shs1Init(&digest);
	    shs1Update(&digest, pre_str, pre_str_len);
	    shs1Update(&digest, (BYTE *)data_str, data_str_len);
	    shs1Output(data_str, 1, c_flag, i_flag, &digest);

	/*
	 * case: digest stdin
	 */
	} else if (optind == argc) {
	    shs1Init(&digest);
	    shs1Stream(pre_str, pre_str_len, stdin, &digest);
	    shs1Output(NULL, 0, c_flag, i_flag, &digest);

	/*
	 * case: digest files
	 */
	} else {
	    for (; optind < argc; optind++) {
		shs1Init(&digest);
		shs1File(pre_str, pre_str_len, argv[optind], i_flag, &digest);
		shs1Output(argv[optind], 0, c_flag, i_flag, &digest);
	    }
	}
    }

    /* all done */
    /* exit(0); */
    return 0;
}
