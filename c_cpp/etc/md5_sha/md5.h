/*
 * md5 - RSA Data Security, Inc. MD5 Message-Digest Algorithm
 *
 * @(#) $Revision: 13.5 $
 * @(#) $Id: md5.h,v 13.5 2010/10/12 21:10:17 chongo Exp $
 * @(#) $Source: /usr/local/src/cmd/hash/RCS/md5.h,v $
 *
 * This file was written by RSA Data Security.
 *
 * This file was modified by Landon Curt Noll.
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
 * See md5drvr.c for version and modification history.
 */

/*
 ***********************************************************************
 ** Copyright (C) 1990, RSA Data Security, Inc. All rights reserved.  **
 **								      **
 ** License to copy and use this software is granted provided that    **
 ** it is identified as the "RSA Data Security, Inc. MD5 Message-     **
 ** Digest Algorithm" in all material mentioning or referencing this  **
 ** software or this function.					      **
 **								      **
 ** License is also granted to make and use derivative works	      **
 ** provided that such works are identified as "derived from the RSA  **
 ** Data Security, Inc. MD5 Message-Digest Algorithm" in all          **
 ** material mentioning or referencing the derived work.	      **
 **								      **
 ** RSA Data Security, Inc. makes no representations concerning       **
 ** either the merchantability of this software or the suitability    **
 ** of this software for any particular purpose.  It is provided "as  **
 ** is" without express or implied warranty of any kind.              **
 **								      **
 ** These notices must be retained in any copies of any part of this  **
 ** documentation and/or software.				      **
 ***********************************************************************
 */

#if !defined(MD5_H)
#define MD5_H

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

/*
 * our version
 */
#define MD5_VERSION "$Revision: 13.5 $"

/*
 * These macros are in common with shs.h, shs1.h and md5.h.  We use
 * HASH_MACROS to gaurd against multiple inclusion by external progs
 * that may want to use multiple hash codes in one module.
 */
#if !defined(HASH_MACROS)
#define HASH_MACROS

/*
 * Useful defines/typedefs
 */
typedef u_int8_t BYTE;		/* must be 1 byte unsigned value */
typedef u_int16_t UINT;		/* must be 2 byte unsigned value */
typedef u_int32_t ULONG;	/* must be 4 byte unsigned value */
typedef u_int64_t ULLONG;	/* must be 8 byte unsigned value */

#endif /* HASH_MACROS */

/* MD5_CHUNKSIZE must be a power of 2 - fixed value defined by the algorithm */
#define MD5_CHUNKSIZE	(1<<6)
#define MD5_CHUNKWORDS (MD5_CHUNKSIZE/sizeof(ULONG))

/* MD5_DIGESTSIZE is a the length of the digest as defined by the algorithm */
#define MD5_DIGESTSIZE	(16)
#define MD5_DIGESTWORDS (MD5_DIGESTSIZE/sizeof(ULONG))

/* MD5_LOW - where low 32 bits of 64 bit count is stored during final */
#define MD5_LOW 14

/* MD5_HIGH - where high 32 bits of 64 bit count is stored during final */
#define MD5_HIGH 15

/*
 * MD5_MAXBLOCK - maximum blocking factor
 *
 * must be power of 2 > MD5_CHUNKSIZE and < MD5_READSIZE and < 2^29
 */
#define MD5_MAXBLOCK (MD5_CHUNKSIZE<<7)

/* MD5_READSIZE must be a multiple of MD5_CHUNKSIZE >= MD5_MAXBLOCK */
#define MD5_READSIZE (MD5_CHUNKSIZE<<9)
#define MD5_READWORDS (MD5_READSIZE/sizeof(ULONG))

/* maximum size of pre_file that is used <= MD5_MAXBLOCK */
#define MD5_MAX_PRE_FILE MD5_MAXBLOCK

/*
 * MD5COUNT(MD5_CTX*, ULONG) - update the 64 bit count in an MD5_CTX
 *
 * We will count bytes and convert to bit count during the final
 * transform.
 */
#define MD5COUNT(md5info, count) {				\
    unsigned long tmp_countLo;					\
    tmp_countLo = (md5info)->countLo;				\
    if (((md5info)->countLo += (count)) < tmp_countLo) {	\
	(md5info)->countHi++;					\
    }								\
}

/*
 * MD5_ROUNDUP(x,y) - round x up to the next multiple of y
 */
#define MD5_ROUNDUP(x,y) ((((x)+(y)-1)/(y))*(y))

/*
 * Data structure for MD5 (Message-Digest) computation
 */
typedef struct {
    BYTE digest[MD5_DIGESTSIZE];	/* actual digest after MD5Final call */
    ULONG countLo;		/* 64 bit count: bits 3-34 */
    ULONG countHi;		/* 64 bit count: bits 35-63 (64-66 ignored) */
    ULONG datalen;		/* length of data in inp.inp_BYTE */
    ULONG sub_block;		/* length of current partial block or 0 */
    union {
	BYTE inp_BYTE[MD5_CHUNKSIZE];	  /* BYTE chunk buffer */
	ULONG inp_ULONG[MD5_CHUNKWORDS];  /* ULONG chunk buffer */
    } inp;
    ULONG buf[MD5_DIGESTWORDS];	       /* scratch buffer */
} MD5_CTX;
#define in_byte inp.inp_BYTE
#define in_ulong inp.inp_ULONG

/*
 * elements of the stat structure that we will process
 */
struct md5_stat {
    dev_t stat_dev;
    ino_t stat_ino;
    mode_t stat_mode;
    nlink_t stat_nlink;
    uid_t stat_uid;
    gid_t stat_gid;
    off_t stat_size;
    time_t stat_mtime;
    time_t stat_ctime;
};

/*
 * Turn off prototypes if requested
 */
#if (defined(NOPROTO) && defined(PROTO))
# undef PROTO
#endif

/* md5.c */
void MD5Init(MD5_CTX*);
void MD5Update(MD5_CTX*, BYTE*, UINT);
void MD5fullUpdate(MD5_CTX*, BYTE*, UINT);
void MD5Final(MD5_CTX*);
void MD5Transform(ULONG *buf, ULONG *in);
extern char *MD5_what;

/* md5io.c */
void MD5Stream(BYTE*, UINT, FILE*, MD5_CTX*);
void MD5File(BYTE*, UINT, char*, int, MD5_CTX*);
extern ULONG md5_zero[];

/*
 * Some external programs use the functions found in md5.c and md5io.c.
 * These routines define MD5_IO so that such external programs will not
 * pick up the following declarations.
 */

#if !defined(MD5_IO)

/* md5dual.c */
void dualMain(int, char**, BYTE*, UINT, char*);
void dualTest(void);
extern char *MD5dual_what;

/* md5drvr.c */
void MD5Print(MD5_CTX*);
extern char *program;
extern int i_flag;
extern int q_flag;
extern int dot_zero;
#endif /* MD5_IO */
extern int debug;
extern char *program;

#endif /* MD5_H */
