/*
 * shs1 - new NIST Secure Hash Standard-1 (SHS1)
 *
 * @(#) $Revision: 13.6 $
 * @(#) $Id: shs1.h,v 13.6 2010/10/12 21:10:17 chongo Exp $
 * @(#) $Source: /usr/local/src/cmd/hash/RCS/shs1.h,v $
 *
 * Written 2 September 1992, Peter C. Gutmann.
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

#if !defined(SHS1_H)
#define SHS1_H

#include <sys/types.h>
#include <sys/stat.h>

/*
 * our version
 */
#define SHS1_VERSION "$Revision: 13.6 $"

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

/* SHS1_CHUNKSIZE must be a power of 2 - fixed value defined by the algorithm */
#define SHS1_CHUNKSIZE (1<<6)
#define SHS1_CHUNKMASK (SHS1_CHUNKSIZE-1)
#define SHS1_CHUNKWORDS (SHS1_CHUNKSIZE/sizeof(ULONG))

/* SHS1_DIGESTSIZE is a the length of the digest as defined by the algorithm */
#define SHS1_DIGESTSIZE (20)
#define SHS1_DIGESTWORDS (SHS1_DIGESTSIZE/sizeof(ULONG))

/* SHS_HIGH - the ULONG where the 64 bit count of bits processed is stored */
#define SHS1_HIGH (SHS1_CHUNKWORDS-2)

/* SHS1_BLOCKSIZE is how large a chunk multiStream processes at a time */
#define SHS1_BLOCKSIZE (SHS1_CHUNKSIZE<<8)

/* SHS1_READSIZE must be a multiple of SHS1_BLOCKSIZE */
#define SHS1_READSIZE (SHS1_BLOCKSIZE<<2)
#define SHS1_READWORDS (SHS1_READSIZE/sizeof(ULONG))

/* maximum part of pre_file used */
#define SHS1_MAX_PRE_FILE 32768

/*
 * SHS1_ROUNDUP(x,y) - round x up to the next multiple of y
 */
#define SHS1_ROUNDUP(x,y) ((((x)+(y)-1)/(y))*(y))

/*
 * SHS1COUNT(SHS1_INFO*, ULONG) - update the 64 bit count in an SHS1_INFO
 *
 * We will count bytes and convert to bit count during the final
 * transform.
 */
#define SHS1COUNT(shs1info, count) {				\
    ULONG tmp_countLo;						\
    tmp_countLo = (shs1info)->countLo;				\
    if (((shs1info)->countLo += (count)) < tmp_countLo) {	\
	(shs1info)->countHi++;					\
    }								\
}

/*
 * The structure for storing SHS1 info
 *
 * We will assume that bit count is a multiple of 8.
 */
typedef struct {
    ULONG digest[SHS1_DIGESTWORDS];	/* message digest */
    ULLONG octets;			/* count of octets processed */
    ULONG datalen;			/* length of data in data */
    ULONG data[SHS1_CHUNKWORDS];	/* SHS1 chunk buffer */
} SHS1_INFO;

/*
 * elements of the stat structure that we will process
 */
struct shs1_stat {
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

/* shs1.c */
void shs1Init(SHS1_INFO*);
void shs1Transform(ULONG *digest, ULONG *W);
void shs1Update(SHS1_INFO*, BYTE*, ULONG);
void shs1fullUpdate(SHS1_INFO*, BYTE*, ULONG);
void shs1Final(SHS1_INFO*);
extern char *shs1_what;

/*
 * Some external programs use the functions found in shs1.c and shs1io.c.
 * These routines define SHS1_IO so that such external programs will not
 * pick up the following declarations.
 */

#if !defined(SHS1_IO)

/* shs1io.c */
void shs1Stream(BYTE*, UINT, FILE*, SHS1_INFO*);
void shs1File(BYTE*, UINT, char*, int, SHS1_INFO*);
extern ULONG shs1_zero[];

/* shs1dual.c */
void multiMain(int, char**, BYTE*, UINT, char*, int, UINT);
void multiTest(void);
extern char *shs1dual_what;

/* shs1drvr.c */
void shs1Print(int, int, SHS1_INFO*);
extern int c_flag;
extern int i_flag;
extern int q_flag;
extern int dot_zero;
#endif /* SHS1_IO */
extern int debug;
extern char *program;

#endif /* SHS1_H */
