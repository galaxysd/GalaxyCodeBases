/*
 * shs - old Secure Hash Standard
 *
 * @(#) $Revision: 13.6 $
 * @(#) $Id: shs.h,v 13.6 2010/10/12 21:10:17 chongo Exp $
 * @(#) $Source: /usr/local/src/cmd/hash/RCS/shs.h,v $
 *
 **************************************************************************
 * This version implements the old Secure Hash Algorithm specified by     *
 * (FIPS Pub 180).  This version is kept for backward compatibility with  *
 * shs version 2.10.1.  See the shs utility for the new standard.         *
 **************************************************************************
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
 * See shsdrvr.c for version and modification history.
 */

#if !defined(SHS_H)
#define SHS_H

#include <sys/types.h>
#include <sys/stat.h>

/*
 * our version
 */
#define SHS_VERSION "$Revision: 13.6 $"

/*
 * These macros are in common with shs.h, shs.h and md5.h.  We use
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

/* SHS_CHUNKSIZE must be a power of 2 - fixed value defined by the algorithm */
#define SHS_CHUNKSIZE (1<<6)
#define SHS_CHUNKMASK (SHS_CHUNKSIZE-1)
#define SHS_CHUNKWORDS (SHS_CHUNKSIZE/sizeof(ULONG))

/* SHS_DIGESTSIZE is a the length of the digest as defined by the algorithm */
#define SHS_DIGESTSIZE (20)
#define SHS_DIGESTWORDS (SHS_DIGESTSIZE/sizeof(ULONG))

/* SHS_HIGH - the ULONG where the 64 bit count of bits processed is stored */
#define SHS_HIGH (SHS_CHUNKWORDS-2)

/* SHS_BLOCKSIZE is how large a chunk multiStream processes at a time */
#define SHS_BLOCKSIZE (SHS_CHUNKSIZE<<8)

/* SHS_READSIZE must be a multiple of SHS_BLOCKSIZE */
#define SHS_READSIZE (SHS_BLOCKSIZE<<2)
#define SHS_READWORDS (SHS_READSIZE/sizeof(ULONG))

/* maximum part of pre_file used */
#define SHS_MAX_PRE_FILE 32768

/*
 * SHS_ROUNDUP(x,y) - round x up to the next multiple of y
 */
#define SHS_ROUNDUP(x,y) ((((x)+(y)-1)/(y))*(y))

/*
 * The structure for storing SHS info
 *
 * We will assume that bit count is a multiple of 8.
 */
typedef struct {
    ULONG digest[SHS_DIGESTWORDS];	/* message digest */
    ULLONG octets;			/* count of octets processed */
    ULONG datalen;			/* length of data in data */
    ULONG data[SHS_CHUNKWORDS];		/* SHS chunk buffer */
} SHS_INFO;

/*
 * elements of the stat structure that we will process
 */
struct shs_stat {
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

/* shs.c */
void shsInit(SHS_INFO*);
void shsTransform(ULONG *digest, ULONG *W);
void shsUpdate(SHS_INFO*, BYTE*, ULONG);
void shsfullUpdate(SHS_INFO*, BYTE*, ULONG);
void shsFinal(SHS_INFO*);
extern char *shs_what;

/*
 * Some external programs use the functions found in shs.c and shsio.c.
 * These routines define SHS_IO so that such external programs will not
 * pick up the following declarations.
 */

#if !defined(SHS_IO)

/* shsio.c */
void shsStream(BYTE*, UINT, FILE*, SHS_INFO*);
void shsFile(BYTE*, UINT, char*, int, SHS_INFO*);
extern ULONG shs_zero[];

/* shsdual.c */
void multiMain(int, char**, BYTE*, UINT, char*, int, UINT);
void multiTest(void);
extern char *shsdual_what;

/* shsdrvr.c */
void shsPrint(int, int, SHS_INFO*);
extern int c_flag;
extern int i_flag;
extern int q_flag;
extern int dot_zero;
#endif /* SHS_IO */
extern int debug;
extern char *program;

#endif /* SHS_H */
