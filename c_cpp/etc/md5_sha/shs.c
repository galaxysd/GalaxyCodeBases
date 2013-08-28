/* $Id: shs.c,v 13.5 2010/10/12 21:10:17 chongo Exp $  */
/*
 * shs - old Secure Hash Standard
 *
 * @(#) $Revision: 13.5 $
 * @(#) $Id: shs.c,v 13.5 2010/10/12 21:10:17 chongo Exp $
 * @(#) $Source: /usr/local/src/cmd/hash/RCS/shs.c,v $
 *
 **************************************************************************
 * This version implements the old Secure Hash Algorithm specified by     *
 * (FIPS Pub 180).  This version is kept for backward compatibility with  *
 * shs version 2.10.1.  See the shs utility for the new standard.	  *
 **************************************************************************
 *
 * Written 2 September 1992, Peter C. Gutmann.
 *
 * This file was Modified/Re-written by Landon Curt Noll.
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

#include <stdio.h>
#include <string.h>
#define SHS_IO
#include "shs.h"
#include "align.h"
#include "endian.h"

/*
 * The SHS f()-functions.  The f1 and f3 functions can be optimized
 * to save one boolean operation each - thanks to Rich Schroeppel,
 * rcs@cs.arizona.edu for discovering this.
 *
 * f1: ((x&y) | (~x&z)) == (z ^ (x&(y^z)))
 * f3: ((x&y) | (x&z) | (y&z)) == ((x&y) | (z&(x|y)))
 */
#define f1(x,y,z)       (z ^ (x&(y^z)))		/* Rounds  0-19 */
#define f2(x,y,z)       (x^y^z)			/* Rounds 20-39 */
#define f3(x,y,z)       ((x&y) | (z&(x|y)))	/* Rounds 40-59 */
#define f4(x,y,z)       (x^y^z)			/* Rounds 60-79 */

/* The SHS Mysterious Constants */
#define K1      0x5A827999L	/* Rounds  0-19 */
#define K2      0x6ED9EBA1L	/* Rounds 20-39 */
#define K3      0x8F1BBCDCL	/* Rounds 40-59 */
#define K4      0xCA62C1D6L	/* Rounds 60-79 */

/* SHS initial values */
#define h0init  0x67452301L
#define h1init  0xEFCDAB89L
#define h2init  0x98BADCFEL
#define h3init  0x10325476L
#define h4init  0xC3D2E1F0L

/* 32-bit rotate left - kludged with shifts */
#define LEFT_ROT(X,n)  (((X)<<(n)) | ((X)>>(32-(n))))

/*
 * The initial expanding function.  The hash function is defined over an
 * 80-word expanded input array W, where the first 16 are copies of the input
 * data, and the remaining 64 are defined by
 *
 *      W[i] = W[i-16] ^ W[i-14] ^ W[i-8] ^ W[i-3]
 *
 * This implementation generates these values on the fly in a circular
 * buffer - thanks to Colin Plumb (colin@nyx10.cs.du.edu) for this
 * optimization.
 */
#define exor(W,i) (W[i&15] ^= (W[(i-14)&15] ^ W[(i-8)&15] ^ W[(i-3)&15]))

/*
 * The prototype SHS sub-round.  The fundamental sub-round is:
 *
 *      a' = e + LEFT_ROT(a,5) + f(b,c,d) + k + data;
 *      b' = a;
 *      c' = LEFT_ROT(b,30);
 *      d' = c;
 *      e' = d;
 *
 * but this is implemented by unrolling the loop 5 times and renaming the
 * variables ( e, a, b, c, d ) = ( a', b', c', d', e' ) each iteration.
 * This code is then replicated 20 times for each of the 4 functions, using
 * the next 20 values from the W[] array each time.
 */
#define subRound(a, b, c, d, e, f, k, data) \
    (e += LEFT_ROT(a,5) + f(b,c,d) + k + data, b = LEFT_ROT(b,30))


/*
 * shsInit - initialize the SHS state
 */
void
shsInit(SHS_INFO *dig)
{
    /* Set the h-vars to their initial values */
    dig->digest[0] = h0init;
    dig->digest[1] = h1init;
    dig->digest[2] = h2init;
    dig->digest[3] = h3init;
    dig->digest[4] = h4init;

    /* Initialise bit count */
    dig->octets = 0;
    dig->datalen = 0;
}


/*
 * shsTransform - perform the SHS transformatio
 *
 * Note that this code, like MD5, seems to break some optimizing compilers.
 * It may be necessary to split it into sections, eg based on the four
 * subrounds.  One may also want to roll each subround into a loop.
 */
void
shsTransform(ULONG *digest, ULONG *W)
{
    ULONG A, B, C, D, E;	/* Local vars */
#if BYTE_ORDER == LITTLE_ENDIAN || defined(DETAILED_DEBUG)
    unsigned int i;
#endif /* BYTE_ORDER == LITTLE_ENDIAN || DETAILED_DEBUG */

#if defined(DETAILED_DEBUG)
    if (debug) {
	ULONG dig;	/* digest in Big Endian order */
	unsigned char *p;
	unsigned int j;

	/* print the input digest */
	fprintf(stderr, "DEBUG: input digest: ");
	for (i=0; i < SHS_DIGESTWORDS; ++i) {
#if BYTE_ORDER == LITTLE_ENDIAN
	    dig = ((digest[i]<<16) | (digest[i]>>16));
	    dig = (((digest[i] & 0xff00ff00UL) >> 8) |
	    	   ((digest[i] & 0x00ff00ffUL) >> 8));
#else /* BYTE_ORDER == LITTLE_ENDIAN */
	    dig = digest[i];
#endif /* BYTE_ORDER == LITTLE_ENDIAN */
	    for (j=0, p=(unsigned char *)&dig; j < sizeof(ULONG); ++j, ++p) {
		fprintf(stderr, "%02x  ", *p);
	    }
	}
	fputc('\n', stderr);

	/* print the input buffer */
	fprintf(stderr, "DEBUG: input buffer: ");
	for (p = (unsigned char *)W, i=0; i < SHS_CHUNKSIZE; ++i, ++p) {
	    fprintf(stderr, "%02x  ", *p);
	}
	fputc('\n', stderr);
    }
#endif /* DETAILED_DEBUG */

    /* Set up first buffer and local data buffer */
    A = digest[0];
    B = digest[1];
    C = digest[2];
    D = digest[3];
    E = digest[4];

    /*
     * The heavy mangling below was encoded for a Big Endian machine.
     * We must byte swap the data on a Little Endian machine unfortunately.
     */
#if BYTE_ORDER == LITTLE_ENDIAN
    for (i=0; i < SHS_CHUNKWORDS; ++i) {
	W[i] = ((W[i]<<16) | (W[i]>>16));
	W[i] = (((W[i] & 0xff00ff00UL) >> 8) |
		((W[i] & 0x00ff00ffUL) << 8));
    }
#endif /* BYTE_ORDER == LITTLE_ENDIAN */

    /* Heavy mangling, in 4 sub-rounds of 20 interations each. */
    subRound(A, B, C, D, E, f1, K1, W[ 0]);
    subRound(E, A, B, C, D, f1, K1, W[ 1]);
    subRound(D, E, A, B, C, f1, K1, W[ 2]);
    subRound(C, D, E, A, B, f1, K1, W[ 3]);
    subRound(B, C, D, E, A, f1, K1, W[ 4]);
    subRound(A, B, C, D, E, f1, K1, W[ 5]);
    subRound(E, A, B, C, D, f1, K1, W[ 6]);
    subRound(D, E, A, B, C, f1, K1, W[ 7]);
    subRound(C, D, E, A, B, f1, K1, W[ 8]);
    subRound(B, C, D, E, A, f1, K1, W[ 9]);
    subRound(A, B, C, D, E, f1, K1, W[10]);
    subRound(E, A, B, C, D, f1, K1, W[11]);
    subRound(D, E, A, B, C, f1, K1, W[12]);
    subRound(C, D, E, A, B, f1, K1, W[13]);
    subRound(B, C, D, E, A, f1, K1, W[14]);
    subRound(A, B, C, D, E, f1, K1, W[15]);
    subRound(E, A, B, C, D, f1, K1, exor(W,16));
    subRound(D, E, A, B, C, f1, K1, exor(W,17));
    subRound(C, D, E, A, B, f1, K1, exor(W,18));
    subRound(B, C, D, E, A, f1, K1, exor(W,19));

    subRound(A, B, C, D, E, f2, K2, exor(W,20));
    subRound(E, A, B, C, D, f2, K2, exor(W,21));
    subRound(D, E, A, B, C, f2, K2, exor(W,22));
    subRound(C, D, E, A, B, f2, K2, exor(W,23));
    subRound(B, C, D, E, A, f2, K2, exor(W,24));
    subRound(A, B, C, D, E, f2, K2, exor(W,25));
    subRound(E, A, B, C, D, f2, K2, exor(W,26));
    subRound(D, E, A, B, C, f2, K2, exor(W,27));
    subRound(C, D, E, A, B, f2, K2, exor(W,28));
    subRound(B, C, D, E, A, f2, K2, exor(W,29));
    subRound(A, B, C, D, E, f2, K2, exor(W,30));
    subRound(E, A, B, C, D, f2, K2, exor(W,31));
    subRound(D, E, A, B, C, f2, K2, exor(W,32));
    subRound(C, D, E, A, B, f2, K2, exor(W,33));
    subRound(B, C, D, E, A, f2, K2, exor(W,34));
    subRound(A, B, C, D, E, f2, K2, exor(W,35));
    subRound(E, A, B, C, D, f2, K2, exor(W,36));
    subRound(D, E, A, B, C, f2, K2, exor(W,37));
    subRound(C, D, E, A, B, f2, K2, exor(W,38));
    subRound(B, C, D, E, A, f2, K2, exor(W,39));

    subRound(A, B, C, D, E, f3, K3, exor(W,40));
    subRound(E, A, B, C, D, f3, K3, exor(W,41));
    subRound(D, E, A, B, C, f3, K3, exor(W,42));
    subRound(C, D, E, A, B, f3, K3, exor(W,43));
    subRound(B, C, D, E, A, f3, K3, exor(W,44));
    subRound(A, B, C, D, E, f3, K3, exor(W,45));
    subRound(E, A, B, C, D, f3, K3, exor(W,46));
    subRound(D, E, A, B, C, f3, K3, exor(W,47));
    subRound(C, D, E, A, B, f3, K3, exor(W,48));
    subRound(B, C, D, E, A, f3, K3, exor(W,49));
    subRound(A, B, C, D, E, f3, K3, exor(W,50));
    subRound(E, A, B, C, D, f3, K3, exor(W,51));
    subRound(D, E, A, B, C, f3, K3, exor(W,52));
    subRound(C, D, E, A, B, f3, K3, exor(W,53));
    subRound(B, C, D, E, A, f3, K3, exor(W,54));
    subRound(A, B, C, D, E, f3, K3, exor(W,55));
    subRound(E, A, B, C, D, f3, K3, exor(W,56));
    subRound(D, E, A, B, C, f3, K3, exor(W,57));
    subRound(C, D, E, A, B, f3, K3, exor(W,58));
    subRound(B, C, D, E, A, f3, K3, exor(W,59));

    subRound(A, B, C, D, E, f4, K4, exor(W,60));
    subRound(E, A, B, C, D, f4, K4, exor(W,61));
    subRound(D, E, A, B, C, f4, K4, exor(W,62));
    subRound(C, D, E, A, B, f4, K4, exor(W,63));
    subRound(B, C, D, E, A, f4, K4, exor(W,64));
    subRound(A, B, C, D, E, f4, K4, exor(W,65));
    subRound(E, A, B, C, D, f4, K4, exor(W,66));
    subRound(D, E, A, B, C, f4, K4, exor(W,67));
    subRound(C, D, E, A, B, f4, K4, exor(W,68));
    subRound(B, C, D, E, A, f4, K4, exor(W,69));
    subRound(A, B, C, D, E, f4, K4, exor(W,70));
    subRound(E, A, B, C, D, f4, K4, exor(W,71));
    subRound(D, E, A, B, C, f4, K4, exor(W,72));
    subRound(C, D, E, A, B, f4, K4, exor(W,73));
    subRound(B, C, D, E, A, f4, K4, exor(W,74));
    subRound(A, B, C, D, E, f4, K4, exor(W,75));
    subRound(E, A, B, C, D, f4, K4, exor(W,76));
    subRound(D, E, A, B, C, f4, K4, exor(W,77));
    subRound(C, D, E, A, B, f4, K4, exor(W,78));
    subRound(B, C, D, E, A, f4, K4, exor(W,79));

    /* Build message digest */
    digest[0] += A;
    digest[1] += B;
    digest[2] += C;
    digest[3] += D;
    digest[4] += E;

#if defined(DETAILED_DEBUG)
    if (debug) {
	int i;
	unsigned char *p;

	/* print the final digest */
	fprintf(stderr, "DEBUG: final digest: ");
	for (p = (unsigned char *)digest, i=0; i < SHS_DIGESTSIZE; ++i, ++p) {
	    fprintf(stderr, "%02x  ", *p);
	}
	fputc('\n', stderr);
    }
#endif /* DETAILED_DEBUG */
}


/*
 * shsUpdate - update SHS with arbitrary length data
 *
 * This code does not assume that the buffer size is a multiple of
 * SHS_CHUNKSIZE bytes long.  This code handles partial chunk between
 * calls to shsUpdate().
 */
void
shsUpdate(SHS_INFO *dig, BYTE *buffer, ULONG count)
{
    ULONG datalen = dig->datalen;

    /*
     * Catch the case of a non-empty data buffer
     */
    if (datalen > 0) {

	/* determine the size we need to copy to fill the buffer */
	ULONG len_to_fill = SHS_CHUNKSIZE - datalen;

	/* case: new data will not fill the buffer */
	if (len_to_fill > count) {
	    memcpy(((BYTE *)dig->data)+datalen, buffer, count);
	    dig->datalen = datalen+count;
	    return;

	/* case: buffer will be filled */
	} else {
	    memcpy(((BYTE *)dig->data)+datalen, buffer, len_to_fill);
	    shsTransform(dig->digest, dig->data);
	    buffer += len_to_fill;
	    count -= len_to_fill;
	    dig->octets += SHS_CHUNKSIZE;
	    dig->datalen = 0;
	}
    }

    /*
     * Process data in SHS_CHUNKSIZE chunks
     */
    if (count >= SHS_CHUNKSIZE) {
	shsfullUpdate(dig, buffer, (count & ~SHS_CHUNKMASK));
	buffer += (count & ~SHS_CHUNKMASK);
	count &= SHS_CHUNKMASK;
    }

    /*
     * Handle any remaining bytes of data.
     * This should only happen once on the final lot of data
     */
    if (count > 0) {
	memcpy((char *)dig->data, (char *)buffer, count);
    }
    dig->datalen = count;
}


/*
 * shsfullUpdate - update SHS with chunk multiple length data
 *
 * This function assumes that count is a multiple of SHS_CHUNKSIZE and that
 * no partial chunk is left over from a previous call.
 */
void
shsfullUpdate(SHS_INFO *dig, BYTE *buffer, ULONG count)
{
#if defined(MUST_ALIGN)
    ULONG in[SHS_CHUNKWORDS];	/* aligned buffer */
#endif /* MUST_ALIGN */

    /*
     * Process data in SHS_CHUNKSIZE chunks
     */
    while (count >= SHS_CHUNKSIZE) {
#if defined(MUST_ALIGN)
	if ((long)buffer & (sizeof(ULONG)-1)) {
	    memcpy((char *)in, (char *)buffer, SHS_CHUNKSIZE);
	    shsTransform(dig->digest, in);
	} else {
	    shsTransform(dig->digest, buffer);
	}
#else /* MUST_ALIGN */
	shsTransform(dig->digest, (ULONG *)buffer);
#endif /* MUST_ALIGN */
	buffer += SHS_CHUNKSIZE;
	count -= SHS_CHUNKSIZE;
	dig->octets += SHS_CHUNKSIZE;
    }
}


/*
 * shsFinal - perform final SHS transforms
 *
 * At this point we have less than a full chunk of data remaining
 * (and possibly no data) in the shs state data buffer.
 *
 * First we append a final 0x80 byte.
 *
 * Next if we have more than 56 bytes, we will zero fill the remainder
 * of the chunk, transform and then zero fill the first 56 bytes.
 * If we have 56 or fewer bytes, we will zero fill out to the 56th
 * chunk byte.  Regardless, we wind up with 56 bytes data.
 *
 * Finally we append the 64 bit length on to the 56 bytes of data
 * remaining.  This final chunk is transformed.
 */
void
shsFinal(SHS_INFO *dig)
{
    int count = dig->datalen;	/* count of actual data in buffer */
    ULLONG bits;		/* number of bits of data processed */

    /*
     * Set the first char of padding to 0x80.
     * This is safe since there is always at least one byte free
     */
    dig->octets += count;	/* add in any remaining data in buffer */
    ((BYTE *)dig->data)[count++] = 0x80;	/* add in guard octet */

    /* Pad out to 56 mod SHS_CHUNKSIZE */
    if (count > 56) {
	/* Two lots of padding:  Pad the first chunk to SHS_CHUNKSIZE bytes */
	memset((BYTE *)dig->data + count, 0, SHS_CHUNKSIZE - count);
	shsTransform(dig->digest, dig->data);

	/* Now fill the next chunk with 56 bytes */
	memset(dig->data, 0, 56);
    } else {
	/* Pad chunk to 56 bytes */
	memset((BYTE *)dig->data + count, 0, 56 - count);
    }

    /*
     * Append length in bits and transform
     *
     * We assume that bit count is a multiple of 8 because we have
     * only processed full bytes.
     */
    bits = (dig->octets << 3);
#if BYTE_ORDER == LITTLE_ENDIAN
    bits = ((bits << 32) | (bits >> 32));
    bits = ((bits & 0xffff0000ffff0000ULL) >> 16) |
	   ((bits & 0x0000ffff0000ffffULL) << 16);
    bits = ((bits & 0xff00ff00ff00ff00ULL) >> 8) |
	   ((bits & 0x00ff00ff00ff00ffULL) << 8);
#endif /* BYTE_ORDER == LITTLE_ENDIAN */
    memcpy(&(dig->data[SHS_HIGH]), &bits, sizeof(ULLONG));
    shsTransform(dig->digest, dig->data);
    dig->datalen = 0;
}
