/*
 * shsfile - old Secure Hash Standard
 *
 * @(#) $Revision: 13.3 $
 * @(#) $Id: shsio.c,v 13.3 2006/08/14 11:25:24 chongo Exp $
 * @(#) $Source: /usr/local/src/cmd/hash/RCS/shsio.c,v $
 *
 **************************************************************************
 * This version implements the old Secure Hash Algorithm specified by     *
 * (FIPS Pub 180).  This version is kept for backward compatibility with  *
 * shs version 2.10.1.  See the shs utility for the new standard.	  *
 **************************************************************************
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#define SHS_IO
#include "shs.h"

/* global variables */
ULONG shs_zero[SHS_CHUNKWORDS];	/* block of zeros */


/*
 * shsStream - digest a open file stream
 *
 * given:
 *	pre_str		data prefix or NULL
 *	pre_len		length of pre_str
 *	stream		the stream to process
 *	dig		current digest
 */
void
shsStream(BYTE *pre_str, UINT pre_len, FILE *stream, SHS_INFO *dig)
{
    BYTE data[SHS_READSIZE];	/* our read buffer */
    unsigned int bytes;		/* bytes last read */
    int ret;			/* partial fread return value */

    /*
     * pre-process prefix if needed
     */
    if (pre_str != NULL && pre_len > 0) {
	shsUpdate(dig, pre_str, pre_len);
    }

    /*
     * Try to read as much as we can (up to SHS_READSIZE bytes).  If there
     * is data in the hash buffer already, try to read so that we can
     * perform as many full hash operations as possible.
     */
    clearerr(stream);
    ret = 0;
    do {

        /* determine what we have so far */
        bytes = dig->datalen;

	/* try to read what we need to fill the chunk */
	while (!feof(stream) && bytes < SHS_READSIZE) {

	    /* try to read what we need */
#if defined(DEBUG)
	    if (debug) {
		fprintf(stderr, "DEBUG: have: %u  will try to read: %u\n",
		    bytes, SHS_READSIZE-bytes);
	    }
#endif /* DEBUG */
	    ret = fread(data+bytes, 1, SHS_READSIZE-bytes, stream);
#if defined(DEBUG)
	    if (debug) {
		fprintf(stderr, "DEBUG: fread returrned: %d\n", ret);
	    }
#endif /* DEBUG */

	    /* look read for errors */
	    if (ret < 0 || ferror(stream)) {
		/* error processing */
		fprintf(stderr, "%s: ", program);
		perror("shs read error");
		exit(1);
	    }

	    /* note that we have more bytes */
	    bytes += ret;
        }

	/* process whatever new data we have in the buffer */
	if (bytes > dig->datalen) {
	    shsUpdate(dig, data+dig->datalen, bytes-dig->datalen);
	}
    } while (ret > 0 && !feof(stream));
    return;
}


/*
 * shsFile - digest a file
 *
 * given:
 *	pre_str		string prefix or NULL
 *	pre_len		length of pre_str
 *	filename	the filename to process
 *	inode		1 => process inode & filename
 *	dig		current digest
 */
void
shsFile(BYTE *pre_str, UINT pre_len, char *filename,
    int inode, SHS_INFO *dig)
{
    FILE *inFile;		/* the open file stream */
    struct stat buf;		/* stat or lstat of file */
    struct shs_stat hashbuf;	/* stat data to digest */
    struct shs_stat hashlbuf;	/* lstat data to digest */
    ULONG filename_len;		/* length of the filename */

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
    if (pre_str == NULL) {
	if (inode) {
	    filename_len = strlen(filename);
	    shsUpdate(dig, (BYTE *)filename, filename_len);
#if defined(DEBUG)
	    if (debug) {
		fprintf(stderr,
		    "DEBUG: filename_len:%lu octets:%llu\n",
		    (unsigned long)filename_len, dig->octets);
	    }
#endif
	}
    } else {
	if (inode) {
	    shsUpdate(dig, pre_str, pre_len);
	    filename_len = strlen(filename);
	    shsUpdate(dig, (BYTE *)filename, filename_len);
#if defined(DEBUG)
	    if (debug) {
		fprintf(stderr,
		    "DEBUG: pre_len:%lu filename_len:%lu octets:%lld\n",
		    (unsigned long)pre_len,
		    (unsigned long)filename_len, dig->octets);
	    }
#endif
	} else {
	    shsUpdate(dig, pre_str, pre_len);
	}
    }

    /*
     * digest file stat and lstat
     */
    if (inode) {
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
#if defined(DEBUG)
	if (debug) {
	    fprintf(stderr,
	      "DEBUG: dev:%ld ino:%ld mode:%o nlink:%d uid:%d gid:%d\n",
	      (unsigned long)hashbuf.stat_dev,
	      (unsigned long)hashbuf.stat_ino,
	      hashbuf.stat_mode, hashbuf.stat_nlink,
	      hashbuf.stat_uid, hashbuf.stat_gid);
	    fprintf(stderr,
	      "DEBUG: size:%llu mtime:%llu ctime:%llu\n",
	      (unsigned long long)hashbuf.stat_size,
	      (unsigned long long)hashbuf.stat_mtime,
	      (unsigned long long)hashbuf.stat_ctime);
	}
#endif
	shsUpdate(dig, (BYTE *)&hashbuf, sizeof(hashbuf));
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
#if defined(DEBUG)
	if (debug) {
	    fprintf(stderr,
	      "DEBUG: ldev:%ld lino:%ld mode:%o lnlink:%d luid:%d lgid:%d\n",
	      (unsigned long)hashbuf.stat_dev,
	      (unsigned long)hashbuf.stat_ino,
	      hashlbuf.stat_mode, hashlbuf.stat_nlink,
	      hashlbuf.stat_uid, hashlbuf.stat_gid);
	    fprintf(stderr,
	      "DEBUG: lsize:%llu lmtime:%llu lctime:%llu\n",
	      (unsigned long long)hashlbuf.stat_size,
	      (unsigned long long)hashlbuf.stat_mtime,
	      (unsigned long long)hashbuf.stat_ctime);
	}
#endif
	shsUpdate(dig, (BYTE *)&hashlbuf, sizeof(hashlbuf));

	/*
	 * pad with zeros to process file data faster
	 */
	if (dig->datalen > 0) {
#if defined(DEBUG)
	    if (debug) {
		fprintf(stderr,
		  "DEBUG: pad_len:%lu\n",
		  (unsigned long)(SHS_CHUNKSIZE - dig->datalen));
	    }
#endif
	    shsUpdate(dig, (BYTE *)shs_zero, SHS_CHUNKSIZE - dig->datalen);
	}
#if defined(DEBUG)
	if (debug) {
	    fprintf(stderr, "DEBUG: datalen:%lu count:%llu\n",
	      (unsigned long)dig->datalen, dig->octets);
	}
#endif
    }

    /*
     * process the data stream
     */
    shsStream(NULL, 0, inFile, dig);
    fclose(inFile);
}
