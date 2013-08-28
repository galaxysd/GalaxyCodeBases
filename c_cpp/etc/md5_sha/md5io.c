/*
 * md5drvr - md5 driver code
 *
 * @(#) $Revision: 13.1 $
 * @(#) $Id: md5io.c,v 13.1 2006/08/14 03:16:33 chongo Exp $
 * @(#) $Source: /usr/local/src/cmd/hash/RCS/md5io.c,v $
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
 * See md5drvr.c for version and modification history.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>
#define MD5_IO
#include "md5.h"

/* global variables */
ULONG md5_zero[MD5_MAXBLOCK/sizeof(ULONG)];	/* block of zeros */


/*
 * MD5Stream - digest a open file stream
 */
void
MD5Stream(pre_str, pre_len, stream, dig)
    BYTE *pre_str;		/* data prefix or NULL */
    UINT pre_len;		/* length of pre_str */
    FILE *stream;		/* the stream to process */
    MD5_CTX *dig;		/* current digest */
{
    ULONG data[MD5_READWORDS];	/* our read buffer */
    int bytes;			/* bytes last read */
    int ret;			/* partial fread return value */

    /*
     * pre-process prefix if needed
     */
    if (pre_str != NULL) {
	MD5Update(dig, pre_str, pre_len);
	MD5COUNT(dig, pre_len);
    }

    /*
     * if we have a partial chunk, try to read until we have a full chunk
     */
    clearerr(stream);
    if (dig->datalen > 0) {

        /* determine what we have so far */
        bytes = dig->datalen;

	/* try to read what we need to fill the chunk */
	while (bytes < MD5_CHUNKSIZE) {

	    /* try to read what we need */
	    ret = fread((char*)data+bytes, 1, MD5_CHUNKSIZE-bytes, stream);

	    /* carefully examine the result */
	    if (ret < 0 || ferror(stream)) {
		/* error processing */
		fprintf(stderr, "%s: ", program);
		perror("read #3 error");
		exit(1);
	    } else if (ret == 0 || feof(stream)) {
		/* EOF processing */
		MD5COUNT(dig, MD5_CHUNKSIZE-dig->datalen);
		MD5Update(dig, (BYTE *)data+dig->datalen,
		  MD5_CHUNKSIZE-dig->datalen);
		return;
	    }

	    /* note that we have more bytes */
	    bytes += ret;
        }
        MD5COUNT(dig, MD5_CHUNKSIZE-dig->datalen);
        MD5Update(dig, (BYTE *)data+dig->datalen, MD5_CHUNKSIZE-dig->datalen);
    }

    /*
     * process the contents of the file
     */
    while ((bytes = fread((char *)data, 1, MD5_READSIZE, stream)) > 0) {

	/*
	 * if we got a partial read, try to read up to a full chunk
	 */
	while (bytes < MD5_READSIZE) {

	    /* try to read more */
	    ret = fread((char *)data+bytes, 1, MD5_READSIZE-bytes, stream);

	    /* carefully examine the result */
	    if (ret < 0 || ferror(stream)) {
	    	/* error processing */
	    	fprintf(stderr, "%s: ", program);
	    	perror("read #1 error");
	    	exit(2);
	    } else if (ret == 0 || feof(stream)) {
	    	/* EOF processing */
	    	MD5Update(dig, (BYTE *)data, bytes);
	    	MD5COUNT(dig, bytes);
	    	return;
	    }

	    /* note that we have more bytes */
	    bytes += ret;
	}

	/*
	 * digest the read
	 */
	MD5fullUpdate(dig, (BYTE *)data, bytes);
	MD5COUNT(dig, bytes);
    }

    /*
     * watch for errors
     */
    if (bytes < 0 || ferror(stream)) {
	/* error processing */
	fprintf(stderr, "%s: ", program);
	perror("read #2 error");
	exit(3);
    }
    return;
}


/*
 * MD5File - digest a file
 */
void
MD5File(pre_str, pre_len, filename, inode, dig)
    BYTE *pre_str;		/* string prefix or NULL */
    UINT pre_len;		/* length of pre_str */
    char *filename;		/* the filename to process */
    int inode;			/*  1 => process inode & filename */
    MD5_CTX *dig;		/* current digest */
{
    FILE *inFile;		/* the open file stream */
    struct stat buf;		/* stat or lstat of file */
    struct md5_stat hashbuf;	/* stat data to digest */
    struct md5_stat hashlbuf;	/* lstat data to digest */
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
	    MD5Update(dig, (BYTE *)filename, filename_len);
	    MD5COUNT(dig, filename_len);
	}
    } else {
	if (inode) {
	    MD5Update(dig, pre_str, pre_len);
	    filename_len = strlen(filename);
	    MD5Update(dig, (BYTE *)filename, filename_len);
	    MD5COUNT(dig, pre_len + filename_len);
	} else {
	    MD5Update(dig, pre_str, pre_len);
	    MD5COUNT(dig, pre_len);
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
	MD5Update(dig, (BYTE *)&hashbuf, sizeof(hashbuf));
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
	MD5Update(dig, (BYTE *)&hashlbuf, sizeof(hashlbuf));

	/*
	 * pad with zeros to process file data faster
	 */
	if (dig->datalen > 0) {
	    MD5COUNT(dig, sizeof(hashbuf) + sizeof(hashlbuf) +
	          MD5_CHUNKSIZE - dig->datalen);
	    MD5Update(dig, (BYTE *)md5_zero, MD5_CHUNKSIZE - dig->datalen);
	} else {
	    MD5COUNT(dig, sizeof(hashbuf) + sizeof(hashlbuf));
	}
    }

    /*
     * process the data stream
     */
    MD5Stream(NULL, 0, inFile, dig);
    fclose(inFile);
}
