/*
 * mega - output 1,000,000 repetitions of "a"
 *
 * @(#) $Revision: 13.1 $
 * @(#) $Id: mega.c,v 13.1 2006/08/14 03:16:33 chongo Exp $
 * @(#) $Source: /usr/local/src/cmd/hash/RCS/mega.c,v $
 *
 * This file was written by Landon Curt Noll.
 *
 * This output is the message specified by FIPS PUB 180-1  APPENDIX C.
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
 */

#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <unistd.h>


/*ARGSUSED*/
int
main(int argc, char *argv[])
{
    int i;

    /* usage */
    if (argc != 1) {
	fprintf(stderr, "usage: %s\n", argv[0]);
	exit(1);
    }

    for (i=0; i < 1000000/50; ++i) {
	write(1, "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", 50);
    }
    exit(0);
}
