/*
 * must_align - determine if longs must be aligned
 *
 * @(#) $Revision: 13.1 $
 * @(#) $Id: must_align.c,v 13.1 2006/08/14 03:16:33 chongo Exp $
 * @(#) $Source: /usr/local/src/cmd/hash/RCS/must_align.c,v $
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
 */

#include <stdio.h>
#include <signal.h>
#include <stdlib.h>

void buserr(int);		/* catch alignment errors */

/*ARGSUSED*/
int
main(int argc, char *argv[])
{
    char byte[2*sizeof(unsigned long)];	/* mis-alignment buffer */
    unsigned long *p;	/* mis-alignment pointer */
    int i;

    /* usage */
    if (argc != 1) {
	fprintf(stderr, "usage: %s\n", argv[0]);
	exit(1);
    }

#if !defined(MUST_ALIGN)
    /* setup to catch alignment bus errors */
    signal(SIGBUS, buserr);
    signal(SIGSEGV, buserr);	/* some systems will generate SEGV instead! */

    /* mis-align our long fetches */
    for (i=0; i < sizeof(long); ++i) {
	p = (unsigned long *)(byte+i);
	*((volatile unsigned long *) p) = i;
	*((volatile unsigned long *) p) += 1;
    }

    /* if we got here, then we can mis-align longs */
    printf("#undef MUST_ALIGN\n");

#else
    /* force alignment */
    printf("#define MUST_ALIGN\n");
#endif
    exit(0);
}


/*
 * buserr - catch an alignment error
 */
/*ARGSUSED*/
void
buserr(int sig)
{
    /* alignment is required */
    printf("#define MUST_ALIGN\n");
    exit(0);
}
