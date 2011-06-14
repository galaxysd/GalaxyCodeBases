#include <locale.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>

// http://c-faq.com/stdio/commaprint.html
char *commaprint(unsigned long n)
{
	static int comma = '\0';
	static char retbuf[4*(sizeof(unsigned long)*CHAR_BIT+2)/3/3+1];
	char *p = &retbuf[sizeof(retbuf)-1];
	int i = 0;

	if(comma == '\0') {
		struct lconv *lcp = localeconv();
		if(lcp != NULL) {
			if(lcp->thousands_sep != NULL &&
				*lcp->thousands_sep != '\0')
				comma = *lcp->thousands_sep;
			else	comma = ',';
		}
	}

	*p = '\0';

	do {
		if(i%3 == 0 && i != 0)
			*--p = comma;
		*--p = '0' + n % 10;
		n /= 10;
		i++;
	} while(n != 0);

    //return p;
    char *pret = malloc(&retbuf[sizeof(retbuf)]-p); // well, just a tiny leak ...
    return strcpy(pret,p);
}

