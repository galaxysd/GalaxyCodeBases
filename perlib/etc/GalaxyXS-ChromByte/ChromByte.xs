#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"


MODULE = GalaxyXS::ChromByte		PACKAGE = GalaxyXS::ChromByte		

long
initchr(len)
	int len
	CODE:
		void *address = malloc( len + 1 );
		memset( address , 0 , len + 1 );
		RETVAL = (long)address;
	OUTPUT:
		RETVAL

void		
setbases( address, begin, end, val )
	long address
	int begin
	int end
	unsigned int val
	CODE:
		char * buf = ( char * ) address ;
		memset( buf + begin , val , end - begin + 1 );

unsigned int
orbase( address, pos, val )
	long address
	int pos
	unsigned int val
	CODE:
		char * buf = ( char * ) address ;
		RETVAL = *( buf + pos );
		RETVAL |= val;
		memset( buf + pos , val , 1 );
	OUTPUT:
		RETVAL

unsigned int
getbase( address, pos )
	 long address
	 int pos
	CODE:
		char * buf = ( char * ) address ;
		RETVAL = *( buf + pos );
	OUTPUT:
		RETVAL

void		
setbasesc( address, pos, val )
	long address
	int pos
	char val
	CODE:
		char * buf = ( char * ) address ;
		memset( buf + pos , val , 1 );

char
getbasec( address, pos )
	 long address
	 int pos
	CODE:
		char * buf = ( char * ) address ;
		RETVAL = *( buf + pos );
	OUTPUT:
		RETVAL

void
freechr( address )
	long address
	CODE:
		void * buf = ( void * ) address ;
		free( buf );
