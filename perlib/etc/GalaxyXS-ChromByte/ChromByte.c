/*
 * This file was generated automatically by xsubpp version 1.9508 from the
 * contents of ChromByte.xs. Do not edit this file, edit ChromByte.xs instead.
 *
 *	ANY CHANGES MADE HERE WILL BE LOST!
 *
 */

#line 1 "ChromByte.xs"
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"


#line 18 "ChromByte.c"
XS(XS_GalaxyXS__ChromByte_initchr); /* prototype to pass -Wmissing-prototypes */
XS(XS_GalaxyXS__ChromByte_initchr)
{
    dXSARGS;
    if (items != 1)
	Perl_croak(aTHX_ "Usage: GalaxyXS::ChromByte::initchr(len)");
    {
	int	len = (int)SvIV(ST(0));
	long	RETVAL;
	dXSTARG;
#line 14 "ChromByte.xs"
		void *address = malloc( len + 1 );
		memset( address , 0 , len + 1 );
		RETVAL = (long)address;
#line 33 "ChromByte.c"
	XSprePUSH; PUSHi((IV)RETVAL);
    }
    XSRETURN(1);
}

XS(XS_GalaxyXS__ChromByte_setbases); /* prototype to pass -Wmissing-prototypes */
XS(XS_GalaxyXS__ChromByte_setbases)
{
    dXSARGS;
    if (items != 4)
	Perl_croak(aTHX_ "Usage: GalaxyXS::ChromByte::setbases(address, begin, end, val)");
    {
	long	address = (long)SvIV(ST(0));
	int	begin = (int)SvIV(ST(1));
	int	end = (int)SvIV(ST(2));
	unsigned int	val = (unsigned int)SvUV(ST(3));
#line 27 "ChromByte.xs"
		char * buf = ( char * ) address ;
		memset( buf + begin , val , end - begin + 1 );
#line 53 "ChromByte.c"
    }
    XSRETURN_EMPTY;
}

XS(XS_GalaxyXS__ChromByte_getbase); /* prototype to pass -Wmissing-prototypes */
XS(XS_GalaxyXS__ChromByte_getbase)
{
    dXSARGS;
    if (items != 2)
	Perl_croak(aTHX_ "Usage: GalaxyXS::ChromByte::getbase(address, pos)");
    {
	long	address = (long)SvIV(ST(0));
	int	pos = (int)SvIV(ST(1));
	unsigned int	RETVAL;
	dXSTARG;
#line 35 "ChromByte.xs"
		char * buf = ( char * ) address ;
		RETVAL = *( buf + pos );
#line 72 "ChromByte.c"
	XSprePUSH; PUSHu((UV)RETVAL);
    }
    XSRETURN(1);
}

XS(XS_GalaxyXS__ChromByte_setbasesc); /* prototype to pass -Wmissing-prototypes */
XS(XS_GalaxyXS__ChromByte_setbasesc)
{
    dXSARGS;
    if (items != 3)
	Perl_croak(aTHX_ "Usage: GalaxyXS::ChromByte::setbasesc(address, pos, val)");
    {
	long	address = (long)SvIV(ST(0));
	int	pos = (int)SvIV(ST(1));
	char	val = (char)*SvPV_nolen(ST(2));
#line 46 "ChromByte.xs"
		char * buf = ( char * ) address ;
		memset( buf + pos , val , 1 );
#line 91 "ChromByte.c"
    }
    XSRETURN_EMPTY;
}

XS(XS_GalaxyXS__ChromByte_getbasec); /* prototype to pass -Wmissing-prototypes */
XS(XS_GalaxyXS__ChromByte_getbasec)
{
    dXSARGS;
    if (items != 2)
	Perl_croak(aTHX_ "Usage: GalaxyXS::ChromByte::getbasec(address, pos)");
    {
	long	address = (long)SvIV(ST(0));
	int	pos = (int)SvIV(ST(1));
	char	RETVAL;
	dXSTARG;
#line 54 "ChromByte.xs"
		char * buf = ( char * ) address ;
		RETVAL = *( buf + pos );
#line 110 "ChromByte.c"
	XSprePUSH; PUSHp((char *)&RETVAL, 1);
    }
    XSRETURN(1);
}

XS(XS_GalaxyXS__ChromByte_freechr); /* prototype to pass -Wmissing-prototypes */
XS(XS_GalaxyXS__ChromByte_freechr)
{
    dXSARGS;
    if (items != 1)
	Perl_croak(aTHX_ "Usage: GalaxyXS::ChromByte::freechr(address)");
    {
	long	address = (long)SvIV(ST(0));
#line 63 "ChromByte.xs"
		void * buf = ( void * ) address ;
		free( buf );
#line 127 "ChromByte.c"
    }
    XSRETURN_EMPTY;
}

#ifdef __cplusplus
extern "C"
#endif
XS(boot_GalaxyXS__ChromByte); /* prototype to pass -Wmissing-prototypes */
XS(boot_GalaxyXS__ChromByte)
{
    dXSARGS;
    char* file = __FILE__;

    XS_VERSION_BOOTCHECK ;

        newXS("GalaxyXS::ChromByte::initchr", XS_GalaxyXS__ChromByte_initchr, file);
        newXS("GalaxyXS::ChromByte::setbases", XS_GalaxyXS__ChromByte_setbases, file);
        newXS("GalaxyXS::ChromByte::getbase", XS_GalaxyXS__ChromByte_getbase, file);
        newXS("GalaxyXS::ChromByte::setbasesc", XS_GalaxyXS__ChromByte_setbasesc, file);
        newXS("GalaxyXS::ChromByte::getbasec", XS_GalaxyXS__ChromByte_getbasec, file);
        newXS("GalaxyXS::ChromByte::freechr", XS_GalaxyXS__ChromByte_freechr, file);
    XSRETURN_YES;
}

