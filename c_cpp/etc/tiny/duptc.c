#include <stdlib.h> //EXIT_FAILURE
#include <stdio.h>
#include <float.h>
#include <math.h>

//#define SQRT2 1.4142135623730950488016887242097
long double Sqrt2=sqrtl(2);

double InsMean, InsSTD, ReadsDep;

/* http://en.wikipedia.org/wiki/Normal_distribution#Numerical_approximations_for_the_normal_CDF
如果用Ncdf(x+1)-Ncdf(x)，则当std<1时，最终结果会减半。
还原回Ncdf(x+0.5)-Ncdf(x-0.5)，则如同预期地趋近SE的结果。
 */
static inline long double NormalCDFdiff(double x, double m, double s) {
	//return 0.5 + 0.5 * erfl( (x - m) / (s * Sqrt2) );
	//return 0.5 * ( erfl( (x+1 - m) / (s * Sqrt2) ) - erfl( (x - m) / (s * Sqrt2) ) );
	return 0.5 * ( erfcl( (m-x-0.5) / (s * Sqrt2) ) - erfcl( (m-x+0.5) / (s * Sqrt2) ) );
}

int main (int argc, char *argv[]) {
	if (argc!=4) {
		fprintf(stderr,"Usage: %s <Reads_Depth> <InsertSize_Mean> <InsertSize_STD>\n",argv[0]);
		exit(EXIT_FAILURE);
	}
	ReadsDep=atof(argv[1]);
	InsMean=atof(argv[2]);
	InsSTD=atof(argv[3]);
	if (ReadsDep<0 || InsMean<0 || InsSTD<=0) {
		fputs("[x]Input Error !\n",stderr);
		exit(EXIT_FAILURE);
	}
	double SEdup = 1-exp(-ReadsDep);
	printf("ReadsDep=(%.9g) InsMean=(%.12g) InsSTD=(%.12g)\nSE theory dup. rate: %.20e\n",
		ReadsDep,InsMean,InsSTD,SEdup);
	long double sumdup=0;
	long double thisdup,Ci;
	long long int ins=1;
	do {
		Ci = (long double)ReadsDep * NormalCDFdiff(ins,InsMean,InsSTD);
		thisdup=Ci*(1-exp(-Ci));
		sumdup += thisdup;
		//printf("%d %.16Le %.16Le\n",ins,thisdup,sumdup);
		++ins;
	} while (ins<=2*InsMean+6*InsSTD || thisdup>LDBL_MIN);
	sumdup /= (long double)ReadsDep;
	printf("PE theory dup. rate: %.20Le (%.6Lf %% of SE)\n",sumdup,100*sumdup/SEdup);
	exit(EXIT_SUCCESS);
}

