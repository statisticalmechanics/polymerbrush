/********************************************
* Box-Muller Gaussian random number generator
*********************************************/

#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

double BoxMuller(double mm, double ss)	/* normal random variate generator */
{				        /* mean mm, standard deviation ss */
	double xx1, xx2, ww, yy1;
	static double yy2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		yy1 = yy2;
		use_last = 0;
	}
	else
	{
		do {
			xx1 = 2.0 * RandomNumber() - 1.0;
			xx2 = 2.0 * RandomNumber() - 1.0;
			ww = xx1 * xx1 + xx2 * xx2;
		} while ( ww >= 1.0 );

		ww = sqrt( (-2.0 * log(ww) ) / ww );
		yy1 = xx1 * ww;
		yy2 = xx2 * ww;
		use_last = 1;
	}

	return( mm + yy1 * ss );
}


