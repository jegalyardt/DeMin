#include <math.h>


/* Returns the value of the function in itself */

Bump(x, obj, NVARS)
 	double x[], *obj;
        int NVARS;
{
	int ii;
	double sumc4, prodc2, sumsq;
 
	sumc4 = 0.0;
	prodc2 = 1.0;
	sumsq = 0.0;
	for (ii = 0; ii<NVARS; ii++) {
		sumc4 = sumc4 + cos(x[ii]) * cos(x[ii]) * cos(x[ii]) * cos(x[ii]);
		prodc2 = prodc2 * cos(x[ii]) * cos(x[ii]);
		sumsq = sumsq + (double)(ii+1) * x[ii] * x[ii];
	}
	sumsq = sqrt(sumsq);
 	*obj = (sumc4 - 2.0 * prodc2)/sumsq;
	if(*obj<0.0) *obj = - *obj;
	return;
}
 
/* Returns the values of the two contraints */

Constaints(x, con1, con2, NVARS)
	double x[], *con1, *con2;
        int NVARS;
{
	int ii;
 
	*con1 = 1.0;
	*con2 = 0.0;
	for (ii = 0; ii<NVARS; ii++) {
		*con1 = *con1 * x[ii];
		*con2 = *con2 + x[ii];
	}
	return;
}

