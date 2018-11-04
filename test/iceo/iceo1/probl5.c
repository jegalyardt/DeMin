#include <math.h>
#include "langerman.h"

/*
 * 1st ICEO test function
 *
 * Problem 5: The Langerman's function 
 *
 */

double min = 500;


double SqrDst( x1, x2, n )
double x1[], x2[];
int n;
{
double dist = 0.0, d;
int i;

dist = 0;

for ( i = 0; i<n; i++ )
  {
    d =x1[i] - x2[i];
    dist += d*d;
  }

return (dist);
}


double
langerman(x, n)  /* Langerman's function */
register double 	x[];
register int		n;

{
 
    	register int 	i;

    	double 		Sum,
	                PI,
	                dist;


	PI  = 3.141592653;
	
	Sum = 0;
	for ( i = 0; i < 5; i++ )
	  {

	    dist = SqrDst( x, a[i], n );
	    Sum -= c[i] * (exp(-dist/PI) * cos( PI * dist ) );
	  }

    	return (Sum);
}


