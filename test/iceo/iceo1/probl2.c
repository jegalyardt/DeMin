/*
 * 1st ICEO test function
 *
 * Problem 2: The Griewank's function
 *
 */




#include <math.h>

#define D	4000.0


double
f_Griew(x, n)
register double 	x[];  
register int		n;    /* The number of Dimensions */

{
 
    	register int 	i;

    	double 		Val1,
			Val2,
			Sum,
			sqrt(),
			cos();  

	for (Val1 = 0.0, Val2 = 1.0, i = 0; i < n; i++) {
		Val1 += (x[i]-100.0) * (x[i]-100.0);
		Val2 *= cos((x[i]-100.0) / sqrt((double) (i + 1)) );
	}

	Sum = Val1 / D - Val2 + 1.0;
    	return (Sum);
}







