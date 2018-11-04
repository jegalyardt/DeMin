/*
 * 1st ICEO test function
 *
 * Problem 4: The Michalewitz's function 
 *
 */

#include <math.h>

#define m  10.0

double Micha(x, n) /* Michalewitz */
double x[]; 
int n;
{ 
double	PI=3.1415927;
double   u;
int     i;

	u=0; 

	for (i=0;i<n;i++) 
		u = u + sin(x[i])
		  * pow(sin((i+1)*x[i]*x[i]/PI),2.0*m); 


	return(-u); 
}

