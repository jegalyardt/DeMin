#include <math.h>

#define m  10.0
#define MaxDim 10

double Micha(y, n) /* Michalewitz */
double y[]; 
int n;
{ 
double x[MaxDim];
double	PI=3.1415927;
double   u, temp;
double cost, sint;
int     i;

	u=0; 
	
	cost = cos(PI/6);
	sint = sin(PI/6);

	for (i = 0; i <n-1; i+=2 )
	  {
	    x[i] = y[i]*cost - y[i+1]*sint;
	    x[i+1] = y[i]*sint + y[i+1]*cost;
	  }
	
	if ( i == n-1 ) x[i] = y[i];

	for (i=0;i<n;i++) 
		u = u + sin(x[i])
		  * pow(sin((i+1)*x[i]*x[i]/PI),2.0*m); 


	return(-u); 
}

