

#include <math.h>


/*
 * 1st ICEO test function
 *
 * Problem 5: The Langerman's function 
 *
 */

double min = 500;

double center_point[] = {1., 1.3, .8, -.4, -1.3, 1.6, -.2, -.6, .5, 1.4, 
                        1., 1.3, .8, -.4, -1.3, 1.6, -.2, -.6, .5, 1.4};

double SqrDst( double *x1, double *x2, int n )
{
double dist = 0.0, d;
int i;

dist = 0;

for ( i = 0; i<n; i++ )
  {
    d =x1[i] - x2[i];
    dist += d*d;
  }

return (sqrt(dist));
}

double MaxSqrDst( double *x1, double *x2, int n )
{
  double dist = 0.0;
  int i;

  dist = 0.0;

  for(i=0; i<n; i++)
    dist = dist>fabs(x1[i] - x2[i])?dist:fabs(x1[i] - x2[i]);

  return (sqrt(n*dist*dist));
}


double
vicious_square(double *x, int n) 

{
 
        register int    i;

        double          Att,
                        PII,
                        dist;


        PII  = 3.141592653;

        
	dist = MaxSqrDst( x, center_point, n );

	Att = (exp(-dist/(2*PII)) * cos( PII * dist ) ) 
	 * ( 1 + 0.2 * (SqrDst( x, center_point, n )/(dist+.01)));

        return (-Att);
}




