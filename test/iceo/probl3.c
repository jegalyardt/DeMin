/*
 * 1st ICEO test function
 *
 * Problem 3: Shekel's Foxholes
 *
 */



#include "shekel.h"




double shekel(x, n)
double *x;
int n;
{
	register	i, j;
	double		sp, h, result = 0.0;

	for (i = 0; i < 30; i++) {
		sp = 0.0;
		for (j = 0; j < n; j++) {
			h   = x[j] - a[i][j];
			sp += h * h;
		}
		result += 1.0 / (sp + c[i]);
	}
	return( - result);
}






