/*
 * 1st ICEO test function
 *
 * Problem 1: The Sphere model
 *
 */



double
f_sphere(x, n)
register double		x[];
register int		n;


{
 
    	register int 	i;
    	double 		Sum;

    	for (Sum = 0.0, i = 0; i < n; i++) {
        	Sum += (x[i]-1.0)*(x[i]-1.0);	
    	}

    	return (Sum);
}


