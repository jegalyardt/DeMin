
double temp1(double d1, double d0)
{
   return (d1 - d0*d0);
}

double temp2(double d0)
{
     return (1. - d0);
}

double RozenBrock(double *x, int N)
{	
   double t0, tt, t1, d=0;
   int i;

   t0=x[0];
   for ( i=1; i < N; i++) 
   {
      t1 = x[i];
      tt = temp2(t0);
      d += tt*tt;
      tt = temp1(t1,t0);
      d += 100*tt*tt;
      
      t0 = t1;
   }
   
   return(d);
}


